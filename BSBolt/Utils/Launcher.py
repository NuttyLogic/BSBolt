import datetime
import subprocess
import time
import pysam
from BSBolt.Align.AlignReads import BisulfiteAlignmentAndProcessing
from BSBolt.CallMethylation.ProcessContigs import ProcessContigs
from BSBolt.Impute.kNN_Impute import ImputeMissingValues
from BSBolt.Index.RRBSGenomeBuild import RRBSGenomeIndexBuild
from BSBolt.Index.WholeGenomeBuild import WholeGenomeIndexBuild
from BSBolt.Matrix.MatrixAggregator import AggregateMatrix
from BSBolt.Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSBolt.Utils.UtilityFunctions import get_external_paths

bwa_path, art_path = get_external_paths()


def launch_index(arguments):
    if arguments.rrbs:
        print(f'Generating RRBS Database at {arguments.DB}: '
              f'lower bound {arguments.rrbs_lower}, upper bound {arguments.rrbs_upper}: '
              f'Cut Format {arguments.rrbs_cut_format}')
        index = RRBSGenomeIndexBuild(reference_file=arguments.G,
                                     genome_database=arguments.DB,
                                     bwa_path=bwa_path,
                                     cut_format=arguments.rrbs_cut_format,
                                     lower_bound=arguments.rrbs_lower,
                                     upper_bound=arguments.rrbs_upper)
        index.generate_rrbs_database()
    else:
        print(f'Generating WGBS Database at {arguments.DB}')
        index = WholeGenomeIndexBuild(reference_file=arguments.G,
                                      genome_database=arguments.DB,
                                      bwa_path=bwa_path,
                                      mappable_regions=arguments.MR)
        index.generate_bsb_database()


def align_bisulfite(bwa_cmd, output_path):
    start = time.time()
    print(' '.join(bwa_cmd))
    bs_alignment = BisulfiteAlignmentAndProcessing(bwa_cmd, output_path)
    bs_alignment.align_reads()
    alignment_time = datetime.timedelta(seconds=round(time.time() - start))
    print(f'Alignment Complete: Time {alignment_time}')
    print('------------------------------')
    mapping_stats = process_mapping_statistics(bs_alignment.mapping_statistics)
    print(mapping_stats)


def process_mapping_statistics(mapping_dict):
    dict(TotalReads=0, TotalAlignments=0, BSAmbiguous=0, C_C2T=0, C_G2A=0,
         W_C2T=0, W_G2A=0, Unaligned=0)
    processed_list = []
    mappability = (mapping_dict['TotalAlignments'] - mapping_dict['Unaligned']) / mapping_dict['TotalAlignments']
    processed_list.append(f'Total Reads: {mapping_dict["TotalReads"]}')
    processed_list.append(f'Mappability: {mappability * 100:.3f} %')
    processed_list.append('------------------------------')
    processed_list.append(f'Reads Mapped to Watson_C2T: {mapping_dict["W_C2T"]}')
    processed_list.append(f'Reads Mapped to Crick_C2T: {mapping_dict["C_C2T"]}')
    processed_list.append(f'Reads Mapped to Watson_G2A: {mapping_dict["W_G2A"]}')
    processed_list.append(f'Reads Mapped to Crick_G2A: {mapping_dict["C_G2A"]}')
    processed_list.append('------------------------------')
    processed_list.append(f'Unmapped Reads (Single / Paired Ends): {mapping_dict["Unaligned"]}')
    processed_list.append(f'\tBisulfite Ambiguous: {mapping_dict["BSAmbiguous"]}')
    return '\n'.join(processed_list)


def launch_alignment(arguments):
    bsb_command_dict = {arg[0]: str(arg[1]) for arg in arguments._get_kwargs()}
    bwa_cmd = [bwa_path, 'mem', '-Y']
    if bsb_command_dict['UN'] == 'True':
        bwa_cmd.extend(['-z', 'true'])
    bool_args = ['M', 'S', 'j', 'p']
    for arg in bool_args:
        if bsb_command_dict[arg] == 'True':
            bwa_cmd.append(f'-{arg}')
    default_args = ['A', 'B', 'D', 'E', 'L', 'T', 'U', 'W', 'c', 'd', 'k', 'm', 'r', 't', 'w', 'y']
    for arg in default_args:
        bwa_cmd.extend([f'-{arg}', bsb_command_dict[arg]])
    if bsb_command_dict['H'] != 'None':
        bwa_cmd.extend([f'-H', bsb_command_dict['H']])
    if bsb_command_dict['I'] != 'None':
        bwa_cmd.extend([f'-I', bsb_command_dict['I']])
    if bsb_command_dict['INDEL']:
        bwa_cmd.extend([f'-O', bsb_command_dict['INDEL']])
    if bsb_command_dict['XA']:
        bwa_cmd.extend([f'-h', bsb_command_dict['XA']])
    database = bsb_command_dict['G']
    if not database.endswith('.fa'):
        if not database.endswith('/'):
            database = f'{database}/BSB_ref.fa'
        else:
            database = f'{database}BSB_ref.fa'
    bwa_cmd.append(database)
    bwa_cmd.append(bsb_command_dict['F1'])
    if bsb_command_dict['F2'] != 'None':
        bwa_cmd.append(bsb_command_dict['F2'])
    align_bisulfite(bwa_cmd, arguments.O)
    if arguments.Sort:
        pysam.sort('-o', f'{arguments.O}.sorted.bam', f'{arguments.O}.bam')
        subprocess.run(['rm', f'{arguments.O}.bam'])
        pysam.index(f'{arguments.O}.sorted.bam')


def launch_methylation_call(arguments):
    if arguments.CG and arguments.ATCG:
        assert False, 'Reporting only CG sites notes compatible with .ATCGmap output'
    methylation_call = ProcessContigs(input_file=arguments.I,
                                      genome_database=arguments.DB,
                                      output_prefix=arguments.O,
                                      remove_ccgg=arguments.remove_ccgg,
                                      text_output=arguments.text,
                                      min_read_depth=arguments.min,
                                      threads=arguments.t,
                                      verbose=arguments.verbose,
                                      min_base_quality=arguments.min_qual,
                                      cg_only=arguments.CG,
                                      ATCGmap=arguments.ATCG,
                                      ignore_orphans=arguments.IO)
    methylation_call.process_contigs()
    methylation_call.watch_pool()


def launch_matrix_aggregation(arguments):
    def get_sample_info(file_path):
        file_list = []
        with open(file_path, 'r') as file:
            for line in file:
                file_list.append(line.replace('\n', ''))
        return file_list

    if len(arguments.F) == 1:
        arguments.F = get_sample_info(arguments.F[0])
    if arguments.S:
        if len(arguments.S) == 1:
            arguments.S = get_sample_info(arguments.S[0])

    matrix = AggregateMatrix(file_list=arguments.F,
                             sample_list=arguments.S,
                             min_site_coverage=arguments.min_coverage,
                             site_proportion_threshold=arguments.min_sample,
                             output_path=arguments.O,
                             cg_only=arguments.CG,
                             verbose=arguments.verbose,
                             threads=arguments.t,
                             count_matrix=arguments.count)
    matrix.aggregate_matrix()


def launch_simulation(arguments):
    read_simulation = SimulateMethylatedReads(reference_file=arguments.G, art_path=art_path,
                                              output_path=arguments.O, paired_end=arguments.PE,
                                              read_length=arguments.RL, read_depth=arguments.RD,
                                              undirectional=arguments.U, methylation_reference_output=arguments.RO,
                                              methylation_reference=arguments.BR, methylation_profile=arguments.RC,
                                              insertion_rate1=arguments.IR1, insertion_rate2=arguments.IR2,
                                              deletion_rate1=arguments.DR1, deletion_rate2=arguments.DR2,
                                              n_base_cutoff=arguments.NF, sequencing_system=arguments.SS,
                                              pe_fragment_size=arguments.M, fragment_size_deviation=arguments.SM,
                                              read1_quality_profile=arguments.Q1,
                                              read2_quality_profile=arguments.Q2)
    read_simulation.run_simulation()


def launch_imputation(arguments):
    output_path = arguments.O
    if not output_path:
        output_path = f'{arguments.M}_imputed.txt'
    impute = ImputeMissingValues(input_matrix_file=arguments.M, batch_size=arguments.B,
                                 imputation_window_size=arguments.W, k=arguments.k, threads=arguments.t,
                                 verbose=arguments.verbose, output_path=output_path, randomize_batch=arguments.R)
    impute.import_matrix()
    impute.impute_values()
    impute.output_imputed_matrix()


bsb_launch = {'Index': launch_index,
              'Align': launch_alignment,
              'CallMethylation': launch_methylation_call,
              'AggregateMatrix': launch_matrix_aggregation,
              'Simulate': launch_simulation,
              'Impute': launch_imputation}
