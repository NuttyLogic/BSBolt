import datetime
import os
import subprocess
import time
import pysam
from BSB.BSB_Align.Align import BisulfiteAlignmentAndProcessing
from BSB.BSB_CallMethylation.ProcessContigs import ProcessContigs
from BSB.BSB_Index.RRBSGenomeBuild import RRBSGenomeIndexBuild
from BSB.BSB_Index.WholeGenomeBuild import WholeGenomeIndexBuild
from BSB.BSB_Matrix.MatrixAggregator import AggregateMatrix
from BSB.BSB_Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSB.BSB_Utils.UtilityFunctions import check_bowtie2_path


def launch_index(arguments):
    check_bowtie2_path(bowtie2_path=arguments.BT2)
    if arguments.rrbs:
        print(f'Generating RRBS Database at {arguments.DB}: '
              f'lower bound {arguments.rrbs_lower}, upper bound {arguments.rrbs_upper}: '
              f'Cut Format {arguments.rrbs_cut_format}')
        index = RRBSGenomeIndexBuild(reference_file=arguments.G,
                                     genome_database=arguments.DB,
                                     bowtie2_path=arguments.BT2,
                                     bowtie2_threads=arguments.BT2_p,
                                     cut_format=arguments.rrbs_cut_format,
                                     lower_bound=arguments.rrbs_lower,
                                     upper_bound=arguments.rrbs_upper)
        index.generate_rrbs_database()
    else:
        print(f'Generating WGBS Database at {arguments.DB}')
        index = WholeGenomeIndexBuild(reference_file=arguments.G,
                                      genome_database=arguments.DB,
                                      bowtie2_path=arguments.BT2,
                                      bowtie2_threads=arguments.BT2_p,
                                      mappable_regions=arguments.MR)
        index.generate_bsb_database()


def align_bisulfite(alignment_kwargs):
    start = time.time()
    print(f'Aligning {alignment_kwargs["fastq1"]} {alignment_kwargs["fastq2"]}')
    bs_alignment = BisulfiteAlignmentAndProcessing(**alignment_kwargs)
    bs_alignment.launch_bisulfite_aligment()
    alignment_time = datetime.timedelta(seconds=round(time.time() - start))
    print(f'Alignment Complete: Time {alignment_time}')
    start = time.time()
    print('Processing Reads')
    bs_alignment.process_reads()
    processing_time = datetime.timedelta(seconds=round(time.time() - start))
    print(f'Read Processing Complete: Time {processing_time}')
    print('Cleaning Temporary Files')
    bs_alignment.clean_temp_files()
    print(f'Total Run Time: {alignment_time + processing_time}')
    print('------------------------------')
    mapping_stats = process_mapping_statistics(bs_alignment.mapping_statistics)
    print(mapping_stats)


def process_mapping_statistics(mapping_dict):
    processed_list = [f'Total Reads: {mapping_dict["total_reads"]}',
                      f'Multi-mapped Reads: {mapping_dict["multimapped_reads"]}',
                      f'Multi-reference Reads: {mapping_dict["multireference_reads"]}',
                      f'Unmapped Reads: {mapping_dict["unmapped_reads"]}']
    mapped_reads = mapping_dict["total_reads"] - mapping_dict["multimapped_reads"] - mapping_dict["unmapped_reads"]
    mappability = mapped_reads / mapping_dict["total_reads"]
    processed_list.append(f'Mappability: {mappability:.3f}')
    processed_list.append('------------------------------')
    processed_list.append(f'Reads Mapped to Watson_C2T: {mapping_dict.get("W_C2T", 0)}')
    processed_list.append(f'Reads Mapped to Crick_C2T: {mapping_dict.get("C_C2T", 0)}')
    if 'W_G2A' in mapping_dict:
        processed_list.append(f'Reads Mapped to Watson_G2A: {mapping_dict.get("W_G2A", 0)}')
        processed_list.append(f'Reads Mapped to Crick_G2A: {mapping_dict.get("C_G2A", 0)}')
    return '\n'.join(processed_list)


def launch_alignment(arguments):
    check_bowtie2_path(bowtie2_path=arguments.BT2)
    bsb_command_dict = {arg[0]: str(arg[1]) for arg in arguments._get_kwargs()}
    arg_order = ['F1', 'F2', 'U', 'BT2', 'NC', 'O', 'DB', 'CP', 'CT', 'M', 'BT2_D', 'BT2_I', 'BT2_L',
                 'BT2_X', 'BT2_k', 'BT2_local', 'BT2_p', 'BT2_score_min']
    bowtie2_commands = ['--quiet', '--norc', '--no-mixed', '--no-discordant', '--sam-nohead', '--reorder',
                        '-k', str(arguments.BT2_k),
                        '-p', str(arguments.BT2_p),
                        '-L', str(arguments.BT2_L),
                        '-D', str(arguments.BT2_D)]
    if arguments.BT2_local:
        bowtie2_commands.append('--local')
        if arguments.BT2_score_min == 'L,-0.6,-0.6':
            arguments.BT2_score_min = 'G,20,8'
    else:
        bowtie2_commands.append('--end-to-end')
    bowtie2_commands.extend(['--score-min', arguments.BT2_score_min])
    if arguments.F2:
        pe_commands = ['-I', str(arguments.BT2_I),
                       '-X', str(arguments.BT2_X)]
        bowtie2_commands.extend(pe_commands)
    if not arguments.DB.endswith('/'):
        arguments.DB = f'{arguments.DB}/'
    command_line_arg = 'BSBolt Align ' + ' '.join([f'-{arg} {bsb_command_dict[arg]}' for arg in arg_order])
    aligment_kwargs = dict(fastq1=arguments.F1, fastq2=arguments.F2, undirectional_library=arguments.U,
                           bowtie2_commands=bowtie2_commands, bsb_database=arguments.DB,
                           bowtie2_path=arguments.BT2, output_path=arguments.O,
                           conversion_threshold=(arguments.CP, arguments.CT), mismatch_threshold=arguments.M,
                           command_line_arg=command_line_arg, non_converted_output=arguments.NC,
                           output_unmapped=arguments.OU)
    align_bisulfite(aligment_kwargs)
    if arguments.S:
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
                                      remove_sx_reads=arguments.remove_sx,
                                      text_output=arguments.text,
                                      min_read_depth=arguments.min,
                                      threads=arguments.t,
                                      verbose=arguments.verbose,
                                      min_base_quality=arguments.min_qual,
                                      cg_only=arguments.CG,
                                      ATCGmap=arguments.ATCG)
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
                             verbose=arguments.verbose)
    matrix.aggregate_matrix()


def launch_simulation(arguments):
    if os.path.isfile(arguments.A):
        read_simulation = SimulateMethylatedReads(reference_file=arguments.G, art_path=arguments.A,
                                                  output_path=arguments.O, paired_end=arguments.PE,
                                                  read_length=arguments.RL, read_depth=arguments.RD,
                                                  undirectional=arguments.U, methylation_reference_output=arguments.RO,
                                                  methylation_reference=arguments.BR, methylation_profile=arguments.RC)
        read_simulation.run_simulation()
    else:
        print('ART Executable Path not Valid')


bsb_launch = {'Index': launch_index,
              'Align': launch_alignment,
              'CallMethylation': launch_methylation_call,
              'AggregateMatrix': launch_matrix_aggregation,
              'Simulate': launch_simulation}
