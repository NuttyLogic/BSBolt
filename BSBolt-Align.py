import argparse
import datetime
import subprocess
import time
import pysam
from BSB_Align.BSB_Align import BisulfiteAlignmentAndProcessing
from BSB_Utils.BSB_UtilityFunctions import check_bowtie2_path


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
    for stat in mapping_stats:
        print(stat)


def process_mapping_statistics(mapping_dict):
    processed_list = [f'Total Reads: {mapping_dict["total_reads"]}',
                      f'Multi-mapped Reads: {mapping_dict["multimapped_reads"]}',
                      f'Unmapped Reads: {mapping_dict["unmapped_reads"]}']
    mapped_reads = mapping_dict["total_reads"] - mapping_dict["multimapped_reads"] - mapping_dict["unmapped_reads"]
    mappability = mapped_reads / mapping_dict["total_reads"]
    processed_list.append(f'Mappability: {mappability:.3f}')
    processed_list.append('------------------------------')
    processed_list.append(f'Reads Mapped to Watson_C2T: {mapping_dict["W_C2T"]}')
    processed_list.append(f'Reads Mapped to Crick_C2T: {mapping_dict["C_C2T"]}')
    if 'W_G2A' in mapping_dict:
        processed_list.append(f'Reads Mapped to Watson_G2A: {mapping_dict["W_G2A"]}')
        processed_list.append(f'Reads Mapped to Crick_G2A: {mapping_dict["C_G2A"]}')
    return processed_list


parser = argparse.ArgumentParser(description='BS-Bolt Sequencing Alignment Tool')

parser.add_argument('-F1', type=str, default=None, help='Path to fastq 1', required=True)
parser.add_argument('-F2', type=str, default=None, help='Path to fastq 2')
parser.add_argument('-U', action="store_false", default=True, help='Library undirectioinal, default=True')
parser.add_argument('-BT2', type=str, default='bowtie2', help='Path to bowtie2 aligner')
parser.add_argument('-O', type=str, default=None, help='Path to Output Prefix', required=True)
parser.add_argument('-DB', type=str, default=None, help='Path to BSSeeker Database', required=True)
parser.add_argument('-CP', type=float, default=0.5, help='Proportion threshold to label read not fully converted')
parser.add_argument('-CT', type=int, default=5, help='Number of mCH that must be observed to label a read unconverted')
parser.add_argument('-M', type=int, default=4, help='Read mismatch threshold, reads with mismatches greater than '
                                                    'threshold will be discarded')
parser.add_argument('-S', action="store_true", default=False, help='Position Sort Output Bam, default=False')
parser.add_argument('-BT2-local', action="store_true", default=False,
                    help='Bowtie2 local alignment or end-to-end, default end-to-end')
parser.add_argument('-BT2-D', type=int, default=40, help='Bowtie2 number of consecutive seed extension attempts '
                                                         'that can fail before Bowtie2 move on')
parser.add_argument('-BT2-k', type=int, default=2, help='Bowtie2 alignment search limit')
parser.add_argument('-BT2-p', type=int, default=2, help='Number of threads for Bowtie2 to use')
parser.add_argument('-BT2-L', type=int, default=20, help='Length of subseeds during alignment')
parser.add_argument('-BT2-score-min', type=str, default='L,-0.6,-0.6', help='Bowtie2 scoring function')
parser.add_argument('-BT2-I', type=int, default=0, help='Bowtie2, minimum fragment length for a valid paired-end '
                                                        'alignment')
parser.add_argument('-BT2-X', type=int, default=500, help='Bowtie2, maximum fragment length for a valid paired-end '
                                                          'alignment')


if __name__ == "__main__":
    arguments = parser.parse_args()
    bsseeker_command_dict = {arg[0]: str(arg[1]) for arg in arguments._get_kwargs()}
    arg_order = ['F1', 'F2', 'U', 'BT2', 'O', 'DB', 'CP', 'CT', 'M',  'BT2_D', 'BT2_I', 'BT2_L',
                 'BT2_X', 'BT2_k', 'BT2_local', 'BT2_p', 'BT2_score_min']
    bowtie2_commands = ['--quiet', '--norc', '--no-mixed', '--no-discordant', '--sam-nohead', '--reorder',
                        '-k', str(arguments.BT2_k),
                        '-p', str(arguments.BT2_p),
                        '-L', str(arguments.BT2_L),
                        '-D', str(arguments.BT2_D)]
    check_bowtie2_path(bowtie2_path=arguments.BT2)
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
    command_line_arg = 'BSBolt-Align.py ' + ' '.join([f'-{arg} {bsseeker_command_dict[arg]}' for arg in arg_order])
    aligment_kwargs = dict(fastq1=arguments.F1, fastq2=arguments.F2, undirectional_library=arguments.U,
                           bowtie2_commands=bowtie2_commands, bsseeker_database=arguments.DB,
                           bowtie2_path=arguments.BT2, output_path=arguments.O,
                           conversion_threshold=(arguments.CP, arguments.CT), mismatch_threshold=arguments.M,
                           command_line_arg=command_line_arg)
    align_bisulfite(aligment_kwargs)
    if arguments.S:
        pysam.sort('-o', f'{arguments.O}.sorted.bam', f'{arguments.O}.bam')
        subprocess.run(['rm', f'{arguments.O}.bam'])
        pysam.index(f'{arguments.O}.sorted.bam')
    else:
        pysam.index(f'{arguments.O}.bam')
