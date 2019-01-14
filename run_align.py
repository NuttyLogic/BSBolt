from BSB_Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSB_Align.BSB_Align import BisulfiteAlignmentAndProcessing
import time
from BSB_Utils.BSB_UtilityFunctions import reverse_complement
import subprocess
from BSB_CallMethylation.ProcessContigs import ProcessContigs
from BSB_CallMethylation.CallMethylation import CallMethylation
bsseeker_base_path = '/Users/colinfarrell/Documents/'
#bsseeker_base_path = '/media/colin/Linux_Data/bsseeker_project/'
#bowtie2_path = '/media/colin/Linux_Data/bsseeker_project/bowtie2-2.3.4.3-linux-x86_64/bowtie2'
bowtie2_path = '/Users/colinfarrell/Desktop/bsseeker_refactor_test_folder/bowtie2-2.3.4.2-macos-x86_64/bowtie2'
#art_path = '/media/colin/Linux_Data/bsseeker_project/art_bin_MountRainier/art_illumina'
art_path = '/Users/colinfarrell/Desktop/bsseeker_refactor_test_folder/art_bin_MountRainier/art_illumina'
'''
test_simulation_files = SimulateMethylatedReads(reference_file=f'{bsb_directory}Tests/TestData/BSB_test.fa',
                               art_path=art_path, paired_end=True,
                               output_path=f'{bsb_directory}Tests/TestSimulations/bsseeker_pe',
                                                undirectional=True)

test_simulation_files.run_simulation()


bt2_command_list = ['--quiet', '--end-to-end', '-D', '50', '--norc', '--no-mixed', '--no-discordant', '-L', '15', '--sam-nohead', '-p', '2', '--score-min', 'L,-0.6,-0.6','-k', '2', '--reorder']
start = time.time()

read_processing_test = BisulfiteAlignmentAndProcessing(bowtie2_path=bowtie2_path,
                                    bsseeker_database=f'{bsb_directory}Tests/TestData/test_bsseeker_db/',
                                    fastq1=f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth_1.fastq',
                                    fastq2=f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth_2.fastq',
                                           bowtie2_commands=bt2_command_list, output_path=f'{bsb_directory}Tests/test_out_pe',
                                                        undirectional_library=False, mismatch_threshold=4)
read_processing_test.launch_bisulfite_aligment()
read_processing_test.process_reads()
print(time.time() - start)
print(read_processing_test.mapping_statistics)

bt2_command_list = ['--quiet', '--end-to-end', '-D', '50', '--norc', '--no-mixed', '--no-discordant', '-L', '15', '--sam-nohead', '-p', '2', '--score-min', 'L,-0.6,-0.6','-k', '2', '--reorder']
start = time.time()

read_processing_test_2 = BisulfiteAlignmentAndProcessing(bowtie2_path=bowtie2_path,
                                    bsseeker_database=f'{bsb_directory}Tests/TestData/test_bsseeker_db/',
                                    fastq1=f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth_1.fastq',
                                           bowtie2_commands=bt2_command_list, output_path=f'{bsb_directory}Tests/test_out_se',
                                                        undirectional_library=False, mismatch_threshold=4)
read_processing_test_2.launch_bisulfite_aligment()
read_processing_test_2.process_reads()
print(time.time() - start)
print(read_processing_test_2.mapping_statistics)

start = time.time()
bt2_command_list = ['--quiet', '--end-to-end', '-D', '50', '--norc', '--no-mixed', '--no-discordant', '-L', '15', '--sam-nohead', '-p', '2', '--score-min', 'L,-0.6,-0.6','-k', '2', '--reorder']

read_processing_test_3 = BisulfiteAlignmentAndProcessing(bowtie2_path=bowtie2_path,
                                    bsseeker_database=f'{bsb_directory}Tests/TestData/test_database_rrbs/',
                                    fastq1=f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth_1.fastq',
                                           bowtie2_commands=bt2_command_list, output_path=f'{bsb_directory}Tests/test_out_rrbs',
                                                        undirectional_library=True, mismatch_threshold=4)
read_processing_test_3.launch_bisulfite_aligment()
read_processing_test_3.process_reads()
print(time.time() - start)
print(read_processing_test_3.mapping_statistics)






bsseeker_2_output = {}

for line in open(f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth.sam', 'r'):
    line_split = line.replace('\n', '').split('\t')
    if '@' not in line_split[0]:
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, *SAM_TAGS = line_split
        bsseeker_2_output[QNAME] = dict(FLAG=FLAG, RNAME=RNAME, POS=POS, MAPQ=MAPQ, CIGAR=CIGAR, RNEXT=RNEXT, PNEXT=PNEXT,
                            TLEN=TLEN, SEQ=SEQ, QUAL=QUAL, SAM_TAGS=SAM_TAGS)

refactor_output = {}

line_count = 0

previous_read = None

for line in open(f'{bsb_directory}Tests/test_out_pe.sam', 'r'):
    line_split = line.replace('\n', '').split('\t')
    if '@' not in line_split[0]:
        line_count += 1
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, *SAM_TAGS = line_split
        if QNAME != previous_read:
            previous_read = QNAME
            refactor_output[f'{QNAME}/2'] = dict(FLAG=FLAG, RNAME=RNAME, POS=POS, MAPQ=MAPQ, CIGAR=CIGAR, RNEXT=RNEXT, PNEXT=PNEXT,
                            TLEN=TLEN, SEQ=SEQ, QUAL=QUAL, SAM_TAGS=SAM_TAGS)
        else:
            refactor_output[f'{QNAME}/1'] = dict(FLAG=FLAG, RNAME=RNAME, POS=POS, MAPQ=MAPQ, CIGAR=CIGAR, RNEXT=RNEXT,
                                                 PNEXT=PNEXT,
                                                 TLEN=TLEN, SEQ=SEQ, QUAL=QUAL, SAM_TAGS=SAM_TAGS)
            previous_read = None

print(len(bsseeker_2_output))
print(len(refactor_output))
missing_count = 0

strand_count = 0

conversion_dict = {}
position_dict = {}
unmatched_seq = {}

for key, value in bsseeker_2_output.items():
    try:
        refactor_value = refactor_output[key]
    except KeyError:
        missing_count += 1
        continue
    if value['SAM_TAGS'][1] != refactor_value['SAM_TAGS'][-2]:
        #print(key, value)
        #print(key, refactor_value)
        conversion_strand = f'Bsseeker {value["SAM_TAGS"][0]} : refactor {refactor_value["SAM_TAGS"][-1]}'
        try:
            conversion_dict[conversion_strand] += 1
        except KeyError:
            conversion_dict[conversion_strand] = 1
    if value['POS'] != refactor_value['POS'] and int(value['POS']) != int(refactor_value['PNEXT']):
        print('Pos Wrong', key, refactor_value["SAM_TAGS"][-1])
        print(value)
        print(refactor_value)
        print(int(value['POS']) - int(refactor_value['POS']), int(value['POS']) - int(refactor_value['PNEXT']))
        try:
            position_dict[f'{refactor_value["FLAG"]}_{refactor_value["SAM_TAGS"][-1]}'] += 1
        except KeyError:
            position_dict[f'{refactor_value["FLAG"]}_{refactor_value["SAM_TAGS"][-1]}'] = 1
    if value['SEQ'] != refactor_value['SEQ'] and value['POS'] == refactor_value['POS']:
        try:
            unmatched_seq[f'{refactor_value["FLAG"]}_{refactor_value["SAM_TAGS"][-1]}'] += 1
        except KeyError:
            unmatched_seq[f'{refactor_value["FLAG"]}_{refactor_value["SAM_TAGS"][-1]}'] = 1
        if reverse_complement(refactor_value['SEQ']) != value['SEQ']:
            print('Sequences Not Equivalent')
            print(value)
            print(refactor_value)


print(missing_count)
print(conversion_dict)
print(position_dict)
print(unmatched_seq)
print('----------------------------------MethylationCalling-------------------------------------')

input_sam = f'{bsb_directory}Tests/test_out_pe.sam'
sorted_bam = f'{bsb_directory}Tests/test_out_pe.sorted.bam'
sort_command = ['samtools', 'sort', input_sam, '-O', 'bam', '-o', sorted_bam]
subprocess.run(args=sort_command)
index_command = ['samtools', 'index', sorted_bam]
subprocess.run(args=index_command)



ProcessContigs(input_file=f'{bsb_directory}Tests/test_out_pe.sorted.bam',
               genome_database=f'{bsb_directory}Tests/TestData/test_bsseeker_db/',
               output_prefix=f'{bsb_directory}Tests/test_pe_methylation_call_2',
               remove_ccgg=False,
               ignore_overlap=False,
               remove_sx_reads=False,
               text_output=True,
               min_read_depth=5,
               threads=6)
'''
bs_align_args = ['python3', 'BSBolt-Align.py', '-F1', f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth_1.fastq',
                 '-F2', f'{bsb_directory}Tests/TestSimulations/bsseeker_pe_meth_2.fastq',
                 '-BT2', bowtie2_path, '-O', f'{bsb_directory}Tests/test_out_pe2', '-DB',
                 f'{bsb_directory}Tests/TestData/test_bsseeker_db/']
print(' '.join(bs_align_args))
subprocess.run(bs_align_args)

test_contigs = ['chr1', 'chr10', 'chr11', 'chr11_KI270721v1_random', 'chr12', 
'chr13', 'chr14', 
                'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 
                'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 
                'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15', 'chr15_KI270727v1_random', 
                'chr16', 'chr16_KI270728v1_random', 'chr17', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 
                'chr17_KI270730v1_random', 'chr18', 'chr19', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random',
                'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 
                'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2', 'chr20', 'chr21',
                'chr22', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 
                'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 
                'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 
                'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3', 'chr3_GL000221v1_random',
                'chr4', 'chr4_GL000008v2_random', 'chr5', 'chr5_GL000208v1_random', 'chr6', 'chr7', 'chr8', 'chr9', 
                'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 
                'chr1_KI270762v1_alt', 'chr1_KI270766v1_alt', 'chr1_KI270760v1_alt', 'chr1_KI270765v1_alt', 
                'chr1_GL383518v1_alt', 'chr1_GL383519v1_alt', 'chr1_GL383520v2_alt', 'chr1_KI270764v1_alt', 
                'chr1_KI270763v1_alt', 'chr1_KI270759v1_alt', 'chr1_KI270761v1_alt', 'chr2_KI270770v1_alt', 
                'chr2_KI270773v1_alt', 'chr2_KI270774v1_alt', 'chr2_KI270769v1_alt', 'chr2_GL383521v1_alt', 
                'chr2_KI270772v1_alt', 'chr2_KI270775v1_alt', 'chr2_KI270771v1_alt', 'chr2_KI270768v1_alt', 
                'chr2_GL582966v2_alt', 'chr2_GL383522v1_alt', 'chr2_KI270776v1_alt', 'chr2_KI270767v1_alt', 
                'chr3_JH636055v2_alt', 'chr3_KI270783v1_alt', 'chr3_KI270780v1_alt', 'chr3_GL383526v1_alt', 
                'chr3_KI270777v1_alt', 'chr3_KI270778v1_alt', 'chr3_KI270781v1_alt', 'chr3_KI270779v1_alt', 
                'chr3_KI270782v1_alt', 'chr3_KI270784v1_alt', 'chr4_KI270790v1_alt', 'chr4_GL383528v1_alt', 
                'chr4_KI270787v1_alt', 'chr4_GL000257v2_alt', 'chr4_KI270788v1_alt', 'chr4_GL383527v1_alt', 
                'chr4_KI270785v1_alt', 'chr4_KI270789v1_alt', 'chr4_KI270786v1_alt', 'chr5_KI270793v1_alt',
                'chr5_KI270792v1_alt', 'chr5_KI270791v1_alt', 'chr5_GL383532v1_alt', 'chr5_GL949742v1_alt', 
                'chr5_KI270794v1_alt', 'chr5_GL339449v2_alt', 'chr5_GL383530v1_alt', 'chr5_KI270796v1_alt', 
                'chr5_GL383531v1_alt', 'chr5_KI270795v1_alt', 'chr6_GL000250v2_alt', 'chr6_KI270800v1_alt', 
                'chr6_KI270799v1_alt', 'chr6_GL383533v1_alt', 'chr6_KI270801v1_alt', 'chr6_KI270802v1_alt', 
                'chr6_KB021644v2_alt', 'chr6_KI270797v1_alt', 'chr6_KI270798v1_alt', 'chr7_KI270804v1_alt', 
                'chr7_KI270809v1_alt', 'chr7_KI270806v1_alt', 'chr7_GL383534v2_alt', 'chr7_KI270803v1_alt', 
                'chr7_KI270808v1_alt', 'chr7_KI270807v1_alt', 'chr7_KI270805v1_alt', 'chr8_KI270818v1_alt', 
                'chr8_KI270812v1_alt', 'chr8_KI270811v1_alt', 'chr8_KI270821v1_alt', 'chr8_KI270813v1_alt', 'chr8_KI270822v1_alt', 'chr8_KI270814v1_alt', 'chr8_KI270810v1_alt', 'chr8_KI270819v1_alt', 'chr8_KI270820v1_alt', 'chr8_KI270817v1_alt', 'chr8_KI270816v1_alt', 'chr8_KI270815v1_alt', 'chr9_GL383539v1_alt', 'chr9_GL383540v1_alt', 'chr9_GL383541v1_alt', 'chr9_GL383542v1_alt', 'chr9_KI270823v1_alt', 'chr10_GL383545v1_alt', 'chr10_KI270824v1_alt', 'chr10_GL383546v1_alt', 'chr10_KI270825v1_alt', 'chr11_KI270832v1_alt', 'chr11_KI270830v1_alt', 'chr11_KI270831v1_alt', 'chr11_KI270829v1_alt', 'chr11_GL383547v1_alt', 'chr11_JH159136v1_alt', 'chr11_JH159137v1_alt', 'chr11_KI270827v1_alt', 'chr11_KI270826v1_alt', 'chr12_GL877875v1_alt', 'chr12_GL877876v1_alt', 'chr12_KI270837v1_alt', 'chr12_GL383549v1_alt', 'chr12_KI270835v1_alt', 'chr12_GL383550v2_alt', 'chr12_GL383552v1_alt', 'chr12_GL383553v2_alt', 'chr12_KI270834v1_alt', 'chr12_GL383551v1_alt', 'chr12_KI270833v1_alt', 'chr12_KI270836v1_alt', 'chr13_KI270840v1_alt', 'chr13_KI270839v1_alt', 'chr13_KI270843v1_alt', 'chr13_KI270841v1_alt', 'chr13_KI270838v1_alt', 'chr13_KI270842v1_alt', 'chr14_KI270844v1_alt', 'chr14_KI270847v1_alt', 'chr14_KI270845v1_alt', 'chr14_KI270846v1_alt', 'chr15_KI270852v1_alt', 'chr15_KI270851v1_alt', 'chr15_KI270848v1_alt', 'chr15_GL383554v1_alt', 'chr15_KI270849v1_alt', 'chr15_GL383555v2_alt', 'chr15_KI270850v1_alt', 'chr16_KI270854v1_alt', 'chr16_KI270856v1_alt', 'chr16_KI270855v1_alt', 'chr16_KI270853v1_alt', 'chr16_GL383556v1_alt', 'chr16_GL383557v1_alt', 'chr17_GL383563v3_alt', 'chr17_KI270862v1_alt', 'chr17_KI270861v1_alt', 'chr17_KI270857v1_alt', 'chr17_JH159146v1_alt', 'chr17_JH159147v1_alt', 'chr17_GL383564v2_alt', 'chr17_GL000258v2_alt', 'chr17_GL383565v1_alt', 'chr17_KI270858v1_alt', 'chr17_KI270859v1_alt', 'chr17_GL383566v1_alt', 'chr17_KI270860v1_alt', 'chr18_KI270864v1_alt', 'chr18_GL383567v1_alt', 'chr18_GL383570v1_alt', 'chr18_GL383571v1_alt', 'chr18_GL383568v1_alt', 'chr18_GL383569v1_alt', 'chr18_GL383572v1_alt', 'chr18_KI270863v1_alt', 'chr19_KI270868v1_alt', 'chr19_KI270865v1_alt', 'chr19_GL383573v1_alt', 'chr19_GL383575v2_alt', 'chr19_GL383576v1_alt', 'chr19_GL383574v1_alt', 'chr19_KI270866v1_alt', 'chr19_KI270867v1_alt', 'chr19_GL949746v1_alt', 'chr20_GL383577v2_alt', 'chr20_KI270869v1_alt', 'chr20_KI270871v1_alt', 'chr20_KI270870v1_alt', 'chr21_GL383578v2_alt', 'chr21_KI270874v1_alt', 'chr21_KI270873v1_alt', 'chr21_GL383579v2_alt', 'chr21_GL383580v2_alt', 'chr21_GL383581v2_alt', 'chr21_KI270872v1_alt', 'chr22_KI270875v1_alt', 'chr22_KI270878v1_alt', 'chr22_KI270879v1_alt', 'chr22_KI270876v1_alt', 'chr22_KI270877v1_alt', 'chr22_GL383583v2_alt', 'chr22_GL383582v2_alt', 'chrX_KI270880v1_alt', 'chrX_KI270881v1_alt', 'chr19_KI270882v1_alt', 'chr19_KI270883v1_alt', 'chr19_KI270884v1_alt', 'chr19_KI270885v1_alt', 'chr19_KI270886v1_alt', 'chr19_KI270887v1_alt', 'chr19_KI270888v1_alt', 'chr19_KI270889v1_alt', 'chr19_KI270890v1_alt', 'chr19_KI270891v1_alt', 'chr1_KI270892v1_alt', 'chr2_KI270894v1_alt', 'chr2_KI270893v1_alt', 'chr3_KI270895v1_alt', 'chr4_KI270896v1_alt', 'chr5_KI270897v1_alt', 'chr5_KI270898v1_alt', 'chr6_GL000251v2_alt', 'chr7_KI270899v1_alt', 'chr8_KI270901v1_alt', 'chr8_KI270900v1_alt', 'chr11_KI270902v1_alt', 'chr11_KI270903v1_alt', 'chr12_KI270904v1_alt', 'chr15_KI270906v1_alt', 'chr15_KI270905v1_alt', 'chr17_KI270907v1_alt', 'chr17_KI270910v1_alt', 'chr17_KI270909v1_alt', 'chr17_JH159148v1_alt', 'chr17_KI270908v1_alt', 'chr18_KI270912v1_alt', 'chr18_KI270911v1_alt', 'chr19_GL949747v2_alt', 'chr22_KB663609v1_alt', 'chrX_KI270913v1_alt', 'chr19_KI270914v1_alt', 'chr19_KI270915v1_alt', 'chr19_KI270916v1_alt', 'chr19_KI270917v1_alt', 'chr19_KI270918v1_alt', 'chr19_KI270919v1_alt', 'chr19_KI270920v1_alt', 'chr19_KI270921v1_alt', 'chr19_KI270922v1_alt', 'chr19_KI270923v1_alt', 'chr3_KI270924v1_alt', 'chr4_KI270925v1_alt', 'chr6_GL000252v2_alt', 'chr8_KI270926v1_alt', 'chr11_KI270927v1_alt', 'chr19_GL949748v2_alt', 'chr22_KI270928v1_alt', 'chr19_KI270929v1_alt', 'chr19_KI270930v1_alt', 'chr19_KI270931v1_alt', 'chr19_KI270932v1_alt', 'chr19_KI270933v1_alt', 'chr19_GL000209v2_alt', 'chr3_KI270934v1_alt', 'chr6_GL000253v2_alt', 'chr19_GL949749v2_alt', 'chr3_KI270935v1_alt', 'chr6_GL000254v2_alt', 'chr19_GL949750v2_alt', 'chr3_KI270936v1_alt', 'chr6_GL000255v2_alt', 'chr19_GL949751v2_alt', 'chr3_KI270937v1_alt', 'chr6_GL000256v2_alt', 'chr19_GL949752v1_alt', 'chr6_KI270758v1_alt', 'chr19_GL949753v2_alt', 'chr19_KI270938v1_alt', 'chrM', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrX', 'chrY', 'chrY_KI270740v1_random', 'lambda']

test_number = 11
for number in range(11):
    print(test_contigs[number])
    test_call_methylation = CallMethylation(input_file=f'{bsseeker_base_path}Gold_Zymo:Covaris:Illumina_TS:0:human.sorted.bam',
                                            genome_database=f'{bsseeker_base_path}hg38_reference',
                                            remove_ccgg=False,
                                            remove_sx_reads=False,
                                            min_read_depth=10,
                                            contig=test_contigs[number])
    test_call_methylation.call_methylation()



bs_align_args = ['python3', 'BSBolt-Call-Methylation.py', '-I',
                 f'{bsseeker_base_path}Gold_Zymo:Covaris:Illumina_TS:0:human.sorted.bam',
                 '-O', f'{bsseeker_base_path}Gold_Zymo:Covaris:Illumina_TS:0:human',
                  '-DB', f'{bsseeker_base_path}hg38_reference', '-t', '11']
print(' '.join(bs_align_args))
subprocess.run(bs_align_args)

bs_align_args = ['python3', 'BSBolt-Call-Methylation.py', '-I',
                 f'{bsseeker_base_path}Gold_Zymo:Covaris:Illumina_TS:0:human.sorted.bam',
                 '-O', f'{bsseeker_base_path}Gold_Zymo:Covaris:Illumina_TS:0:human',
                  '-DB', f'{bsseeker_base_path}hg38_reference', '-t', '4']
print(' '.join(bs_align_args))
subprocess.run(bs_align_args)
