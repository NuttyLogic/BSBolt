import gzip
import os
import pickle
import subprocess
import unittest
import numpy as np
import time
from BSB.BSB_Utils.UtilityFunctions import get_external_paths
from BSB.BSB_Align.Align import BisulfiteAlignmentAndProcessing

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'
bt2, art = get_external_paths()
print('Building Methylation Index')
bsb_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_S']
#subprocess.run(bsb_index_commands)

bsb_database = f'{bsb_directory}Tests/TestData/BSB_Test_DB_S'
fastq1 = f'{bsb_directory}Tests/TestData/wgbs_pe_meth_1.fastq'
fastq2 = f'{bsb_directory}Tests/TestData/wgbs_pe_meth_2.fastq'
bt2_commands = {'fastq1': f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_1.fastq',
                'fastq2': f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_2.fastq',
                'undirectional_library': True,
                'bowtie2_commands':
                ['--quiet', '--no-mixed', '--no-discordant', '--sam-nohead', '--reorder',
                 '-k', '10', '-p', '12', '-L', '20', '-D', '40', '--end-to-end', '--score-min',
                 'L,-0.6,-0.6', '-I', '0', '-X', '500'],
                'bsb_database': f'{bsb_directory}Tests/TestData/BSB_Test_DB_S/',
                'bowtie2_path': f'{bsb_directory}External/BT2/bowtie2-2.3.4.3-linux-x86_64/bowtie2',
                'output_path': f'{bsb_directory}Tests/BSB_pe_test'}
start = time.time()
test_align = BisulfiteAlignmentAndProcessing(**bt2_commands)
test_align.align_reads()
print(time.time() - start)