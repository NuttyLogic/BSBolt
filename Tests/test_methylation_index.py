
import os
import subprocess
import unittest
from BSB.BSB_Index.ProcessCutSites import ProcessCutSites
from BSB.BSB_Index.WholeGenomeBuild import WholeGenomeIndexBuild
from BSB.BSB_Index.RRBSGenomeBuild import RRBSGenomeIndexBuild

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'


test_masked_build = WholeGenomeIndexBuild(genome_database=f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs',
                                          reference_file=f'{bsb_directory}Tests/TestData/BSB_test.fa',
                                          mappable_regions=f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed',
                                          bowtie2_path='bowtie2')
# generate methylation indices for rrbs, wgbs, and wgbs masked

bsb_wgbs_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                           '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs']

#subprocess.run(bsb_wgbs_index_commands)

bsb_wgbs_masked_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                                  '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs_masked', '-MR',
                                  f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed']

subprocess.run(bsb_wgbs_masked_index_commands)


bsb_rrbs_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                           '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_rrbs', '-rrbs', '-MR',
                           f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed']

subprocess.run(bsb_rrbs_index_commands)


# get unmasked sequence ranges