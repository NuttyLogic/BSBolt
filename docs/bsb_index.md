Following bisulfite conversion the Watson (sense) and  Crick (anti-sense) DNA strands are no longer complementary. 
Correctly aligning bisulfite converted reads requires the construction of an alignment index containing a bisulfite 
converted sequence for each reference strand; resulting in a larger index. Not all regions of the genome will be 
bisulfite unique, ie reads from these regions will not align uniquely to either the Watson or Crick strand. 
Reads mapping to different reference strands are handled during the alignment step.  

BSB support generation of 3 types of alignment indices. 

1. Whole genome bisulfite sequencing indices 
    - Index is generated for complete Watson and Crick strands
2. Masked alignment indices
    - Sequence outside the reference region is masked in both the Watson and Crick strands before index generation
3. Reduced representation bisulfite sequencing (RRBS) indices
    - The reference sequence is *In silico* restriction enzyme digested to give a list of plausible target regions

**BSB Index Commands**
```shell
  -h, --help            show this help message and exit
  -G                    Path to reference genome fasta file, fasta file should
                        contain all contigs
  -DB                   Path to index directory, will create directory if
                        folder does not exist
  -MR                   Path to bed file of mappable regions. Index will be
                        built using using masked contig sequence
  -rrbs                 Generate a Reduced Representative Bisulfite Sequencing
                        (RRBS) index
  -rrbs-cut-format      Cut format to use for generation of RRBS database,
                        default= C-CGG (MSPI), input multiple enzymes as a
                        comma seperate string, C-CGG,C-CGG,...
  -rrbs-lower           Lower bound fragment size to consider RRBS
                        indexgeneration, default = 40
  -rrbs-upper           Upper bound fragment size to consider RRBS
                        indexgeneration, default = 500
 
```
**WGBS Index Generation Example**
```shell
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB 
```

**Masked Alignment Index Generation Example**
```shell
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -MR /Tests/TestData/test_wgbs_madking.bed
```

**RRBS Index Generation Example**
```shell
# RRBS Index, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400
```