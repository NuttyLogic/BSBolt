
Following bisulfite conversion, the sense (Watson) and anti-sense (Crick) reference strands are no longer complementary. 
Correctly aligning bisulfite converted reads requires the construction of an alignment index containing a bisulfite 
converted sequence for each reference strand.

BSB support generation of 3 types of alignment indices.

1. Whole genome bisulfite sequencing indices
    - Index is generated for complete Watson and Crick strands
2. Masked alignment indices
    - Sequence outside the reference region is masked in both the Watson and Crick strands before index generation
3. Reduced representation bisulfite sequencing (RRBS) indices
    - The reference sequence is *In silico* restriction enzyme digested to give a list of plausible target regions

#### **BSB Index Commands**

```shell
BSBolt Index -G {fasta reference} -DB {database output}

-h, --help     show this help message and exit

Index Options:
  -G File      path to reference genome fasta file, fasta file should contain all contigs
  -DB File     path to index directory, alignment index generated inside existing directory
               or new directory made if directory doesn't exist
  -B Int       block size for bwtsw algorithm,
               increasing block will speed up indexing and increase memory consumption [10000000]
  -MR File     path to bed file of mappable regions.
               Index will be built using using masked contig sequence [null]
  -IA          ignore alt contigs when constructing alignment index
  -rrbs        generate a Reduced Representative Bisulfite Sequencing (RRBS) index

  -rrbs-cut-format Str    Cut format to use for generation of RRBS database, [C-CGG] MSPI,
                          input multiple enzymes as a comma separate string, C-CGG,C-CGG,...
  -rrbs-lower Int      lower bound fragment size to consider RRBS index generation [40]
  -rrbs-upper Int      upper bound fragment size to consider RRBS index generation [500]
```

#### **WGBS Index Generation Example**

```shell
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB 
```

#### **Masked Alignment Index Generation Example**

```shell
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -MR /Tests/TestData/test_wgbs_madking.bed
```

#### **RRBS Index Generation Example**

```shell
# RRBS Index, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400
```
