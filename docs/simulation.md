## Bisulfite Sequencing Simulation

Simulating bisuflite reads is performed by assiging each cytosine and guanine in the user provided reference file a 
methylation value. This done based on an input methylation reference file (.CGmap or BSBolt Reference Directory) or 
sampling from a binomial distribution for CpG and CH sites. After a methylation value is assigned for cytosine and 
guanine individual reads simulated using ART are modified according to the methylation value of the nucleotides in the 
read. The methylation value represent the probability a base with be methylated, and as a result there random noise is 
introduced to more accurately represent bisulfite sequencing reads. 

**BSB Simulate Commands**
```
  -h, --help  show this help message and exit
  -G          Path for reference genome fasta file, fasta file should contain
              all contigs
  -A          Path to ART executable, default = bundled ART
  -O          Output prefix
  -PE         Simulate Paired End Reads, default Single End
  -RL         Simulated Read Lenghth
  -RD         Simulated Read Depth
  -U          Simulate Undirectional Reads, default=Directional
  -RC         Path to CGmap file to generate simulation reference profile
  -RO         Methylation reference output directory, default = output path
  -BR         Path to previously generate BSB simulation reference
  -IR1        Read 1 insertion rate
  -IR2        Read 2 insertion rate
  -DR1        Read 1 deletion rate
  -DR2        Read 2 deletion rate
  -NF         Cutoff theshold for a read with gaps, -, or ambiguous bases, N.
              Reads below the threshold will not be output
  -M          Mean paired end fragment size
  -SM         Paired end fragment length distribution standard deviation
  -SS         Sequencing system, default HS25
              ART Sequencing System Codes
              * HS10 - HiSeq 1000 (100bp) 
              * HS20 - HiSeq 2000 (100bp) 
              * HS25 - HiSeq 2500 (125bp, 150bp) 
              * HSXn - HiSeqX PCR free (150bp) 
              * HSXt - HiSeqX TruSeq (150bp) 
              * MinS - MiniSeq TruSeq (50bp) 
              * MSv1 - MiSeq v1 (250bp)
              * MSv3 - MiSeq v3 (250bp) 
              * NS50 - NextSeq500 v2 (75bp)
  -Q1         Optional read 1 quality profile, generated using ART Illumina
              read profiler
  -Q2         Path to read 2 quality profile

```
**Simulate Paired End, Undirectional Methylation Reads**
```bash
python3 BSBolt.py Simulate -G ~/Tests/TestData/BSB_test.fa -O ~/Tests/TestSimulations/BSB_pe -U -PE
```