Read alignments are reported in the standard SAM file format. The resulting alignment file contains [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
SAM tags and bisulfite specific SAM tags.

### **Bisulfite Sequencing SAM Fields**

Alignments to the Watson and Crick strands are reported as alignments to the sense and anti-sense strands respectively.
As a result checking for expected Illumina read orientation is not possible. Additionally, the the mapping quality
field conveys both the alignment mapping uniqueness and alignment bisulfite uniqueness. A read that is not bisulfite unique
will always be assigned a mapping quality of zero.  

### **Bisulfite Alignment Tags**

| Tag | Type    | Description | Example |
| :---: | :---:    | :---: | :---: |
|MD    |Z (string)| Mismatch position and bases, excluding bisulfite conversion mismatches| 10T90|
|XB    |Z (string)| Read bisulfite conversion position and context | 5zz1x3z11zZ15z8zZ6z2z11z3xz3z13Xz2|
|XC    |i (int)| Bisulfite conversion status, 1 if not fully converted, 0 if converted | XC:i:0|  
|YC    |i (integer)| Bisulfite ambiguous, 1 if ambiguous, not reported if non-ambiguous| YC:i:1 |
|YS    |Z (string)| Mapping strand (C=Crick, W=Watson) and alignment conversion pattern (C2T or G2A) | YS:Z:W_C2T (ie. Watson_Cytosine.to.Thymine)|
|XG | Z (string) | Mapping strand for MethylDackel compatibility (Crick=GA, Watson=CT)| XG:Z:CT |

### **MD and XB Tags**

Each alignment is compared to an unconverted reference sequence to generated the mismatch (MD) and bisulfite conversion (XB) tags.
The MD tag reports alignment mismatches if the mismatch can not be explained by bisulfite conversion. The read base
is reported for a mismatch or insertion, *^* is used for a deletion, and bases between mismatches are reported as counts.
The XB tag reports the bisulfite status of all methylatable bases within a read and the methylation context relative to the reference (see below),
with non-methylatable bases reported as counts.

| Character | Type (H=CTA) |
| :---: | :---:    |
|X| methylated CG|
|x| unmethylated CG|
|Y| methylated CHG|
|y| unmethylated CHG|
|Z| methylated CHH|
|z| unmethylated CHH|

### **XC Tag**

Read bisulfite conversion status is assessed by counting the number of methylated CH sites to total number of observed CH sites
within a read or read pair. A high proportion of methylated CH to total CH sites suggests the read has not been fully bisulfite converted.
The parameters used to assess bisulfite conversion status are set through the *-CT* and *-CP* alignment options. If not enough CH sites are
observed the XC tag is not reported.
