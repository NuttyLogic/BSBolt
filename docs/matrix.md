
BSBolt AggregrateMatrix takes a list of CGmap files, compressed or uncompressed, and assembles a consensus methylation matrix. Methylated sites that 
pass a read depth threshold and are present in a set proportion of samples are included in the matrix. 

**BSBolt AggregateMatrix Commands**
```shell
BSBolt AggregateMatrix -F {file1.CGmap,file2.CGmap,...} -O {output_matrix.txt}

-h, --help              show this help message and exit

Options:
  -F File,File,.        comma separated list of CGmap files, 
                        or path to text file with list of line separated CGmap files
  -S Str,Str,.          comma separated list of samples labels. If sample labels are not provided sample labels 
                        are extracted from CGmap files. Can also pass path to txt for line separated sample labels.
  -min-coverage Int     minimum site read depth coverage 
  -min-sample Flaot     proportion of samples that must have a valid site (above minimum coverage threshold)
  -O File               Aggregate matrix output path
  -CG                   Only output CG sites
  -verbose              Verbose aggregation
  -t Int                Number of threads to use when assembling matrix
  -count                Output a count matrix with count of methylated cytosines and total observed cytosines
```
**Aggregate Matrix Default Settings**

```shell
python3 -m BSBolt AggregateMatrix -F cgmap_1,cgmap_2,cgmap_3 -O ~/test_matrix.txt
```
**Aggregate Matrix Default Settings - File List**

```shell
python3 -m BSBolt AggregateMatrix -F cgmap_file_list.txt -O ~/test_matrix.txt
```

**Aggregate Matrix Default Settings - File List, Sample Labels, Verbose**

```shell
python3 -m BSBolt AggregateMatrix -F cgmap_file_list.txt -S sample1,sample2,sample3 -O ~/test_matrix.txt -verbose
```


