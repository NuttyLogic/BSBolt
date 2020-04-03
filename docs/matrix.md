
BSBolt AggregrateMatrix takes a list of CGmap files, compressed or uncompressed, and assembles a consensus methylation matrix. Methylated sites that 
pass a read depth threshold and are present in a set proportion of samples are included in the matrix. 

**BSBolt AggregateMatrix Commands**
```shell
optional arguments:
  -h, --help            show this help message and exit
  -F                    Comma separated list of CGmap file paths, or path to
                        text file with list of line separated CGmap file paths
  -S                    Comma separated list of samples labels. If sample
                        labels are not provided sample labels are extracted
                        from CGmap file paths. Can also pass path to txt for
                        line separated sample labels.
  -min-coverage         Minimum site read depth coverage for a site to be
                        included in the aggregate matrix
  -min-sample           Proportion of samples that must have a valid site
                        (above minimum coverage threshold), for a site to
                        beincluded in the aggregate matrix.
  -O                    Aggregate matrix output path
  -CG                   Only output CG sites
  -verbose              Verbose aggregation
  -t                    Number of threads to use when assembling matrix
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


