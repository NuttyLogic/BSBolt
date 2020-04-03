BSBolt Impute leverages the correlation structure between neighboring CpG sites to impute missing values through the use of a kNN sliding window.  
Within each window the nearest neighbors are calculated using Euclidean distance for non-null sites. The average value of k nearest neighbors is used to impute the null methylation value. 
To efficiently scale the algorithm, imputation can be performed in batches. 

![](img/kNN_graphic.png)

```shell
  -h, --help  show this help message and exit
  -M          Path to BSB matrix file
  -B          Imputation sample batch size kNN imputation, by default the all
              of the samples will be processed as a single batch
  -W          Sliding window size for imputation
  -k          Number of neighbors to use for imputation, default = 5
  -t          Number of threads available for imputation
  -verbose    Verbose output
  -O          Output path for imputed matrix
  -R          Randomize batches
```  

**Impute No Batches**
```shell
python3 -m BSBolt Impute -M ~/test_matrix.txt -W 100000 -k 3 -t 4 -O ~/test_matrix.impute.txt
```

**Batch Imputation**
```shell
python3 -m BSBolt Impute -M ~/test_matrix.txt -W 100000 -k 3 -t 4 -O ~/test_matrix.impute.txt -B 10 -R
```