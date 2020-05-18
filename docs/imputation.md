BSBolt Impute leverages the correlation structure between neighboring CpG sites to impute missing values through the use of a kNN sliding window.  
Within each window the nearest neighbors are calculated using Euclidean distance for non-null sites. The average value of k nearest neighbors is used to impute the null methylation value. 
To efficiently scale the algorithm, imputation can be performed in batches. 

![](img/kNN_graphic.png)

```shell
BSBolt Impute -M {BSBolt_matrix.txt} -O {imputed_matrix.txt}

-h, --help  show this help message and exit

Options:
  -M File     path to BSBolt matrix file
  -B Int      imputation sample batch size kNN imputation, by default the all of the samples will be processed as a single batch
  -W Int      sliding window size for imputation [3000000]
  -k Int      number of neighbors to use for imputation [5]
  -t Int      number of threads available for imputation [1]
  -verbose    verbose imputaiton
  -O File     output path for imputed matrix
  -R          randomize batches
```  

**Impute No Batches**
```shell
python3 -m BSBolt Impute -M ~/test_matrix.txt -W 100000 -k 3 -t 4 -O ~/test_matrix.impute.txt
```

**Batch Imputation**
```shell
python3 -m BSBolt Impute -M ~/test_matrix.txt -W 100000 -k 3 -t 4 -O ~/test_matrix.impute.txt -B 10 -R
```