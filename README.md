# DE_rpy2

Differential expression analysis: DESeq2, edgeR, limma. 

Realized in python based on rpy2

## Download R and required libraries

R: https://www.r-project.org/

Libraries: DESeq2, edgeR, limma 

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("limma")
```

## Download python packages

rpy2, matplotlib_venn, matplotlib

```
pip install -r requirements.txt
```

## Input data format

**count matrix:**

pandas.DataFrame object. 'id' is a column not index.

| id    | sample1 | sample2 | sample3 | sample4 | sample5 | sample6 |
| ----- | ------- | ------- | ------- | ------- | ------- | ------- |
| gene1 | 808.77  | 878.44  | 1017.81 | 626.82  | 535.66  | 569.03  |
| gene2 | 26.72   | 39.76   | 41.99   | 61.06   | 54.86   | 57.74   |
| gene3 | 46.1    | 64.35   | 74.96   | 67.84   | 63.18   | 66.24   |

**design matrix:**

sample name as index. 
It's OK not to contain this index as you give the conditions in the same order of the samples in count matrix's columns.

|         | condition |
| ------- | --------- |
| sample1 | treated   |
| sample2 | treated   |
| sample3 | treated   |
| sample4 | untreated |
| sample5 | untreated |
| sample6 | untreated |

**design formula:**

string. Default= column name of design matrix.

'~ condition'

## example
see [example]('./example.ipynb')

