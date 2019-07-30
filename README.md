# A comparison of automatic cell identification methods for single-cell RNA-sequencing data
We present a comprehensive evaluation of the performance of state-of-the-art classification methods, in addition to general-purpose classifiers, for automatic cell identification single cell RNA-sequencing datasets. Our goal is to provide the community with a fair evaluation of all available methods to facilitate the users’ choice as well as direct further developments to focus on the challenging aspects of automated cell type identification.

### Repository description
We provide all the scripts to run and evaluate all classifiers, and to reproduce the results introduced in the paper.


1. 'Scripts' folder contains a wrapper function to read the data and apply certain classification method.
2. ```Cross_Validation``` R script can be used to produce training and test indices for cross validation.
3. ```rank_gene_dropouts``` Python script can be used to apply feature selection using the dropout method, and rank genes accordingly.
4. ```evaluate``` R script can be used to evaluate the prediction of a certain classifier and obtain scores such as accuracy, median F1-score and % unlabeld cells.

For more details, please check function documentations.

### General Usage

To benchmark and fairly evaluate the performance of different classifiers using benchmark-datasets (Filtered datasets can be downloaded from https://zenodo.org/record/2877646#.XN8l__kzapo), apply the following steps:

#### Step 1

Apply the ```Cross_Validation``` R function on the corresponding dataset to obtain fixed training and test cell indices, straitified across different cell types. For example, using the Tabula Muris (TM) dataset

```R
Cross_Validation('~/TM/Labels.csv', 1, '~/TM/')
```

This command will create a ```CV_folds.RData``` file used as input in Step 2.

#### Step 2

Run each classifier wrapper. For example, running scPred on TM dataset

```R
run_scPred('~/TM/Filtered_TM_data.csv','~/TM/Labels.csv','~/TM/CV_folds.RData','~/Results/TM/')
```

This command will output the true and predicted cell labels as csv files, as well as the classifier computation time.

#### Step 3

Evaluate the classifier prediction by 

```R
result <- evaluate('~/Results/TM/scPred_True_Labels.csv', '~/Results/TM/scPred_Pred_Labels.csv')
```

This command will return the corresponding accuracy, median F1-score, F1-scores for all cell populations, % unlabeled cells, and confusion matrix.

### Usage with feature selection

#### Step 1

Apply the ```Cross_Validation``` R function on the corresponding dataset to obtain fixed training and test cell indices, straitified across different cell types. For example, using the Tabula Muris (TM) dataset

```R
Cross_Validation('~/TM/Labels.csv', 1, '~/TM/')
```

This command will create a ```CV_folds.RData``` file used as input in Step 2 and 3.

#### Step 2

Apply the ```rank_gene_dropouts``` Python script to get the genes ranking for each training fold using the dropout criteria

```
rank_gene_dropouts('~/TM/Filtered_TM_data.csv', '~/TM/CV_folds.RData', '~/TM/')
```

This command will create a ```rank_genes_dropouts.csv``` file used as input in Step 3.

#### Step 3

Run each classifier wrapper. For example, running scPred on TM dataset with 1000 genes

```R
run_scPred('~/TM/Filtered_TM_data.csv','~/TM/Labels.csv','~/TM/CV_folds.RData','~/Results/TM/',
GeneOrderPath = '~/TM/rank_genes_dropouts.csv',NumGenes = 1000)
```

This command will output the true and predicted cell labels as csv files, as well as the classifier computation time.

#### Step 4

Evaluate the classifier prediction by 

```R
result <- evaluate('~/Results/TM/scPred_True_Labels.csv', '~/Results/TM/scPred_Pred_Labels.csv')
```

This command will return the corresponding accuracy, median F1-score, F1-scores for all cell populations, % unlabeled cells, and confusion matrix.

### Evaluate Marker-based methods using DE genes

To evaluate the marker-based methods SCINA, DigitalCellSorter and Garnett using DE genes learned from the data, you may follow these steps:

#### Step 1

Apply the ```Cross_Validation``` R function on the corresponding dataset to obtain fixed training and test cell indices, straitified across different cell types. For example, using the Zheng_sorted dataset

```R
Cross_Validation('~/TM/Labels.csv', 1, '~/Zheng_sorted/')
```

This command will create a ```CV_folds.RData``` file used as input in Step 2 and 3.

#### Step 2

For each fold use the training data to get the DE genes using the ```DEgenesMAST``` R function, and pass these DE genes to the corresponding method, for example here we use SCINA, to obtain cell prediction for the test data.

```R
load('CV_folds.RData')
Data <- read.csv('~/Zheng_sorted/Filtered_DownSampled_SortedPBMC_data',row.names = 1)
Labels <- as.matrix(read.csv('~/Zheng_sorted/Labels.csv'))
Labels <- as.vector(Labels[,col_Index])
Data <- Data[Cells_to_Keep,]
Labels <- Labels[Cells_to_Keep]

for (i in c(1:n_folds))
{
    MarkerGenes <-  DEgenesMAST(t(Data[Train_Idx[[i]],]), Labels[Train_Idx[[i]]], Normalize = TRUE, LogTransform = TRUE)
    ## write the MarkerGenes into a marker genes file format, depending on the tested method, for example for SCINA
    write.csv(MarkerGenes, 'MarkerGenes.csv')
    ## run the SCINA wrapper using these DE marker genes
    run_SCINA(Data[Test_Idx[[i]],], Labels[Test_Idx[[i]]], 'MarkerGenes.csv', '~/Results/Zheng_sorted/')
}
```

### Snakemake

To support future extension of this benchmarking work with new classifiers and datasets, we provide a Snakemake workflow to automate the performed benchmarking analyses (https://github.com/tabdelaal/scRNAseq_Benchmark/tree/snakemake_and_docker).
