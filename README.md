# scRNAseq_Benchmark
Benchmarking classification tools for scRNA-seq data

## To Add
1. List of used R and Python tools
2. Link to filtered datasets used

## How to use
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and
[singularity](https://www.sylabs.io/docs/) need to be available on your system.

From the root of this repository:
```
snakemake \
  --configfile <configfile> \
  --use-singularity
```

If your data or output directory is not located under the root of this
repository, be sure to tell snakemake to mount the appropriate directories
in singularity:
```
snakemake \
  --configfile <configfile> \
  --use-singularity \
  --singularity-args '--bind <location of inputs>:<location of inputs> --bind <output directory>:<output directory>'
```

#### The config file
```YML
output_dir: <path to outputs directory>
datafile: <path to csv file with counts per cell>
labfile: <csv with true labels per cell>
column: <The index of the column in the labels file which ought to be used, defaults to 1>
tools_to_run: # List of tools to run
  - <tool 1>
  - <tool 2>
  - <...>

```
<!-- TODO explain these input files -->

## Included tools/methods
- kNN
- LDA
- NMC
- RF
- SVM
- [singleCellNet](https://github.com/pcahan1/singleCellNet)

## Adding new tools
TBD
