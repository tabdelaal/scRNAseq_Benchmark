# scRNAseq_Benchmark
Benchmarking classification tools for scRNA-seq data

## To Add
1. List of used R and Python tools
2. Link to filtered datasets used

## How to use

#### Using docker
This will require docker to be available on your system.
```
docker run \
  -v <where the output should end up>:<Output folder specified in the config> \
  -v <where the inputs are>:<inputs folder specified in the config> \
  benchmark-sc-classifiers \
  --configfile <config file> \
  <optionally additional snakemake flags>
```

#### Outside of docker
This will require
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and
[conda](https://conda.pydata.org/miniconda.html) to be available on your system.
```
snakemake \
  --snakefile <path to this repo's Snakefile> \
  --configfile <configfile> \
  --use-conda \
  <optionally additional snakemake flags>
```

#### The config file
```YML
input_dir: <path to inputs directory>
output_dir: <path to outputs directory>
datafile: <path to csv file with counts per cell>
labfile: <csv with true labels per cell>
folds: <Rdata file defining cross-validation folds>
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

## Adding new tools
#### Write a wrapper script
1. The script should take the following input:
   - TBD
1. The script should produce the following output:
   - TBD
1. The script must be either in python or R.

#### Adding the tools
1. Add a wrapper script to the `Scripts/` folder.
   1. Make sure the wrapper script is named exactly the same way as the tool
      will be in listed in the config file. Ending in the appropriate extension
      for the used language: `.R` for R and `.py` for python.
1. Add a conda environment YAML to the `Environments/` folder.
   1. Make sure environment YAML is named exactly the same way as the tool
      will be in listed in the config file. Ending with the `.yml` extension.
   1. Make sure that the environment contains all dependencies for both the
      tool itself and the wrapper script.

#### Rebuild the Docker image
From the root directory of this project (the directory with `run.sh`):
```
docker build . --tag benchmark-sc-classifiers:my-custom-tool --network host
```

> It may take a while to build the image. As all tools and there dependencies
need to get installed.
