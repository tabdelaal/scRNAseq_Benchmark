# scRNAseq_Benchmark
Benchmarking classification tools for scRNA-seq data

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
- [CHETAH](https://github.com/jdekanter/CHETAH)

## Adding new tools
In order to add a tool to this benchmarking workflow, a rule for this tool
needs to be added to the `Snakefile`. This rule should produce as output:
- a table of predicted label (`<output directory/<tool>/<tool>_pred.csv`).
- a table of true labels (`<output directory/<tool>/<tool>_true.csv`).
- a tables of testing, prediction and/or total time:
  - `<output directory>/<tool>/<tool>_test_time.csv`
  - `<output directory>/<tool>/<tool>_training_time.csv`
  - `<output directory>/<tool>/<tool>_total_time.csv`

The input to this rule should be:
- a count table (specified as the `datafile` in the config).
- a true labels file (specified as the `labfile` in the config).

You will likely want to write a wrapper script for the tool you want to
add to facilitate this. The `"{output_dir}/CV_folds.RData"` input may be
used to provide your wrapper script with premade folds for cross_validation.
It is recommended to make a docker image containing all dependencies for both
the tool and any wrappers for the tool.

The following can be used as a template for new rules. Replace everything
surrounded by (and including the) `<>` with appropriate values.
```
rule <tool name>:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData"
  output:
    pred = "{output_dir}/<tool name>/<tool name>_pred.csv",
    true = "{output_dir}/<tool name>/<tool name>_true.csv",
    test_time = "{output_dir}/<tool name>/<tool name>_test_time.csv",
    training_time = "{output_dir}/<tool name>/<tool name>_training_time.csv",
    total_time = "{output_dir}/<tool name>/<tool name>_total_time.csv"
  log: "{output_dir}/<tool name>/<tool name>.log"
  singularity: "docker://<docker image>"
  shell:
    "<python or Rscript> <wrapper script> "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/<tool name> &> {log}"
```
