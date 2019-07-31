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
number_of_features: <number of features to be used as input for the classification methods, 0 means all, defaults to 0>
genes: <path to gene name list, only needed for garnett_CV and Garnett_Pretrained>
human: <whether or not the data is human, true means human, false means mouse, defaults to true>
tools_to_run: # List of tools to run
  - <tool 1>
  - <tool 2>
  - <...>
```

##### Tool specific inputs
Some tools require specific inputs. Add the following to your config file when
one of these tools:
- Garnett_CV
  ```YML
  Garnett_CV:
    markers: <path to Gernett marker gene file>
  ```
- Garnett_Pretrained
  ```YML
  Garnett_Pretrained:
    classifier: <path to Gernett classifier>
  ```

<!-- TODO explain these input files -->

## Included tools/methods
- kNN50
- kNN9
- LDA
- NMC
- RF
- SVM
- [singleCellNet](https://github.com/pcahan1/singleCellNet)
- [CHETAH](https://github.com/jdekanter/CHETAH)
- [scmap](https://github.com/hemberg-lab/scmap)
  - scmapcell
  - scmapcluster
- [SingleR](https://github.com/dviraran/SingleR)
- [scID](https://github.com/BatadaLab/scID)
- [scVI](https://github.com/YosefLab/scVI)
- [Cell_BLAST](https://github.com/gao-lab/Cell_BLAST)
- [Garnett](https://cole-trapnell-lab.github.io/garnett/)
  - Garnett_CV (without pretrained classifier)
  - Garnett_Pretrained (with pretrained classifier)

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

You will want to write a wrapper script for the tool you want to
add to facilitate this. The `"{output_dir}/CV_folds.RData"` input may be
used to provide your wrapper script with folds for cross_validation.
It is recommended to make a docker image containing all dependencies for both
the tool and any wrappers for the tool.  
This wrapper script should also make a selection of the features to be used.
This selection should be based on ranking which can be accessed by providing
`feature ranking` as input to the wrapper script. The number of features to be
used should be configurable and settable through the 'number_of_features' field
in the config.

The following can be used as a template for new rules. Replace everything
surrounded by (and including the) `<>` with appropriate values.
```
rule SVM:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/<tool name>/<tool name>_pred.csv",
    true = "{output_dir}/<tool name>/<tool name>_true.csv",
    test_time = "{output_dir}/<tool name>/<tool name>_test_time.csv",
    training_time = "{output_dir}/<tool name>/<tool name>_training_time.csv"
  log: "{output_dir}/<tool name>/<tool name>.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://<docker image>"
  shell:
    "<python or Rscript> <wrapper script> "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/<tool name> "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"
```
