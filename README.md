# scRNAseq_Benchmark
Benchmarking classification tools for scRNA-seq data

## To Add
1. List of used R and Python tools
2. Link to filtered datasets used

## How to use

```
docker run \
  -v <where the output should be place>:/benchmark/workspace \
  -v <where the inputs are>:/inputs \
  benchmark-sc-classifiers <config_file>
```

config_file:
```YML
dataset: <path to dataset>
tools_to_use: # tools to use
  - <tool 1>
  - <tool 2>
```

## Included tools
- TBD

## Adding new tools
#### Write a wrapper script
1. The script should take the following input:
   - TBD
1. The script should produce the following output:
   - TBD
1. The script must be either in python or R.

#### Adding the tools
1. Add a wrapper script to the `Scripts/` folder.
   1. Make sure the wrapper script is named exactly the same as the tool
      will be in the input tool_list. Ending in the appropriate extension for
      the used language: `.R` for R and `.py` for python.
1. Add a conda environment YAML to the `environments/` folder.
   1. Make sure the YAML file contains the `name` field, the value should be
      exactly the same as the name of the tool will be in the input tool_list.
   2. Make sure the environment YAML end in the `.yml` extension.
   3. Make sure that the environment contains all dependencies for both the
      tool itself and the wrapper script.

#### Rebuild the Docker image
From the root directory of this project (the directory with `run.sh`):
```
docker build . --tag benchmark-sc-classifiers:my-custom-tool --network host
```

> It may take a while to build the image. As all tools and there dependencies
need to get installed.
