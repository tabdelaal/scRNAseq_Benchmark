rule all:
  input:
    "{}/final_report".format(config["output_dir"])

"""
Rule for creating cross validation folds
"""
rule folds:
  input: "{}/{}".format(config["input_dir"], config["labfile"]),
  output: "{output_dir}/CV_folds.RData"
  log: "{output_dir}/CV_folds.log"
  params:
    column = config["column"]
  conda: "Environments/Cross_Validation.yml"
  shell:
    "Rscript Scripts/Cross_Validation.R "
    "{input} "
    "{params.column} "
    "{wildcards.output_dir} &> {log}"


"""
Rule for R based tools.
"""
rule run_R_wrapper:
  input:
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
    wrapper="Scripts/Run_{tool}.R"
  output:
    pred = "{output_dir}/{tool}/{tool}_pred.csv",
    true = "{output_dir}/{tool}/{tool}_true.csv",
    test_time = "{output_dir}/{tool}/{tool}_test_time.csv",
    training_time = "{output_dir}/{tool}/{tool}_training_time.csv"
  conda: "Environments/{tool}.yml"
  log: "{output_dir}/{tool}{tool}.log"
  shell:
    "Rscript Scripts/Run_{wildcards.tool}.R "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/{wildcards.tool} &> {log}"


"""
Rules for python based tools.
"""
rule kNN:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
  output:
    pred = "{output_dir}/kNN/kNN_pred.csv",
    true = "{output_dir}/kNN/kNN_true.csv",
    test_time = "{output_dir}/kNN/kNN_test_time.csv",
    training_time = "{output_dir}/kNN/kNN_training_time.csv"
  log: "{output_dir}/kNN/kNN.log"
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_kNN.py "
    "{input.input_dir} "
    "{wildcards.output_dir}/kNN "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} &> {log}"

rule LDA:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
  output:
    pred = "{output_dir}/LDA/LDA_pred.csv",
    true = "{output_dir}/LDA/LDA_true.csv",
    test_time = "{output_dir}/LDA/LDA_test_time.csv",
    training_time = "{output_dir}/LDA/LDA_training_time.csv"
  log: "{output_dir}/LDA/LDA.log"
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_LDA.py "
    "{input.input_dir} "
    "{wildcards.output_dir}/LDA "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} &> {log}"

rule NMC:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
  output:
    pred = "{output_dir}/NMC/NMC_pred.csv",
    true = "{output_dir}/NMC/NMC_true.csv",
    test_time = "{output_dir}/NMC/NMC_test_time.csv",
    training_time = "{output_dir}/NMC/NMC_training_time.csv"
  log: "{output_dir}/NMC/NMC.log"
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_NMC.py "
    "{input.input_dir} "
    "{wildcards.output_dir}/NMC "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} &> {log}"

rule RF:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
  output:
    pred = "{output_dir}/RF/RF_pred.csv",
    true = "{output_dir}/RF/RF_true.csv",
    test_time = "{output_dir}/RF/RF_test_time.csv",
    training_time = "{output_dir}/RF/RF_training_time.csv"
  log: "{output_dir}/RF/RF.log"
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_RF.py "
    "{input.input_dir} "
    "{wildcards.output_dir}/RF "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} &> {log}"

rule SVM:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
  output:
    pred = "{output_dir}/SVM/SVM_pred.csv",
    true = "{output_dir}/SVM/SVM_true.csv",
    test_time = "{output_dir}/SVM/SVM_test_time.csv",
    training_time = "{output_dir}/SVM/SVM_training_time.csv"
  log: "{output_dir}/SVM/SVM.log"
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_SVM.py "
    "{input.input_dir} "
    "{wildcards.output_dir}/SVM "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} &> {log}"

"""
Rule for making the final report.
"""  #TODO
rule make_final_report:
  input:
    tool_outputs = expand("{output_dir}/{tool}/{tool}_true.csv",
        tool=config["tools_to_run"], output_dir=config["output_dir"])
  params:
    output_dir = config["output_dir"]
  output:
    "{}/final_report".format(config["output_dir"])
  shell:
    "touch {params.output_dir}/final_report"
