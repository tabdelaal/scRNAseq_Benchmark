docker_images = {
    "LDA": "baseline"
}

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
Rule for python based tools.
"""
rule run_python_wrapper:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{output_dir}/CV_folds.RData",
    wrapper="Scripts/run_{tool}.py"
  output:
    pred = "{output_dir}/{tool}/{tool}_pred.csv",
    true = "{output_dir}/{tool}/{tool}_true.csv",
    test_time = "{output_dir}/{tool}/{tool}_test_time.csv",
    training_time = "{output_dir}/{tool}/{tool}_training_time.csv"
  log: "{output_dir}/{tool}/{tool}.log"
  singularity: "docker://"
  shell:
    "python Scripts/run_{wildcards.tool}.py "
    "{input.input_dir} "
    "{wildcards.output_dir}/{wildcards.tool} "
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
