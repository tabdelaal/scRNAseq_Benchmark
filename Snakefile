rule all:
  input:
    "{}/final_report".format(config["output_dir"])

"""
Rule for R based tools.
""" #TODO
rule run_R_wrapper:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    wrapper="Scripts/{tool}.R"
  output:
    "workspace/{tool}/output"
  conda: "Environments/{tool}.yml"
  shell:
    "Rscript Scripts/{wildcards.tool}.R"


"""
Rule for python based tools.
"""
rule run_python_wrapper:
  input:
    input_dir = config["input_dir"],
    datafile = "{}/{}".format(config["input_dir"], config["datafile"]),
    labfile = "{}/{}".format(config["input_dir"], config["labfile"]),
    folds = "{}/{}".format(config["input_dir"], config["folds"]),
    wrapper="Scripts/run_{tool}.py"
  output:
    pred = "{output_dir}/{tool}_[0]_pred.csv",
    true = "{output_dir}/{tool}_[0]_true.csv",
    test_time = "{output_dir}/{tool}_[0]_test_time.csv",
    training_time = "{output_dir}/{tool}_[0]_training_time.csv"
  conda: "Environments/{tool}.yml"
  shell:
    "python Scripts/run_{wildcards.tool}.py "
    "{input.input_dir} "
    "{wildcards.output_dir} "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds}"


"""
Rule for making the final report.
"""  #TODO
rule make_final_report:
  input:
    tool_outputs = expand("{output_dir}/{tool}/{tool}_[0]_true.csv",
        tool=config["tools_to_run"], output_dir=config["output_dir"])
  params:
    output_dir = config["output_dir"]
  output:
    "{}/final_report".format(config["output_dir"])
  shell:
    "touch {params.output_dir}/final_report"
