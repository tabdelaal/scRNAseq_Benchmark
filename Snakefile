rule all:
  input:
    "workspace/final_report"


"""
Rule for R based tools.
"""
rule run_R_wrapper:
  input:
    dataset=config["dataset"],
    wrapper="Scripts/{tool}.R"
  output:
    "workspace/{tool}/output"
  conda: "environments/{tool}.yml"
  shell:
    "Rscript Scripts/{wildcards.tool}.R"


"""
Rule for python based tools.
"""
rule run_python_wrapper:
  input:
    dataset=config["dataset"],
    wrapper="Scripts/{tool}.py"
  output:
    "workspace/{tool}/output"
  conda: "environments/{tool}.yml"
  shell:
    "python Scripts/{wildcards.tool}.py"


"""
Rule for making the final report.
"""
rule make_final_report:
  input:
    tool_outputs = expand("workspace/{tool}/output",
        tool=config["tools_to_run"])
  output:
    "workspace/final_report"
  shell:
    "touch workspace/final_report" #TODO
