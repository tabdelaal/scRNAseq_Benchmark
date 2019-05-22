#TODO replace all 'latest' tags with actual verions

def feature_ranking(w):
    if "feature_ranking" in config.keys():
        return config["feature_ranking"]
    else:
        return "{output_dir}/rank_genes_dropouts.csv".format(
            output_dir=w.output_dir)


"""
Rule for making the final report.
"""  #TODO
rule make_final_report:
  input:
    tool_outputs = expand("{output_dir}/evaluation/Confusion/{tool}.csv",
        tool=config["tools_to_run"], output_dir=config["output_dir"])
  params:
    output_dir = config["output_dir"]
  output:
    "{}/final_report".format(config["output_dir"])
  shell:
    "touch {params.output_dir}/final_report"
    #TODO this should aggregate the evaluation results

rule evaluate:
  input:
    true="{output_dir}/{tool}/{tool}_true.csv",
    pred="{output_dir}/{tool}/{tool}_pred.csv"
  output:
    "{output_dir}/evaluation/Confusion/{tool}.csv",
    "{output_dir}/evaluation/F1/{tool}.csv",
    "{output_dir}/evaluation/PopSize/{tool}.csv",
    "{output_dir}/evaluation/Summary/{tool}.csv",
  log: "{output_dir}/evaluation/{tool}.log"
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "Rscript evaluate.R "
    "{input.true} "
    "{input.pred} "
    "{wildcards.output_dir}/evaluation "
    "{wildcards.tool} "
    "&> {log}"


"""
Rule for creating cross validation folds
"""
rule generate_CV_folds:
  input: config["labfile"],
  output: "{output_dir}/CV_folds.RData"
  log: "{output_dir}/CV_folds.log"
  params:
    column = config.get("column", 1) # default to 1
  singularity: "docker://scrnaseqbenchmark/cross_validation:latest"
  shell:
    "Rscript Cross_Validation.R "
    "{input} "
    "{params.column} "
    "{wildcards.output_dir} "
    "&> {log}"


"""
Rule for creating feature rank lists
"""
rule generate_dropouts_feature_rankings:
    input:
        datafile = config["datafile"],
        folds = "{output_dir}/CV_folds.RData"
    output: "{output_dir}/rank_genes_dropouts.csv"
    log: "{output_dir}/rank_genes_dropouts.log"
    singularity: "docker://scrnaseqbenchmark/baseline:latest"
    shell:
        "echo test > {wildcards.output_dir}/test\n"
        "python3 rank_gene_dropouts.py "
        "{input.datafile} "
        "{input.folds} "
        "{wildcards.output_dir} "
        "&> {log}"


"""
Rule for R based tools.
"""
rule singleCellNet:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/singleCellNet/singleCellNet_pred.csv",
    true = "{output_dir}/singleCellNet/singleCellNet_true.csv",
    test_time = "{output_dir}/singleCellNet/singleCellNet_test_time.csv",
    training_time = "{output_dir}/singleCellNet/singleCellNet_training_time.csv"
  log: "{output_dir}/singleCellNet/singleCellNet.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/singlecellnet:latest"
  shell:
    "Rscript Scripts/run_singleCellNet.R "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/singleCellNet "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

rule CHETAH:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/CHETAH/CHETAH_pred.csv",
    true = "{output_dir}/CHETAH/CHETAH_true.csv",
    total_time = "{output_dir}/CHETAH/CHETAH_total_time.csv"
  log: "{output_dir}/CHETAH/CHETAH.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/chetah:latest"
  shell:
    "Rscript Scripts/run_CHETAH.R "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/CHETAH "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

rule SingleR:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/SingleR/SingleR_pred.csv",
    true = "{output_dir}/SingleR/SingleR_true.csv",
    total_time = "{output_dir}/SingleR/SingleR_total_time.csv"
  log: "{output_dir}/SingleR/SingleR.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/singler:latest"
  shell:
    "Rscript Scripts/run_SingleR.R "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/SingleR "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"


"""
Rules for python based tools.
"""
rule kNN50:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/kNN50/kNN50_pred.csv",
    true = "{output_dir}/kNN50/kNN50_true.csv",
    test_time = "{output_dir}/kNN50/kNN50_test_time.csv",
    training_time = "{output_dir}/kNN50/kNN50_training_time.csv"
  log: "{output_dir}/kNN50/kNN50.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_kNN50.py "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/kNN50 "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

rule LDA:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/LDA/LDA_pred.csv",
    true = "{output_dir}/LDA/LDA_true.csv",
    test_time = "{output_dir}/LDA/LDA_test_time.csv",
    training_time = "{output_dir}/LDA/LDA_training_time.csv"
  log: "{output_dir}/LDA/LDA.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_LDA.py "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/LDA "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

rule NMC:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/NMC/NMC_pred.csv",
    true = "{output_dir}/NMC/NMC_true.csv",
    test_time = "{output_dir}/NMC/NMC_test_time.csv",
    training_time = "{output_dir}/NMC/NMC_training_time.csv"
  log: "{output_dir}/NMC/NMC.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_NMC.py "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/NMC "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

rule RF:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/RF/RF_pred.csv",
    true = "{output_dir}/RF/RF_true.csv",
    test_time = "{output_dir}/RF/RF_test_time.csv",
    training_time = "{output_dir}/RF/RF_training_time.csv"
  log: "{output_dir}/RF/RF.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_RF.py "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/RF "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"

rule SVM:
  input:
    datafile = config["datafile"],
    labfile = config["labfile"],
    folds = "{output_dir}/CV_folds.RData",
    ranking = feature_ranking
  output:
    pred = "{output_dir}/SVM/SVM_pred.csv",
    true = "{output_dir}/SVM/SVM_true.csv",
    test_time = "{output_dir}/SVM/SVM_test_time.csv",
    training_time = "{output_dir}/SVM/SVM_training_time.csv"
  log: "{output_dir}/SVM/SVM.log"
  params:
    n_features = config.get("number_of_features", 0)
  singularity: "docker://scrnaseqbenchmark/baseline:latest"
  shell:
    "python3 Scripts/run_SVM.py "
    "{input.datafile} "
    "{input.labfile} "
    "{input.folds} "
    "{wildcards.output_dir}/SVM "
    "{input.ranking} "
    "{params.n_features} "
    "&> {log}"
