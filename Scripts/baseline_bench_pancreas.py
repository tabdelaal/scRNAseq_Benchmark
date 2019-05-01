# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:23:12 2019

@author: Lieke
"""

from run_baseline import run_baseline

R = "CV_folds.RData"
input_dir = "/tudelft.net/staff-bulk/ewi/insy/DBL/Tamim/scRNAseq_Benchmark/Datasets/pancreatic_data/Human/Xin/Filtered Data/"
datafile = "Filtered_Xin_HumanPancreas_data.csv"
labfile = "Labels.csv"
output_dir = "/tudelft.net/staff-bulk/ewi/insy/DBL/Tamim/scRNAseq_Benchmark/Results/pancreatic_data/Human/Xin/"

run_baseline(input_dir,output_dir,datafile,labfile,R)



