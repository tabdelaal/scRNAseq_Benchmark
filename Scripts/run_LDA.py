#!/usr/bin/python3

from sys import argv

from run_baseline_all import run_baseline

run_baseline(argv[1], argv[2], argv[3], argv[4], ["LDA"], int(argv[5]), argv[6])
