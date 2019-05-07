FROM debian:9-slim

COPY Scripts/run_baseline.py \
     Scripts/run_kNN.py \
     Scripts/run_LDA.py \
     Scripts/run_NMC.py \
     Scripts/run_RF.py \
     Scripts/run_SVM.py \
     /Scripts/

RUN apt-get update && \
    apt-get install --no-install-recommends --yes r-base python3 python3-pip && \
    pip3 install setuptools && \
    pip3 install pandas rpy2 scikit-learn && \
    rm -rf /var/lib/apt/lists/*
