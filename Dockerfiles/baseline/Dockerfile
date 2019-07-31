FROM debian:9.9-slim

# Install newest R version
RUN apt-get update && \
    apt-get install --no-install-recommends --yes wget gnupg apt-transport-https && \
    wget -qO - http://keys.gnupg.net/pks/lookup?op=get\&search=0xAD5F960A256A04AF | apt-key add - && \
    echo 'deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install --no-install-recommends --yes r-base && \
    apt-get purge --yes wget gnupg apt-transport-https && \
    apt-get autoremove --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install python
RUN apt-get update && \
    apt-get install --no-install-recommends --yes python3 python3-pip && \
    pip3 --no-cache-dir install setuptools && \
    pip3 --no-cache-dir install pandas rpy2 scikit-learn statsmodels && \
    rm -rf /var/lib/apt/lists/*

COPY Scripts/run_kNN50.py \
     Scripts/run_kNN9.py \
     Scripts/run_LDA.py \
     Scripts/run_LDA_rejection.py \
     Scripts/run_NMC.py \
     Scripts/run_RF.py \
     Scripts/run_SVM.py \
     Scripts/run_SVM_rejection.py \
     rank_gene_dropouts.py \
     /Scripts/
