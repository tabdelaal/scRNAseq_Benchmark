FROM python:3.7-slim-stretch

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

# Install python and pip deps
RUN apt-get update && apt-get install --no-install-recommends --yes build-essential python3-dev libxml2 libxml2-dev zlib1g-dev && \
    pip3 --no-cache-dir install --upgrade pip && \
    pip3 --no-cache-dir install --upgrade setuptools && \
    pip3 --no-cache-dir install pandas rpy2 scikit-learn statsmodels tensorflow Cell-BLAST && \
    apt-get autoremove --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY Scripts/run_Cell_BLAST.py /Scripts/
