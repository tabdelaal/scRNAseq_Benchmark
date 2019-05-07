FROM debian:9-slim

COPY Scripts/Cross_Validation.R \
     Dockerfiles/cross_validation/install_packages.R \
     /Scripts/

RUN apt-get update && \
    apt-get install --no-install-recommends --yes wget gnupg apt-transport-https && \
    wget -qO - http://keys.gnupg.net/pks/lookup?op=get\&search=0xAD5F960A256A04AF | apt-key add - && \
    echo 'deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install --no-install-recommends --yes r-base make gcc libc6-dev g++ && \
    Rscript --vanilla /Scripts/install_packages.R && \
    apt-get purge --yes make gcc wget gnupg apt-transport-https g++ && \
    apt-get autoremove --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
