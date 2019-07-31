FROM r-base:3.6.0

COPY Scripts/run_scID.R \
     Dockerfiles/scid/install_packages.R \
     /Scripts/

# Install R packages
RUN apt-get update && \
    apt-get install --no-install-recommends --yes make gcc g++ r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev && \
    Rscript --vanilla /Scripts/install_packages.R && \
    apt-get purge --yes make gcc g++ r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev && \
    apt-get autoremove --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
