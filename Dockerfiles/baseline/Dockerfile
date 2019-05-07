FROM debian:9-slim

COPY run_*.py /scripts/

RUN apt-get update && \
    apt-get install --no-install-recommends --yes r-base python3 python3-pip && \
    pip3 install setuptools && \
    pip3 install pandas rpy2 scikit-learn && \
    rm -rf /var/lib/apt/lists/*
