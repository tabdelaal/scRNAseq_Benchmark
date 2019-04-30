FROM continuumio/miniconda3:4.5.12
WORKDIR /benchmark
COPY . ./

RUN conda install snakemake=5.4.5 -c bioconda -c conda-forge

## The following would likely decrease runtime, as all packages will be cached
## so they won't have to be downloaded. It will also, however, greatly increase
## the size of the image and build time.
#RUN for ENVYML in $(ls environments/*.yml); do \
#      conda env create -f "$ENVYML"; \
#    done && \
#    conda env list | cut -f1 -d' ' | tail -n +3 > installed-tools

ENTRYPOINT ["snakemake"]
