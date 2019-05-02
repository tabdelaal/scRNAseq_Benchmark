FROM continuumio/miniconda3:4.5.12
WORKDIR /benchmark
COPY Environments/ Environments/
COPY Scripts/ Scripts/
COPY Snakefile Snakefile
COPY Build_scripts/ Build_scripts/

RUN mkdir workspace
RUN apt-get update && apt-get install libgfortran3 --yes
RUN conda install snakemake=5.4.5 -c bioconda -c conda-forge --yes

# Create the environments
RUN ENVIRONMENTS=$(ls Environments/ | sed 's/$/_environment_created/' | sed 's/^/Build_scripts\/\./') && \
    snakemake --use-conda --snakefile Build_scripts/Snakefile $ENVIRONMENTS

ENTRYPOINT ["snakemake", "--use-conda"]
