FROM condaforge/miniforge3:latest AS base

LABEL maintainer="Andrew Page"
LABEL description="socru - genome structural typing via ribosomal operon arrangement"

# Install bioinformatics tools
RUN mamba install -y -c bioconda -c conda-forge \
    barrnap=0.9 \
    blast=2.15.0 \
    python=3.11 \
    && mamba clean -afy

# Install socru
COPY . /opt/socru
WORKDIR /opt/socru
RUN pip install --no-cache-dir .

ENTRYPOINT ["socru"]
CMD ["--help"]
