FROM continuumio/miniconda3:24.1.2-0
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -y python=3.11 barrnap blast git
RUN conda install -y pip
RUN pip install git+https://github.com/quadram-institute-bioscience/socru.git
