############################################################
# Dockerfile to build CRISPResso2
############################################################

# Set the base image to anaconda python 2.7
#FROM continuumio/anaconda:5.1.0
FROM continuumio/anaconda:5.0.0.1

# File Author / Maintainer
MAINTAINER Kendell Clement

ENV SHELL bash

RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

#RUN conda install -c bioconda trimmomatic
#RUN conda install -c bioconda flash
#RUN conda install -c anaconda seaborn

#Add build tools
RUN apt-get update && apt-get install build-essential default-jre samtools bowtie2 make gcc g++ zlib1g-dev zlib1g unzip -y 

RUN conda install biopython
RUN conda install -c bioconda flash 
RUN conda update matplotlib
RUN conda install -c bioconda trimmomatic

COPY . /CRISPResso

ENV PYTHONPATH /opt/conda/lib/python2.7/site-packages/numpy/core/include

WORKDIR /CRISPResso/CRISPResso
RUN /opt/conda/bin/python setup.py build_ext --inplace
RUN mv CRISPResso/*.so .
RUN rmdir CRISPResso/

WORKDIR /CRISPResso


# Reroute to enable the STREAM CLI and STREAM webapp
ENTRYPOINT ["/opt/conda/bin/python", "/CRISPResso/CRISPResso2_router.py"]
