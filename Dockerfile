############################################################
# Dockerfile to build CRISPResso2
############################################################

FROM continuumio/miniconda:4.5.12

# File Author / Maintainer
MAINTAINER Kendell Clement
RUN apt-get update && apt-get install gcc g++ python-matplotlib python-numpy bowtie2 samtools \
  -y --no-install-recommends \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
# conda installs
  && conda config --add channels defaults \
  && conda config --add channels conda-forge \
  && conda config --add channels bioconda \
  && conda install --debug biopython \
  && conda install --debug -c bioconda trimmomatic flash \
  && conda clean -ay

# install crispresso
COPY . /CRISPResso2
WORKDIR /CRISPResso2
RUN python setup.py install

ENTRYPOINT ["python","/CRISPResso2/CRISPResso2_router.py"]
