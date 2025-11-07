############################################################
# Dockerfile to build CRISPResso2
############################################################

#FROM continuumio/miniconda3
FROM mambaorg/micromamba:2.3.3

USER root

LABEL org.opencontainers.image.authors="support@edilytics.com"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN apt-get update && apt-get install gcc g++ bowtie2 samtools libsys-hostname-long-perl \
  -y --no-install-recommends \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && echo "deb http://deb.debian.org/debian trixie main contrib" > /etc/apt/sources.list.d/contrib.list \
  && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo


RUN micromamba install -c conda-forge -c bioconda -y -n base --debug fastp numpy cython jinja2 tbb=2020.2 pyparsing=2.3.1 setuptools scipy matplotlib-base seaborn pandas plotly upsetplot\
  && micromamba clean --all --yes

# install crispresso
COPY . /CRISPResso2
WORKDIR /CRISPResso2
RUN pip install . \
  && CRISPResso -h \
  && CRISPRessoBatch -h \
  && CRISPRessoPooled -h \
  && CRISPRessoWGS -h \
  && CRISPRessoCompare -h

ENTRYPOINT ["python","/CRISPResso2/CRISPResso2_router.py"]