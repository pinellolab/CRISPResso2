############################################################
# Dockerfile to build CRISPResso2
############################################################

FROM mambaorg/micromamba:0.13.1

# File Author / Maintainer
MAINTAINER Kendell Clement
RUN apt-get update && apt-get install gcc g++ bowtie2 samtools \
  -y --no-install-recommends \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && micromamba install -c defaults -c conda-forge -c bioconda -y -n base --debug -c bioconda trimmomatic flash numpy cython jinja2 tbb=2020.2 \
  && micromamba clean --all --yes

#install ms fonts
RUN echo "deb http://httpredir.debian.org/debian buster main contrib" > /etc/apt/sources.list \
  && echo "deb http://security.debian.org/ buster/updates main contrib" >> /etc/apt/sources.list \
  && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo

# install crispresso
COPY . /CRISPResso2
WORKDIR /CRISPResso2
RUN python setup.py install \
  && CRISPResso -h \
  && CRISPRessoBatch -h \
  && CRISPRessoPooled -h \
  && CRISPRessoWGS -h \
  && CRISPRessoCompare -h


ENTRYPOINT ["python","/CRISPResso2/CRISPResso2_router.py"]
