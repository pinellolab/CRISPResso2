############################################################
# Dockerfile to build CRISPResso2
############################################################

FROM ghcr.io/prefix-dev/pixi:0.43.0 AS build

USER root

LABEL org.opencontainers.image.authors="support@edilytics.com"

# Install system dependencies and MS core fonts for plots
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y --no-install-recommends gcc g++ ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo

# Copy Pixi and project config first for layer caching
WORKDIR /CRISPResso2
COPY pixi.toml pyproject.toml ./
RUN pixi install

# Copy source and install the package
COPY . .
RUN pixi run install \
  && pixi run -- CRISPResso -h \
  && pixi run -- CRISPRessoBatch -h \
  && pixi run -- CRISPRessoPooled -h \
  && pixi run -- CRISPRessoWGS -h \
  && pixi run -- CRISPRessoCompare -h

ENTRYPOINT ["pixi", "run", "--", "python", "/CRISPResso2/CRISPResso2_router.py"]
