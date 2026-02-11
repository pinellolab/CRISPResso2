############################################################
# Dockerfile to build CRISPResso2 (multi-stage)
############################################################

# --- Build stage: resolve deps and install package ---
FROM ghcr.io/prefix-dev/pixi:0.43.0 AS build

USER root

RUN apt-get update \
  && apt-get install -y --no-install-recommends gcc g++ \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /CRISPResso2
COPY pixi.toml pyproject.toml ./
RUN pixi install

COPY . .
RUN pixi run install

# --- Runtime stage: slim image with just the environment ---
FROM ubuntu:24.04

LABEL org.opencontainers.image.authors="support@edilytics.com"

# Install MS core fonts for plots
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y --no-install-recommends ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo

# Copy pixi binary from build stage
COPY --from=build /usr/local/bin/pixi /usr/local/bin/pixi

# Copy the project (includes .pixi/ environment) from build stage
COPY --from=build /CRISPResso2 /CRISPResso2

WORKDIR /CRISPResso2

# Verify
RUN pixi run -- CRISPResso -h \
  && pixi run -- CRISPRessoBatch -h \
  && pixi run -- CRISPRessoPooled -h \
  && pixi run -- CRISPRessoWGS -h \
  && pixi run -- CRISPRessoCompare -h

ENTRYPOINT ["pixi", "run", "--", "python", "/CRISPResso2/CRISPResso2_router.py"]
