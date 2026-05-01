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
RUN pixi run -e default install

# --- Prune stage: trim non-runtime artifacts before final copy ---
FROM build AS runtime_env

RUN rm -rf /CRISPResso2/.pixi/envs/default/include \
  && rm -rf /CRISPResso2/.pixi/envs/default/conda-meta \
  && rm -rf /CRISPResso2/.pixi/envs/default/man \
  && rm -rf /CRISPResso2/.pixi/envs/default/share/man \
  && rm -rf /CRISPResso2/.pixi/envs/default/share/doc \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/cmake \
  && find /CRISPResso2/.pixi/envs/default -type d -name '__pycache__' -prune -exec rm -rf '{}' + \
  && find /CRISPResso2/.pixi/envs/default -type f -name '*.a' -delete \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/pip \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/pip-* \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/setuptools \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/setuptools-* \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/wheel \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/wheel-* \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/Cython \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/cython* \
  && rm -rf /CRISPResso2/.pixi/envs/default/lib/python*/site-packages/pyximport \
  && rm -f /CRISPResso2/.pixi/envs/default/bin/pip* \
  && rm -f /CRISPResso2/.pixi/envs/default/bin/cython* \
  && rm -f /CRISPResso2/.pixi/envs/default/bin/wheel

# --- Runtime stage: slim image with runtime environment only ---
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

# Keep env path identical to build stage to avoid shebang/entrypoint path issues.
COPY --from=runtime_env /CRISPResso2/.pixi/envs/default /CRISPResso2/.pixi/envs/default

# Copy only minimal app files required at runtime.
COPY --from=runtime_env /CRISPResso2/CRISPResso2 /CRISPResso2/CRISPResso2
COPY --from=runtime_env /CRISPResso2/CRISPResso2.egg-info /CRISPResso2/CRISPResso2.egg-info
COPY --from=runtime_env /CRISPResso2/CRISPResso2_router.py /CRISPResso2/CRISPResso2_router.py
COPY --from=runtime_env /CRISPResso2/LICENSE.txt /CRISPResso2/LICENSE.txt

WORKDIR /CRISPResso2

ENV PATH="/CRISPResso2/.pixi/envs/default/bin:${PATH}"

# Verify (disable .pyc generation to avoid adding cache files in this layer)
RUN PYTHONDONTWRITEBYTECODE=1 CRISPResso -h \
  && PYTHONDONTWRITEBYTECODE=1 CRISPRessoBatch -h \
  && PYTHONDONTWRITEBYTECODE=1 CRISPRessoPooled -h \
  && PYTHONDONTWRITEBYTECODE=1 CRISPRessoWGS -h \
  && PYTHONDONTWRITEBYTECODE=1 CRISPRessoCompare -h

ENTRYPOINT ["python", "/CRISPResso2/CRISPResso2_router.py"]
