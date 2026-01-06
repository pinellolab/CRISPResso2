[![Docker Image Version (tag)](https://img.shields.io/docker/v/pinellolab/crispresso2/latest?logo=docker&label=Docker)](https://hub.docker.com/r/pinellolab/crispresso2/tags)
[![CircleCI branch](https://img.shields.io/circleci/project/github/pinellolab/CRISPResso2/master.svg)](https://circleci.com/gh/pinellolab/CRISPResso2)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crispresso2/README.html)
[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://docs.crispresso.com)

# CRISPResso2

CRISPResso2 is a software pipeline designed to enable rapid and intuitive interpretation of genome editing experiments. A limited web implementation is available at: https://crispresso2.pinellolab.org/.

Briefly, CRISPResso2:

- aligns sequencing reads to a reference sequence
- quantifies insertions, mutations and deletions to determine whether a read is modified or unmodified by genome editing
- summarizes editing results in intuitive plots and datasets

Access the full documentation at <https://docs.crispresso.com>.

In addition, CRISPResso can be run as part of a larger tool suite:

- [CRISPRessoBatch](https://docs.crispresso.com/suite/batch/tool.html) - for analyzing and comparing multiple experimental conditions at the same site
- [CRISPRessoPooled](https://docs.crispresso.com/suite/pooled/tool.html) - for analyzing multiple amplicons from a pooled amplicon sequencing experiment
- [CRISPRessoWGS](https://docs.crispresso.com/suite/wgs/tool.html) - for analyzing specific sites in whole-genome sequencing samples
- [CRISPRessoCompare](https://docs.crispresso.com/suite/compare/tool.html) - for comparing editing between two samples (e.g., treated vs control)
- [CRISPRessoAggregate](https://docs.crispresso.com/suite/aggregate/tool.html) - for aggregating results from previously-run CRISPResso analyses

![CRISPResso2 Schematic](https://raw.githubusercontent.com/pinellolab/CRISPResso2/master/crispresso_schematic.png "CRISPResso2 Schematic")

## Installation

CRISPResso2 can be [installed](https://docs.crispresso.com/installation.html) in the following ways:

- [Bioconda](https://docs.crispresso.com/installation.html#bioconda)
  - [Bioconda on Apple Silicon](https://docs.crispresso.com/installation.html#bioconda-for-apple-silicon)
- [Docker](https://docs.crispresso.com/installation.html#docker)

## Examples

- [CRISPResso example runs](https://docs.crispresso.com/suite/core/examples.html)
- [CRISPRessoBatch example runs](https://docs.crispresso.com/suite/batch/examples.html)
- [CRISPRessoPooled example runs](https://docs.crispresso.com/suite/pooled/examples.html)
- [CRISPRessoWGS example runs](https://docs.crispresso.com/suite/wgs/examples.html)
- [CRISPRessoCompare example runs](https://docs.crispresso.com/suite/compare/examples.html)
- [CRISPRessoPooledWGCompare example runs](https://docs.crispresso.com/suite/pooledwgscompare/examples.html)
- [CRISPRessoAggregate example runs](https://docs.crispresso.com/suite/aggregate/examples.html)

## Troubleshooting

If you run into any issues, check out the [Troubleshooting page](https://docs.crispresso.com/troubleshooting.html) or submit a [new discussion](https://github.com/pinellolab/CRISPResso2/discussions/new?category=troubleshooting).
