# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CRISPResso2 is a bioinformatics pipeline for analyzing CRISPR/Cas9 genome editing outcomes from deep sequencing data. It aligns reads to reference amplicons, quantifies indels/mutations, and generates reports and visualizations.

**Documentation:** https://docs.crispresso.com

## Common Commands

### Installation (Development)
```bash
pip install -e .
```

### Running Tests

Unit tests (pytest):
```bash
pytest tests/unit_tests/

# With coverage
pytest tests --cov CRISPResso2

# Single test file
pytest tests/unit_tests/test_CRISPRessoCORE.py

# Single test
pytest tests/unit_tests/test_CRISPRessoCORE.py::test_function_name
```

Integration tests (from sibling `../CRISPResso2_tests` directory):
```bash
cd ../CRISPResso2_tests

make install       # Installs CRISPResso2 from ../CRISPResso2
make basic         # Run basic test only
make basic test    # Run basic test and diff against expected results
make all           # Run all integration tests
make all test      # Run all tests and diff against expected results
make clean         # Remove test output directories
```

Individual integration test targets:
- `make basic` - Basic core analysis (CRISPResso_on_FANC.Cas9)
- `make params` - Analysis with advanced parameters
- `make batch` - Batch processing
- `make pooled` - Pooled amplicons
- `make wgs` - WGS analysis
- `make compare` - Sample comparison
- `make aggregate` - Aggregate analysis
- `make prime-editor` - Prime editing analysis
- `make base_editor` - Base editing analysis

Add `test` target to diff results against expected: `make basic test`

### Running CRISPResso Tools

```bash
CRISPResso -r1 reads.fastq -a AMPLICON_SEQUENCE -g GUIDE_SEQUENCE
CRISPRessoBatch -bs batch_file.txt -a AMPLICON_SEQUENCE -g GUIDE_SEQUENCE
CRISPRessoPooled -r1 reads.fastq -f amplicons.txt
CRISPRessoWGS -b aligned.bam -r reference.fa -f regions.txt
CRISPRessoCompare sample1_dir/ sample2_dir/
CRISPRessoAggregate -p 'CRISPResso_on_*'
```

## Architecture

### Entry Points (Console Scripts)

Each tool has a corresponding `*CORE.py` module with a `main()` function:

| Command | Module |
|---------|--------|
| `CRISPResso` | `CRISPRessoCORE.py` |
| `CRISPRessoBatch` | `CRISPRessoBatchCORE.py` |
| `CRISPRessoPooled` | `CRISPRessoPooledCORE.py` |
| `CRISPRessoWGS` | `CRISPRessoWGSCORE.py` |
| `CRISPRessoCompare` | `CRISPRessoCompareCORE.py` |
| `CRISPRessoPooledWGSCompare` | `CRISPRessoPooledWGSCompareCORE.py` |
| `CRISPRessoAggregate` | `CRISPRessoAggregateCORE.py` |

### Core Modules

- **`CRISPRessoCORE.py`** (~8,600 lines) - Main analysis engine: read alignment, indel quantification, result aggregation
- **`CRISPRessoShared.py`** - Exception classes, logging utilities, version info, shared helper functions
- **`writers/vcf.py`** - VCF writing, alternate allele mapping, edit processing
- **`CRISPRessoPlot.py`** (~6,000 lines) - All matplotlib/seaborn visualizations
- **`CRISPRessoMultiProcessing.py`** - Parallel processing orchestration

### Cython Modules (Performance-Critical)

- **`CRISPResso2Align.pyx`** - Custom sequence alignment algorithms
- **`CRISPRessoCOREResources.pyx`** - Data structures including `ResultsSlotsDict` for efficient result storage

Pre-compiled `.so` files exist for macOS (x86_64, arm64) and Linux (x86_64). Rebuild with `pip install -e .` if modifying `.pyx` files.

### Report Generation

- **`CRISPRessoReports/CRISPRessoReport.py`** - Jinja2-based HTML report generation
- **Templates:** `CRISPRessoReports/templates/` - HTML templates for each tool type

### Parameter System

- **`args.json`** - Central parameter definitions for all tools. Contains argument names, types, defaults, help text, and which tools each parameter applies to.

## External Dependencies

Required system tools (for pooled/WGS analysis):
- `bowtie2` - Read alignment
- `samtools` - BAM file processing
- `fastp` - Quality filtering (optional)

## Design Documents

See `design_docs/` for detailed write-ups on specific subsystems and past debugging decisions:

- **`LEFT_NORMALIZATION.md`** - VCF indel left-normalization in `writers/vcf.py`: why it's needed, how the fix works, key data structures

## Key Constraints

- **Python 3 only**
- **numpy < 2** required (see test_env.yml)
- Cython build requires numpy headers at compile time
