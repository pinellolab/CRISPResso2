"""Lightweight CLI entry points.

The ASCII header is printed before CORE (pandas/numpy) is imported, and
``--help`` / ``--version`` are served without importing CORE at all.
"""

import sys
from importlib import import_module

from CRISPResso2.startup import (
    __version__,
    getCRISPRessoArgParser,
    print_startup_header,
)

# Per-tool banner and argparse metadata. Header art matches the CORE modules.
_TOOLS = {
    'CRISPResso': {
        'module': 'CRISPResso2.CRISPRessoCORE',
        'parser_tool': 'Core',
        'parser_title': 'CRISPResso Parameters',
        'description': [
            '~~~CRISPResso 2~~~',
            '-Analysis of genome editing outcomes from deep sequencing data-',
        ],
        'header_str': None,
    },
    'CRISPRessoBatch': {
        'module': 'CRISPResso2.CRISPRessoBatchCORE',
        'parser_tool': 'Batch',
        'parser_title': 'CRISPRessoBatch Parameters',
        'description': [
            '~~~CRISPRessoBatch~~~',
            '-Analysis of CRISPR/Cas9 outcomes from batch deep sequencing data-',
        ],
        'header_str': r'''
 _________________
| __    ___ __    |
||__) /\ | /  |__||
||__)/--\| \__|  ||
|_________________|
        ''',
    },
    'CRISPRessoPooled': {
        'module': 'CRISPResso2.CRISPRessoPooledCORE',
        'parser_tool': 'Pooled',
        'parser_title': 'CRISPRessoPooled Parameters',
        'description': [
            '~~~CRISPRessoPooled~~~',
            '-Analysis of CRISPR/Cas9 outcomes from POOLED deep sequencing data-',
        ],
        'header_str': r'''
 _______________________
| __  __  __     __ __  |
||__)/  \/  \|  |_ |  \ |
||   \__/\__/|__|__|__/ |
|_______________________|
        ''',
    },
    'CRISPRessoWGS': {
        'module': 'CRISPResso2.CRISPRessoWGSCORE',
        'parser_tool': 'WGS',
        'parser_title': 'CRISPRessoWGS Parameters',
        'description': [
            '~~~CRISPRessoWGS~~~',
            '-Analysis of CRISPR/Cas9 outcomes from WGS data-',
        ],
        'header_str': r'''
 ____________
|     __  __ |
||  |/ _ (_  |
||/\|\__)__) |
|____________|
        ''',
    },
    'CRISPRessoCompare': {
        'module': 'CRISPResso2.CRISPRessoCompareCORE',
        'parser_tool': 'Compare',
        'parser_title': 'CRISPRessoCompare Parameters',
        'description': [
            '~~~CRISPRessoCompare~~~',
            '-Comparison of two CRISPResso analyses-',
        ],
        'header_str': r'''
 ___________________________
| __ __      __      __  __ |
|/  /  \|\/||__) /\ |__)|_  |
|\__\__/|  ||   /--\| \ |__ |
|___________________________|
        ''',
    },
    'CRISPRessoPooledWGSCompare': {
        'module': 'CRISPResso2.CRISPRessoPooledWGSCompareCORE',
        'parser_tool': None,
        'parser_title': 'CRISPRessoPooledWGSCompare Parameters',
        'description': [
            '~~~CRISPRessoPooledWGSCompare~~~',
            '-Comparison of two CRISPRessoPooled or CRISPRessoWGS analyses-',
        ],
        'header_str': r'''
 ____________________________________
| __  __  __     __ __        __  __ |
||__)/  \/  \|  |_ |  \ /|  |/ _ (_  |
||   \__/\__/|__|__|__// |/\|\__)__) |
|   __ __      __      __  __        |
|  /  /  \|\/||__) /\ |__)|_         |
|  \__\__/|  ||   /--\| \ |__        |
|____________________________________|
        ''',
    },
    'CRISPRessoAggregate': {
        'module': 'CRISPResso2.CRISPRessoAggregateCORE',
        'parser_tool': None,
        'parser_title': 'Aggregate CRISPResso2 Runs',
        'description': [
            '~~~CRISPRessoAggregate~~~',
            '-Aggregation of CRISPResso Run Data-',
        ],
        'header_str': r'''
___________________________________
|      __  __  _   _  __     ___ _ |
| /\  /__ /__ |_) |_ /__  /\  | |_ |
|/--\ \_| \_| | \ |_ \_| /--\ | |_ |
|__________________________________|
        ''',
    },
}


def _wants_help(args):
    return any(a in ('-h', '--help') for a in args)


def _wants_version(args):
    return args == ['--version'] or (args and args[0] == '--version')


def _run(prog):
    spec = _TOOLS[prog]
    args = sys.argv[1:]

    if _wants_version(args):
        print(f'{prog} {__version__}')
        return 0

    print_startup_header(spec['description'], spec['header_str'])

    if _wants_help(args) and spec['parser_tool'] is not None:
        parser = getCRISPRessoArgParser(spec['parser_tool'], parser_title=spec['parser_title'])
        parser.parse_args()
        return 0

    return import_module(spec['module']).main()


def main():
    """Run CRISPResso."""
    return _run('CRISPResso')


def main_batch():
    """Run CRISPRessoBatch."""
    return _run('CRISPRessoBatch')


def main_pooled():
    """Run CRISPRessoPooled."""
    return _run('CRISPRessoPooled')


def main_wgs():
    """Run CRISPRessoWGS."""
    return _run('CRISPRessoWGS')


def main_compare():
    """Run CRISPRessoCompare."""
    return _run('CRISPRessoCompare')


def main_pooled_wgs_compare():
    """Run CRISPRessoPooledWGSCompare."""
    return _run('CRISPRessoPooledWGSCompare')


def main_aggregate():
    """Run CRISPRessoAggregate."""
    return _run('CRISPRessoAggregate')
