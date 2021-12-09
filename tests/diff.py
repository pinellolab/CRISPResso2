import argparse
import re
import sys
from difflib import unified_diff
from pathlib import Path
from os.path import basename, join


FLOAT_REGEXP = re.compile(r'\d+\.\d+')
IGNORE_FILES = frozenset([
    'CRISPResso_RUNNING_LOG.txt',
    'CRISPRessoBatch_RUNNING_LOG.txt',
    'CRISPRessoPooled_RUNNING_LOG.txt',
    'CRISPRessoWGS_RUNNING_LOG.txt',
    'CRISPRessoCompare_RUNNING_LOG.txt',
])


def round_float(f):
    return str(round(float(f.group(0)), 3))


def diff(file_a, file_b):
    with open(file_a) as fh_a, open(file_b) as fh_b:
        lines_a = [FLOAT_REGEXP.sub(round_float, line) for line in fh_a]
        lines_b = [FLOAT_REGEXP.sub(round_float, line) for line in fh_b]
        return list(unified_diff(lines_a, lines_b))


def diff_dir(dir_a, dir_b):
    files_a = {basename(f): f for f in Path(dir_a).rglob('*.txt')}
    files_b = {basename(f): f for f in Path(dir_b).rglob('*.txt')}
    diff_exists = False
    for file_basename_a, file_path_a in files_a.items():
        if file_basename_a in IGNORE_FILES:
            continue
        if file_basename_a in files_b:
            diff_results = diff(file_path_a, files_b[file_basename_a])
            if diff_results:
                print('Comparing {0} to {1}'.format(
                    file_path_a, files_b[file_basename_a],
                ))
                for result in diff_results:
                    print(result, end='')
                diff_exists |= True
        else:
            print('{0} is not in {1}'.format(file_basename_a, dir_b))
            diff_exists |= True

    for file_basename_b in files_b.keys():
        if file_basename_b not in files_a:
            print('{0} is not in {1}'.format(file_basename_b, dir_a))
            diff_exists |= True

    return diff_exists


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dir_a', help='Directory of text files to compare.')
    parser.add_argument('--dir_b', help='Other directory of text files to compare.')
    parser.add_argument('--expected_prefix', default='expectedResults')

    args = parser.parse_args()

    if args.dir_b is None:
        dir_b = join(args.expected_prefix, args.dir_a)
    else:
        dir_b = args.dir_b

    if diff_dir(args.dir_a, dir_b):
        sys.exit(1)
    sys.exit(0)
