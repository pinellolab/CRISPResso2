'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''

import argparse
import datetime
import errno
import gzip
import json
import textwrap
import importlib.util
import importlib.metadata
from pathlib import Path

import numpy as np
import os
import pandas as pd
import re
import string
import shutil
import shlex
import signal
import subprocess as sb
import unicodedata
import logging
from inspect import getmodule, stack

from CRISPResso2 import CRISPResso2Align
from CRISPResso2 import CRISPRessoCOREResources

def read_version():
    return importlib.metadata.version('CRISPResso2')

__version__ = read_version()

###EXCEPTIONS############################
class FastpException(Exception):
    pass


class NoReadsAlignedException(Exception):
    pass


class AlignmentException(Exception):
    pass


class SgRNASequenceException(Exception):
    pass


class NTException(Exception):
    pass


class ExonSequenceException(Exception):
    pass


class DuplicateSequenceIdException(Exception):
    pass


class NoReadsAfterQualityFilteringException(Exception):
    pass


class BadParameterException(Exception):
    pass


class AutoException(Exception):
    pass


class OutputFolderIncompleteException(Exception):
    pass


class InstallationException(Exception):
    pass


class InputFileFormatException(Exception):
    pass


class PlotException(Exception):
    pass



#########################################

class StatusFormatter(logging.Formatter):
    def format(self, record):
        record.percent_complete = ''
        if record.args and 'percent_complete' in record.args:
            record.percent_complete = float(record.args['percent_complete'])
            self.last_percent_complete = record.percent_complete
        elif hasattr(self, 'last_percent_complete'): # if we don't have a percent complete, use the last one
            record.percent_complete = self.last_percent_complete
        else:
            record.percent_complete = 0.0
        record.json_message = record.getMessage().replace('\\', r'\\').replace('\n', r'\n').replace('"', r'\"')
        return super().format(record)


class StatusHandler(logging.FileHandler):
    def __init__(self, filename):
        super().__init__(filename, 'w')
        self.setFormatter(StatusFormatter('{\n  "message": "%(json_message)s",\n  "percent_complete": %(percent_complete)s\n}'))

    def emit(self, record):
        """Overwrite the existing file and write the new log."""
        if self.stream is None:  # log file is empty
            self.stream = self._open()
        else:  # log file is not empty, overwrite
            self.stream.seek(0)
        logging.StreamHandler.emit(self, record)
        self.stream.truncate()


class LogStreamHandler(logging.StreamHandler):
    def __init__(self, stream=None):
        super().__init__(stream)
        self.setFormatter(StatusFormatter(
            '%(levelname)-5s @ %(asctime)s (%(percent_complete).1f%% done):\n\t %(message)s \n',
            datefmt='%a, %d %b %Y %H:%M:%S',
        ))
        self.setLevel(logging.INFO)


def set_console_log_level(logger, level, debug=False):
    for handler in logger.handlers:
        if isinstance(handler, LogStreamHandler):
            if level == 4 or debug:
                handler.setLevel(logging.DEBUG)
            elif level == 3:
                handler.setLevel(logging.INFO)
            elif level == 2:
                handler.setLevel(logging.WARNING)
            elif level == 1:
                handler.setLevel(logging.ERROR)
            break


class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return list(map(
                lambda x: textwrap.fill(x, width, subsequent_indent=' ' * 24),
                text[2:].splitlines(),
            ))
        return argparse.HelpFormatter._split_lines(self, text, width)


def getCRISPRessoArgParser(tool, parser_title="CRISPResso Parameters"):
    parser = argparse.ArgumentParser(description=parser_title, formatter_class=CustomHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    # Getting the directory of the current script
    current_dir = Path(__file__).parent

    # Adjusting the path to point directly to the args.json file
    json_path = current_dir / 'args.json'

    with open(json_path, 'r') as json_file:
        args_dict = json.load(json_file)
        args_dict = args_dict["CRISPResso_args"]
    type_mapper = {
        "str": str,
        "int": int,
        "float": float,
    }

    for key, value in args_dict.items():
        # print(key, value)
        tools = value.get('tools', [])  # Default to empty list if 'tools' is not found
        if tool in tools:
            action = value.get('action')  # Use None as default if 'action' is not found
            required = value.get('required', False)  # Use False as default if 'required' is not found
            default = value.get('default')  # Use None as default if 'default' is not found
            type_value = value.get('type', 'str')  # Assume 'str' as default type if 'type' is not specified
            arg_help = value.get('help', '') if value.get('help') != "SUPPRESS" else argparse.SUPPRESS

            # Determine the correct function based on conditions
            if action:
                parser.add_argument(*value['keys'], help=arg_help, action=action)
            elif required and default is None:  # Checks if 'required' is true and 'default' is not provided
                parser.add_argument(*value['keys'], help=arg_help, type=type_mapper[type_value], required=True)
            elif required:  # Checks if 'required' is true (default is provided, as checked above)
                parser.add_argument(*value['keys'], help=arg_help, default=default, type=type_mapper[type_value], required=True)
            else:  # Case when neither 'action' nor 'required' conditions are met
                # Here, it handles the case where default might be None, which is a valid scenario
                kwargs = {'help': arg_help, 'type': type_mapper[type_value]}
                if default is not None: kwargs['default'] = default  # Add 'default' only if it's specified
                parser.add_argument(*value['keys'], **kwargs)
    return parser


def get_core_crispresso_options():
    parser = getCRISPRessoArgParser("Core")
    crispresso_options = set()
    d = parser.__dict__['_option_string_actions']
    for key in d.keys():
        d2 = d[key].__dict__['dest']
        crispresso_options.add(d2)

    return crispresso_options


def get_crispresso_options_lookup(tool):
    ##dict to lookup abbreviated params
    #    crispresso_options_lookup = {
    #    'r1':'fastq_r1',
    #    'r2':'fastq_r2',
    #    'a':'amplicon_seq',
    #    'an':'amplicon_name',
    #    .....
    # }
    crispresso_options_lookup = {}
    parser = getCRISPRessoArgParser(tool)
    d = parser.__dict__['_option_string_actions']
    for key in d.keys():
        d2 = d[key].__dict__['dest']
        key_sub = re.sub("^-*", "", key)
        if key_sub != d2:
            crispresso_options_lookup[key_sub] = d2
    return crispresso_options_lookup

def overwrite_crispresso_options(cmd, option_names_to_overwrite, option_values, paramInd=None, set_default_params=False, tool='Core'):
    """
    Updates a given command (cmd) by setting parameter options with new values in option_values.

    Parameters
    ----------
    cmd : str
        The command to run including original parameters
    option_names_to_overwrite : list
        List of options to overwrite e.g. crispresso options
    option_values : dict or Pandas DataFrame
        Values for the options to overwrite.
    paramInd : int, optional
        Index in dict - this is the run number in case of multiple runs.
        If paramInd is specified, option_values should be a DataFrame and the function will look up the value in the row with index paramInd.
        The default is None.
    set_default_params : bool, optional
        If True, for add values in option_values that are the same as the default values
        If False, default values will not be added to the command
    tool : str
        The CRISPResso tool to create the argparser - the params from this tool will be used
    Returns
    -------
    str
        The updated command with the new options set.
    -------
    """
    parser = getCRISPRessoArgParser(tool)
    this_program = cmd.split()[0]
    cmd = ' '.join(cmd.split()[1:])  # remove the program name from the command
    args = parser.parse_args(shlex.split(cmd)) # shlex split keeps quoted parameters together

    for option in option_names_to_overwrite:
        if option:
            if option in option_values:
                if paramInd is None:
                    if type(option_values) == dict:
                        val = option_values[option]
                    else:
                        val = getattr(option_values, option)
                else:
                    val = option_values.loc[paramInd, option]
                if val is None or str(val) == 'None':
                    pass
                elif str(val) == "True":
                    setattr(args, option, True)
                elif str(val) == "False":
                    setattr(args, option, False)
                elif isinstance(val, str):
                    if val != "":
                        setattr(args, option, str(val))
                elif isinstance(val, bool):
                    setattr(args, option, val)
                else:
                    setattr(args, option, str(val))

    # reconstruct the command
    new_cmd = this_program
    for action in parser._actions:
        if action.dest in args:
            val = getattr(args, action.dest)
            if not set_default_params and (val == action.default) or (str(val) == str(action.default)):
                continue
            if val is None or str(val) == "None":
                continue
            # argparse conveniently doesn't set type for bools - those action types are None
            if action.nargs == 0:
                if val: # if value is true
                    new_cmd += ' --%s' % action.dest
            elif action.type == bool: # but just in case...
                if val:
                    new_cmd += ' --%s' % action.dest
            elif action.type == str:
                if val != "":
                    if re.fullmatch(r"[a-zA-Z0-9\._]*", val): # if the value is alphanumeric, don't have to quote it
                        new_cmd += ' --%s %s' % (action.dest, val)
                    elif val.startswith('"') and val.endswith('"'):
                        new_cmd += ' --%s %s' % (action.dest, val)
                    else:
                        new_cmd += ' --%s "%s"' % (action.dest, val)
            elif action.type == int:
                new_cmd += ' --%s %s' % (action.dest, val)

    return new_cmd





def propagate_crispresso_options(cmd, options, params, paramInd=None):
    """
    Updates a given command (cmd) by setting parameter options with new values in params.
    This is used to propagate options from the command line to the command that is run.

    Parameters
    ----------
    cmd : str
        The command to run including original parameters
    options : list
        List of options to propagate e.g. crispresso options
    params : dict or Pandas DataFrame
        Values for the options to propagate.
    paramInd : int, optional
        Index in dict - this is the run number in case of multiple runs.
        If paramInd is specified, params should be a DataFrame and the function will look up the value in the row with index paramInd.
        The default is None.
    Returns
    -------
    str
        The updated command with the new options set.
    -------
    """

    for option in options:
        if option:
            if option in params:
                if paramInd is None:
                    if type(params) == dict:
                        val = params[option]
                    else:
                        val = getattr(params, option)
                else:
                    val = params.loc[paramInd, option]
                if val is None:
                    pass
                elif str(val) == "True":
                    cmd += ' --%s' % option
                elif str(val) == "False":
                    pass
                elif isinstance(val, str):
                    if val != "":
                        if re.match(r'-\d+$', val):
                            cmd += ' --%s %s' % (option, str(val))
                        elif " " in val or "-" in val:
                            cmd += ' --%s "%s"' % (option, str(val))  # quotes for options with spaces
                        else:
                            cmd += ' --%s %s' % (option, str(val))
                elif isinstance(val, bool):
                    if val:
                        cmd += ' --%s' % option
                else:
                    cmd += ' --%s %s' % (option, str(val))
    return cmd


#######
# Sequence functions
#######
nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-'})


def reverse_complement(seq):
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])


def reverse(seq):
    return "".join(c for c in seq.upper()[-1::-1])


def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference({'A', 'T', 'C', 'G', 'N'}))


def capitalize_sequence(x):
    return str(x).upper() if not pd.isnull(x) else x


def slugify(value):
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = re.sub(rb'[\s\'*"/\\\[\]:;|,<>?]', b'_', value).strip()
    value = re.sub(rb'_{2,}', b'_', value)

    return value.decode('utf-8')


CIGAR_LOOKUP = {
    ('A', 'A'): 'M', ('A', 'C'): 'M', ('A', 'T'): 'M', ('A', 'G'): 'M', ('A', 'N'): 'M',
    ('C', 'A'): 'M', ('C', 'C'): 'M', ('C', 'T'): 'M', ('C', 'G'): 'M', ('C', 'N'): 'M',
    ('T', 'A'): 'M', ('T', 'C'): 'M', ('T', 'T'): 'M', ('T', 'G'): 'M', ('T', 'N'): 'M',
    ('G', 'A'): 'M', ('G', 'C'): 'M', ('G', 'T'): 'M', ('G', 'G'): 'M', ('G', 'N'): 'M',
    ('N', 'A'): 'M', ('N', 'C'): 'M', ('N', 'T'): 'M', ('N', 'G'): 'M', ('N', 'N'): 'M',
    ('A', '-'): 'I', ('T', '-'): 'I', ('C', '-'): 'I', ('G', '-'): 'I', ('N', '-'): 'I',
    ('-', 'A'): 'D', ('-', 'T'): 'D', ('-', 'C'): 'D', ('-', 'G'): 'D', ('-', 'N'): 'D',
    }

cigarUnexplodePattern = re.compile(r'((\w)\2{0,})')


def unexplode_cigar(exploded_cigar_string):
    """Make a CIGAR string from an exploded cigar string.
    Exploded cigar: IIMMMMMMM
    cigar_els: ['2I','7M']
    CIGAR: 2I7M

    Parameters
    ----------
    exploded_cigar_string : str
        Exploded cigar string, e.g. IIMMMMMMM

    Returns
    -------
    cigar_els : list
        List of CIGAR elements


    """
    cigar_els = []
    for (cigar_str, cigar_char) in re.findall(cigarUnexplodePattern, exploded_cigar_string):
        cigar_els.append(str(len(cigar_str)) + cigar_char)
    return cigar_els


def get_ref_length_from_cigar(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    reference sequence.
    """
    read_consuming_ops = ("M", "D", "N", "=", "X")
    result = 0
    ops = re.findall(r'(\d+)(\w)', cigar_string)
    for c in ops:
        length, op = c
        if op in read_consuming_ops:
            result += int(length)
    return result


######
# File functions
######

def clean_filename(filename):
    # get a clean name that we can use for a filename
    validFilenameChars = "+-_.%s%s" % (string.ascii_letters, string.digits)
    filename = slugify(str(filename).replace(' ', '_'))
    cleanedFilename = unicodedata.normalize('NFKD', filename)
    return (''.join(c for c in cleanedFilename if c in validFilenameChars))

def check_file(filename):
    try:
        with open(filename):
            pass
    except IOError:
        files_in_curr_dir = os.listdir('.')
        if len(files_in_curr_dir) > 15:
            files_in_curr_dir = files_in_curr_dir[0:15]
            files_in_curr_dir.append("(Complete listing truncated)")
        dir_string = ""
        file_dir = os.path.dirname(filename)
        if file_dir == "":
            dir_string = ""
        elif os.path.isdir(file_dir):
            files_in_file_dir = os.listdir(file_dir)
            if len(files_in_file_dir) > 15:
                files_in_file_dir = files_in_file_dir[0:15]
                files_in_file_dir.append("(Complete listing truncated)")
            dir_string = "\nAvailable files in " + file_dir + ":\n\t" + "\n\t".join(files_in_file_dir)
        else:
            dir_string = "\nAdditionally, the folder '" + os.path.dirname(filename) + "' does not exist"

        raise BadParameterException(
            "The specified file '" + filename + "' cannot be opened.\nAvailable files in current directory:\n\t" + "\n\t".join(
                files_in_curr_dir) + dir_string)


def force_symlink(src, dst):
    if os.path.exists(dst) and os.path.samefile(src, dst):
        return

    try:
        os.symlink(src, dst)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)
        elif exc.errno == errno.EPROTO:
            # in docker on windows 7, symlinks don't work so well, so we'll just copy the file.
            shutil.copyfile(src, dst)


def parse_count_file(fileName):
    if os.path.exists(fileName):
        with open(fileName) as infile:
            lines = infile.readlines()
            ampSeq = lines[0].rstrip().split("\t")
            ampSeq.pop(0)  # get rid of 'Amplicon' at the beginning of line
            ampSeq = "".join(ampSeq)
            lab_freqs = {}
            for i in range(1, len(lines)):
                line = lines[i].rstrip()
                lab_freq_arr = line.split()
                lab = lab_freq_arr.pop(0)
                lab_freqs[lab] = lab_freq_arr
        return ampSeq, lab_freqs
    else:
        print("Cannot find output file '%s'" % fileName)
        return None, None


def parse_alignment_file(fileName):
    if os.path.exists(fileName):
        with open(fileName) as infile:
            lines = infile.readlines()
            ampSeq = lines[0].rstrip().split("\t")
            ampSeq.pop(0)  # get rid of 'Amplicon' at the beginning of line
            ampSeq = "".join(ampSeq)
            lab_freqs = {}
            for i in range(1, len(lines)):
                line = lines[i].rstrip()
                lab_freq_arr = line.split()
                lab = lab_freq_arr.pop(0)
                lab_freqs[lab] = lab_freq_arr
        return ampSeq, lab_freqs
    else:
        print("Cannot find output file '%s'" % fileName)
        return None, None

def assert_fastq_format(file_path, max_lines_to_check=100):
    """
    Checks to see that the fastq file is in the correct format
    Accepts files in gzipped format if they end in .gz.
    Raises a InputFileFormatException if the file is not in the correct format.
    params:
        file_path: path to fastq file
        max_lines_to_check: number of lines to check in the file
    returns:
        True if the file is in the correct format
    """

    try:
        if file_path.endswith('.gz'):
            # Read gzipped file
            with gzip.open(file_path, 'rt') as file:
                for line_num, line in enumerate(file):
                    if line_num >= max_lines_to_check:
                        break

                    if line_num % 4 == 0:
                        if not line.startswith('@'):
                            raise InputFileFormatException('File %s is not in fastq format! Line %d does not start with @ \n%s: %s' % (file_path, line_num, line_num, line))
                    elif line_num % 4 == 1 or line_num % 4 == 3:
                        if len(line.strip()) == 0:
                            raise InputFileFormatException('File %s is not in fastq format! Line %d is empty \n%s: %s' % (file_path, line_num, line_num, line))
                    elif line_num % 4 == 2:
                        if not line.startswith('+'):
                            raise InputFileFormatException('File %s is not in fastq format! Line %d does not start with + \n%s: %s' % (file_path, line_num, line_num, line))
        else:
            # Read uncompressed file
            with open(file_path, 'r') as file:
                for line_num, line in enumerate(file):
                    if line_num >= max_lines_to_check:
                        break

                    if line_num % 4 == 0:
                        if not line.startswith('@'):
                            raise InputFileFormatException('File %s is not in fastq format! Line %d does not start with @ \n%s: %s' % (file_path, line_num, line_num, line))
                    elif line_num % 4 == 1 or line_num % 4 == 3:
                        if len(line.strip()) == 0:
                            raise InputFileFormatException('File %s is not in fastq format! Line %d is empty \n%s: %s' % (file_path, line_num, line_num, line))
                    elif line_num % 4 == 2:
                        if not line.startswith('+'):
                            raise InputFileFormatException('File %s is not in fastq format! Line %d does not start with + \n%s: %s' % (file_path, line_num, line_num, line))

        return True

    except UnicodeDecodeError as e:
        raise InputFileFormatException('File %s is not in fastq format! Perhaps it is a gzipped file but does not end in .gz?' % (file_path)) from e
    except Exception as e:
        raise InputFileFormatException('File %s is not in fastq format!' % (file_path)) from e


def get_n_reads_fastq(fastq_filename):
    if not os.path.exists(fastq_filename) or os.path.getsize(fastq_filename) == 0:
        return 0

    p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | grep -c ." % fastq_filename, shell=True, stdout=sb.PIPE)
    n_reads = int(float(p.communicate()[0])/4.0)
    return n_reads


def check_output_folder(output_folder):
    """
    Checks to see that the CRISPResso run has completed, and gathers the amplicon info for that run
    returns:
    - quantification file = CRISPResso_quantification_of_editing_frequency.txt for this run
    - amplicons = a list of amplicons analyzed in this run
    - amplicon_info = a dict of attributes found in quantification_file for each amplicon
    """
    run_file = os.path.join(output_folder, 'CRISPResso2_info.json')
    if not os.path.exists(run_file):
        raise OutputFolderIncompleteException(
            'The folder %s is not a valid CRISPResso2 output folder. Cannot find summary file %s.' % (
            output_folder, run_file))
    with open(run_file) as fh:
        run_data = json.load(fh)

    amplicon_info = {}
    amplicons = run_data['results']['ref_names']

    quantification_file = os.path.join(output_folder, run_data['running_info']['quant_of_editing_freq_filename'])
    if os.path.exists(quantification_file):
        with open(quantification_file) as quant_file:
            head_line = quant_file.readline()
            head_line_els = head_line.split("\t")
            for line in quant_file:
                line_els = line.split("\t")
                amplicon_name = line_els[0]
                amplicon_info[amplicon_name] = {}
                amplicon_quant_file = os.path.join(output_folder, run_data['results']['refs'][amplicon_name][
                    'combined_pct_vector_filename'])
                if not os.path.exists(amplicon_quant_file):
                    raise OutputFolderIncompleteException(
                        'The folder %s is not a valid CRISPResso2 output folder. Cannot find quantification file %s for amplicon %s.' % (
                        output_folder, amplicon_quant_file, amplicon_name))
                amplicon_info[amplicon_name]['quantification_file'] = amplicon_quant_file

                amplicon_mod_count_file = os.path.join(output_folder, run_data['results']['refs'][amplicon_name][
                    'quant_window_mod_count_filename'])
                if not os.path.exists(amplicon_mod_count_file):
                    raise OutputFolderIncompleteException(
                        'The folder %s  is not a valid CRISPResso2 output folder. Cannot find modification count vector file %s for amplicon %s.' % (
                        output_folder, amplicon_mod_count_file, amplicon_name))
                amplicon_info[amplicon_name]['modification_count_file'] = amplicon_mod_count_file

                if 'allele_frequency_files' in run_data['results']['refs'][amplicon_name]:
                    amplicon_info[amplicon_name]['allele_files'] = [os.path.join(output_folder, x) for x in
                                                                    run_data['results']['refs'][amplicon_name][
                                                                        'allele_frequency_files']]
                else:
                    amplicon_info[amplicon_name]['allele_files'] = []

                for idx, el in enumerate(head_line_els):
                    amplicon_info[amplicon_name][el] = line_els[idx]

        return quantification_file, amplicons, amplicon_info
    else:
        raise OutputFolderIncompleteException(
            "The folder %s  is not a valid CRISPResso2 output folder. Cannot find quantification file '%s'." % (
            output_folder, quantification_file))


# Thanks https://gist.github.com/simonw/7000493 for this idea
class CRISPRessoJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, CRISPRessoCOREResources.ResultsSlotsDict):
            return {
                '_type': 'ResultsSlotsDict',
                'value': obj.__dict__,
            }
        if isinstance(obj, np.ndarray):
            return {
                '_type': 'np.ndarray',
                'value': obj.tolist(),
            }
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, pd.DataFrame):
            return {
                '_type': 'pd.DataFrame',
                'value': obj.to_json(orient='split'),
            }
        if isinstance(obj, datetime.datetime):
            return {
                '_type': 'datetime.datetime',
                'value': str(obj),
            }
        if isinstance(obj, datetime.timedelta):
            return {
                '_type': 'datetime.timedelta',
                'value': {
                    'days': obj.days,
                    'seconds': obj.seconds,
                    'microseconds': obj.microseconds,
                },
            }
        if isinstance(obj, set):
            return {
                '_type': 'set',
                'value': repr(obj),
            }
        if isinstance(obj, range):
            return {
                '_type': 'range',
                'value': repr(obj),
            }
        if isinstance(obj, argparse.Namespace):
            return {
                '_type': 'argparse.Namespace',
                'value': vars(obj),
            }
        return json.JSONEncoder.default(self, obj)


class CRISPRessoJSONDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(
            self,
            object_hook=self.object_hook,
            *args,
            **kwargs,
        )

    def object_hook(self, obj):
        if '_type' in obj:
            if obj['_type'] == 'ResultsSlotsDict':
                return CRISPRessoCOREResources.ResultsSlotsDict(**obj['value'])
            if obj['_type'] == 'np.ndarray':
                return np.array(obj['value'])
            if obj['_type'] == 'pd.DataFrame':
                return pd.read_json(obj['value'], orient='split')
            if obj['_type'] == 'datetime.datetime':
                return datetime.datetime.fromisoformat(obj['value'])
            if obj['_type'] == 'datetime.timedelta':
                return datetime.timedelta(
                    days=obj['value']['days'],
                    seconds=obj['value']['seconds'],
                    microseconds=obj['value']['microseconds'],
                )
            if obj['_type'] == 'set':
                return eval(obj['value'])
            if obj['_type'] == 'range':
                start, end, step = re.match(
                    r'range\((\d+), (\d+)(?:, (\d+))?\)', obj['value'],
                ).groups()
                if step is not None:
                    return range(int(start), int(end), int(step))
                return range(int(start), int(end))
            if obj['_type'] == 'argparse.Namespace':
                return argparse.Namespace(**obj['value'])
        return obj


def load_crispresso_info(
    crispresso_output_folder="",
    crispresso_info_file_name='CRISPResso2_info.json',
    crispresso_info_file_path=None
):
    """Load the CRISPResso2 info for a CRISPResso run.

        If crispresso_info_file_path is given, attempt to read info file at that location
        Otherwise, read file at crispresso_output_folder/crispresso_info_file_name

    Parameters
    ----------
    crispresso_output_folder : string
        Path to CRISPResso folder
    crispresso_info_file_name : string
        Name of info file in CRISPResso folder
    crispresso_info_path: string
        Path to info file

    Returns
    -------
    dict
        Dict of relevant information for a CRISPResso run

    """
    crispresso_info_file = os.path.join(
        crispresso_output_folder, crispresso_info_file_name,
    )
    if crispresso_info_file_path is not None:
        crispresso_info_file = crispresso_info_file_path

    if not os.path.isfile(crispresso_info_file):
        raise Exception('Cannot open CRISPResso info file at ' + crispresso_info_file)
    try:
        with open(crispresso_info_file) as fh:
            crispresso2_info = json.load(fh, cls=CRISPRessoJSONDecoder)
            return crispresso2_info
    except json.JSONDecodeError as e:
        raise Exception('Cannot open CRISPResso info file at ' + crispresso_info_file + "\n" + str(e))
    except (AttributeError, EOFError, ImportError, IndexError) as e:
        # secondary errors
        raise Exception('Cannot open CRISPResso info file at ' + crispresso_info_file + "\n" + str(e))
    except Exception as e:
        raise Exception('Cannot open CRISPResso info file at ' + crispresso_info_file + "\n" + str(e))


def write_crispresso_info(crispresso_output_file, crispresso2_info):
    """Write info hash to crispresso info output file.

    Parameters
    ----------
    crispresso_output_file : string
        File path to write info to
    crispresso2_info : dict
        Dict of relevant run properties

    Returns
    -------
    Nothing

    """
    with open(crispresso_output_file, 'w') as fh:
        json.dump(crispresso2_info, fh, cls=CRISPRessoJSONEncoder, indent=2)


def get_command_output(command):
    """
    Run a shell command and returns an iter to read the output.

    param:
        command: shell command to run

    returns:
        iter to read the output
    """
    p = sb.Popen(command,
                 stdout=sb.PIPE,
                 stderr=sb.STDOUT, shell=True,
                 #  encoding='utf-8',universal_newlines=True)
                 universal_newlines=True,
                 bufsize=-1)  # bufsize system default
    while True:
        retcode = p.poll()
        line = p.stdout.readline()
        yield line
        if retcode is not None:
            break


def get_most_frequent_reads(fastq_r1, fastq_r2, number_of_reads_to_consider, fastp_command, min_paired_end_reads_overlap, split_interleaved_input=False, debug=False):
    """
    Get the most frequent amplicon from a fastq file (or after merging a r1 and r2 fastq file).

    Note: only works on paired end or single end reads (not interleaved)

    input:
    fastq_r1: path to fastq r1 (can be gzipped)
    fastq_r2: path to fastq r2 (can be gzipped)
    number_of_reads_to_consider: number of reads from the top of the file to examine
    min_paired_end_reads_overlap: min overlap in bp for merging r1 and r2

    returns:
    list of amplicon strings sorted by order in format:
    12345 AATTCCG
    124 ATATATA
    5 TTATA
    """

    if split_interleaved_input:
        output_r1 = fastq_r1 + ".tmp.r1.gz"
        output_r2 = fastq_r2 + ".tmp.r2.gz"
        if fastq_r1.endswith('.gz'):
            fastq_handle = gzip.open(fastq_r1, 'rt')
        else:
            fastq_handle=open(fastq_r1)

        try:
            o1 = gzip.open(output_r1, 'wt')
            o2 = gzip.open(output_r2, 'wt')
            lines_read = 1
            r1_id = fastq_handle.readline()
            r1_seq = fastq_handle.readline()
            r1_plus = fastq_handle.readline()
            r1_qual = fastq_handle.readline()

            r2_id = fastq_handle.readline()
            r2_seq = fastq_handle.readline()
            r2_plus = fastq_handle.readline()
            r2_qual = fastq_handle.readline()
            while (r1_id and lines_read <= number_of_reads_to_consider):
                o1.write("%s%s%s%s"%(r1_id,r1_seq,r1_plus,r1_qual))
                o2.write("%s%s%s%s"%(r2_id,r2_seq,r2_plus,r2_qual))
                r1_id = fastq_handle.readline()
                r1_seq = fastq_handle.readline()
                r1_plus = fastq_handle.readline()
                r1_qual = fastq_handle.readline()

                r2_id = fastq_handle.readline()
                r2_seq = fastq_handle.readline()
                r2_plus = fastq_handle.readline()
                r2_qual = fastq_handle.readline()
            o1.close()
            o2.close()
            fastq_handle.close()
        except:
            raise BadParameterException('Error in splitting read pairs from a single file')
        fastq_r1 = output_r1
        fastq_r2 = output_r2

    view_cmd_1 = 'cat'
    if fastq_r1.endswith('.gz'):
        view_cmd_1 = 'gunzip -c'
    file_generation_command = "%s %s | head -n %d " % (view_cmd_1, fastq_r1, number_of_reads_to_consider * 4)

    if fastq_r2:
        view_cmd_2 = 'cat'
        if fastq_r2.endswith('.gz'):
            view_cmd_2 = 'gunzip -c'
        min_overlap_param = ""
        if min_paired_end_reads_overlap:
            min_overlap_param = "--overlap_len_require {0}".format(min_paired_end_reads_overlap)
        file_generation_command = "{paste} | {head} | paste - - - - | {awk} | {fastp}".format(
            paste="bash -c 'paste <({view_cmd_1} \"{fastq_r1}\") <({view_cmd_2} \"{fastq_r2}\")'".format(
                view_cmd_1=view_cmd_1, fastq_r1=fastq_r1, view_cmd_2=view_cmd_2, fastq_r2=fastq_r2,
            ),
            head='head -n {num_reads}'.format(
                num_reads=number_of_reads_to_consider * 4,
            ),
            awk="awk -v OFS=\"\\n\" -v FS=\"\\t\" '{{print($1,$3,$5,$7,$2,$4,$6,$8)}}'",
            fastp='{fastp_command} --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --stdin --interleaved_in --merge {min_overlap_param} --stdout 2>/dev/null'.format(
                fastp_command=fastp_command,
                min_overlap_param=min_overlap_param,
            ),
        )
    count_frequent_cmd = file_generation_command + " | awk '((NR-2)%4==0){print $1}' | sort | uniq -c | sort -nr "

    def default_sigpipe():
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    if (debug):
        print('command used: ' + count_frequent_cmd)

    piped_commands = count_frequent_cmd.split("|")
    pipes = [None] * len(piped_commands)
    pipes[0] = sb.Popen(piped_commands[0], stdout=sb.PIPE, preexec_fn=default_sigpipe, shell=True)
    for pipe_i in range(1, len(piped_commands)):
        pipes[pipe_i] = sb.Popen(piped_commands[pipe_i], stdin=pipes[pipe_i - 1].stdout, stdout=sb.PIPE,
                                 preexec_fn=default_sigpipe, shell=True)
    top_unaligned = pipes[-1].communicate()[0]

    if pipes[-1].poll() != 0:
        raise AutoException('Cannot retrieve most frequent amplicon sequences. Got nonzero return code.')
    seq_lines = top_unaligned.decode('utf-8').strip().split("\n")
    if len(seq_lines) == 0 or seq_lines == ['']:
        raise AutoException('Cannot parse any frequent amplicons sequences.')

    if split_interleaved_input:
        os.remove(output_r1)
        os.remove(output_r2)

    return seq_lines

def check_if_failed_run(folder_name, info):
    """
    Check the output folder for a info.json file and a status.txt file to see if the run completed successfully or not

    Parameters
    ----------
    folder_name: path to output folder
    info: logger

    Returns
    -------
    bool True if run failed, False otherwise
    string describing why it failed
    """

    run_data_file = os.path.join(folder_name, 'CRISPResso2_info.json')
    status_info = os.path.join(folder_name, 'CRISPResso_status.json')
    if not os.path.isfile(run_data_file) or not os.path.isfile(status_info):
        if not os.path.isfile(run_data_file):
            info("Skipping folder '%s'. Cannot find run data file at '%s'."%(folder_name, run_data_file))
        if not os.path.isfile(status_info):
            info("Skipping folder '%s'. Cannot find status file at '%s'."%(folder_name, status_info))
        if "CRISPRessoPooled" in folder_name:
            unit = "amplicon"
        elif "CRISPRessoWGS" in folder_name:
            unit = "region"
        else:
            unit = "sample"

        return True, f"CRISPResso failed for this {unit}! Please check your input files and parameters."
    else:
        with open(status_info) as fh:
            try:
                status_dict = json.load(fh)
                if status_dict['percent_complete'] != 100.0:
                    info("Skipping folder '%s'. Run is not complete (%s)." % (folder_name, status_dict['status']))
                    return True, str(status_dict['message'])
                else:
                    return False, ""
            except Exception as e:
                pass

        with open(status_info) as fh:
            try:
                file_contents = fh.read()
                search_result = re.search(r'(\d+\.\d+)% (.+)', file_contents)
                if search_result:
                    percent_complete, status = search_result.groups()
                    if percent_complete != '100.00':
                        info("Skipping folder '%s'. Run is not complete (%s)." % (folder_name, status))
                        return True, status
                else:
                    return True, file_contents
            except Exception as e:
                print(e)
                info("Skipping folder '%s'. Cannot parse status file '%s'." % (folder_name, status_info))
                return True, "Cannot parse status file '%s'." % (status_info)
        return False, ""


def guess_amplicons(fastq_r1, fastq_r2, number_of_reads_to_consider, fastp_command, min_paired_end_reads_overlap, aln_matrix, needleman_wunsch_gap_open, needleman_wunsch_gap_extend, split_interleaved_input=False, min_freq_to_consider=0.2, amplicon_similarity_cutoff=0.95):
    """
    guesses the amplicons used in an experiment by examining the most frequent read (giant caveat -- most frequent read should be unmodified)
    input:
    fastq_r1: path to fastq r1 (can be gzipped)
    fastq_r2: path to fastq r2 (can be gzipped)
    number_of_reads_to_consider: number of reads from the top of the file to examine
    fastp_command: command to call fastp
    min_paired_end_reads_overlap: min overlap in bp for merging r1 and r2
    aln_matrix: matrix specifying alignment substitution scores in the NCBI format
    needleman_wunsch_gap_open: alignment penalty assignment used to determine similarity of two sequences
    needleman_wunsch_gap_extend: alignment penalty assignment used to determine similarity of two sequences
    split_interleaved_input: if true, split interleaved input into two files
    min_freq_to_consider: selected ampilcon must be frequent at least at this percentage in the population
    amplicon_similarity_cutoff: if the current amplicon has similarity of greater than this cutoff to any other existing amplicons, it won't be added

    returns:
    list of putative amplicons
    """
    seq_lines = get_most_frequent_reads(fastq_r1, fastq_r2, number_of_reads_to_consider, fastp_command, min_paired_end_reads_overlap, split_interleaved_input=split_interleaved_input)

    curr_amplicon_id = 1

    amplicon_seq_arr = []

    # add most frequent amplicon to the list
    count, seq = seq_lines[0].strip().split()
    amplicon_seq_arr.append(seq)
    curr_amplicon_id += 1

    # for the remainder of the amplicons, test them before adding
    for i in range(1, len(seq_lines)):
        count, seq = seq_lines[i].strip().split()
        last_count, last_seq = seq_lines[i - 1].strip().split()
        # if this allele is present in at least XX% of the samples
        #        print('debug 509 testing ' + str(seq_lines[i]) + ' with ' + str(count) + ' out of consididered ' + str(number_of_reads_to_consider) + ' min freq: ' + str(min_freq_to_consider))
        if float(last_count) / float(number_of_reads_to_consider) > min_freq_to_consider:
            this_amplicon_seq_arr = amplicon_seq_arr[:]
            this_amplicon_max_pct = 0  # keep track of similarity to most-similar already-found amplicons
            for amp_seq in this_amplicon_seq_arr:
                ref_incentive = np.zeros(len(amp_seq) + 1, dtype=int)
                fws1, fws2, fwscore = CRISPResso2Align.global_align(seq, amp_seq, matrix=aln_matrix,
                                                                    gap_incentive=ref_incentive,
                                                                    gap_open=needleman_wunsch_gap_open,
                                                                    gap_extend=needleman_wunsch_gap_extend, )
                rvs1, rvs2, rvscore = CRISPResso2Align.global_align(reverse_complement(seq), amp_seq, matrix=aln_matrix,
                                                                    gap_incentive=ref_incentive,
                                                                    gap_open=needleman_wunsch_gap_open,
                                                                    gap_extend=needleman_wunsch_gap_extend, )
                # if the sequence is similar to a previously-seen read, don't add it
                min_len = min(len(last_seq), len(seq))
                max_score = max(fwscore, rvscore)
                if max_score / float(min_len) > this_amplicon_max_pct:
                    this_amplicon_max_pct = max_score / float(min_len)
            # if this amplicon was maximally-similar to all other chosen amplicons by less than amplicon_similarity_cutoff, add to the list
            if this_amplicon_max_pct < amplicon_similarity_cutoff:
                amplicon_seq_arr.append(seq)
                curr_amplicon_id += 1
        else:
            break

    return amplicon_seq_arr


def guess_guides(amplicon_sequence, fastq_r1, fastq_r2, number_of_reads_to_consider, fastp_command,
            min_paired_end_reads_overlap, exclude_bp_from_left, exclude_bp_from_right,
            aln_matrix, needleman_wunsch_gap_open, needleman_wunsch_gap_extend,
            min_edit_freq_to_consider=0.1, min_edit_fold_change_to_consider=3,
            pam_seq="NGG", min_pct_subs_in_base_editor_win=0.8, split_interleaved_input=False):
    """
    guesses the guides used in an experiment by identifying the most-frequently edited positions, editing types, and PAM sites
    input:
    ampilcon_sequence - amplicon to analyze
    fastq_r1: path to fastq r1 (can be gzipped)
    fastq_r2: path to fastq r2 (can be gzipped)
    number_of_reads_to_consider: number of reads from the top of the file to examine
    fastp_command: command to call fastp
    min_paired_end_reads_overlap: min overlap in bp for flashing (merging) r1 and r2
    exclude_bp_from_left: number of bp to exclude from the left side of the amplicon sequence for the quantification of the indels
    exclude_bp_from_right: number of bp to exclude from the right side of the amplicon sequence for the quantification of the indels
    aln_matrix: matrix specifying alignment substitution scores in the NCBI format
    needleman_wunsch_gap_open: alignment penalty assignment used to determine similarity of two sequences
    needleman_wunsch_gap_extend: alignment penalty assignment used to determine similarity of two sequences
    min_edit_freq_to_consider: edits must be at least this frequency for consideration
    min_edit_fold_change_to_consider: edits must be at least this fold change over background for consideration
    pam_seq: pam sequence to look for (can be regex or contain degenerate bases)
    min_pct_subs_in_base_editor_win: if at least this percent of substitutions happen in the predicted base editor window, return base editor flag
    split_interleaved_input: if true, interleaved fastq will be split into r1 and r2

    returns:
    tuple of (putative guide, boolean is_base_editor)
    or (None, None)
    """
    seq_lines = get_most_frequent_reads(fastq_r1, fastq_r2, number_of_reads_to_consider, fastp_command, min_paired_end_reads_overlap,split_interleaved_input=split_interleaved_input)

    amp_len = len(amplicon_sequence)
    gap_incentive = np.zeros(amp_len + 1, dtype=int)
    include_idxs = range(amp_len)
    exclude_idxs = []

    if exclude_bp_from_left:
        exclude_idxs += range(exclude_bp_from_left)
    if exclude_bp_from_right:
        exclude_idxs += range(amp_len)[-exclude_bp_from_right:]
    include_idxs = np.ravel(include_idxs)
    exclude_idxs = np.ravel(exclude_idxs)

    include_idxs = set(np.setdiff1d(include_idxs, exclude_idxs))

    all_indel_count_vector = np.zeros(amp_len)
    all_sub_count_vector = np.zeros(amp_len)
    tot_count = 0;
    for i in range(len(seq_lines)):
        count, seq = seq_lines[i].strip().split()
        count = int(count)
        tot_count += count
        fws1, fws2, fwscore = CRISPResso2Align.global_align(seq, amplicon_sequence, matrix=aln_matrix,
                                                            gap_incentive=gap_incentive,
                                                            gap_open=needleman_wunsch_gap_open,
                                                            gap_extend=needleman_wunsch_gap_extend, )
        payload = CRISPRessoCOREResources.find_indels_substitutions(fws1, fws2, include_idxs)
        all_indel_count_vector[payload['all_insertion_positions']] += count
        all_indel_count_vector[payload['all_deletion_positions']] += count
        all_sub_count_vector[payload['all_substitution_positions']] += count

    background_val = np.mean(all_indel_count_vector)
    if len(exclude_idxs) > 0:
        all_indel_count_vector[exclude_idxs] = 0
    max_loc = np.argmax(all_indel_count_vector)
    max_val = all_indel_count_vector[max_loc]

    # return nothing if the max edit doesn't break threshold
    if max_val / float(tot_count) < min_edit_freq_to_consider:
        return (None, None)

    # return nothing if the max edit doesn't break threshold over background
    if max_val / background_val < min_edit_fold_change_to_consider:
        return (None, None)

    pam_regex_string = pam_seq.upper()
    pam_regex_string = pam_regex_string.replace('I', '[ATCG]')
    pam_regex_string = pam_regex_string.replace('N', '[ATCG]')
    pam_regex_string = pam_regex_string.replace('R', '[AG]')
    pam_regex_string = pam_regex_string.replace('Y', '[CT]')
    pam_regex_string = pam_regex_string.replace('S', '[GC]')
    pam_regex_string = pam_regex_string.replace('W', '[AT]')
    pam_regex_string = pam_regex_string.replace('K', '[GT]')
    pam_regex_string = pam_regex_string.replace('M', '[AC]')
    pam_regex_string = pam_regex_string.replace('B', '[CGT]')
    pam_regex_string = pam_regex_string.replace('D', '[AGT]')
    pam_regex_string = pam_regex_string.replace('H', '[ACT]')
    pam_regex_string = pam_regex_string.replace('V', '[ACG]')

    is_base_editor = False
    # offset from expected position
    for offset in (0, +1, -1, +2, +3, +4, -2):
        # forward direction
        # find pam near max edit loc
        pam_start = max_loc + 4 + offset
        pam_end = max_loc + 7 + offset
        guide_start = max_loc - 16 + offset
        guide_end = max_loc + 4 + offset
        base_edit_start = max_loc - 16 + offset
        base_edit_end = max_loc - 6 + offset
        if pam_start > 0 and guide_end < amp_len:
            if re.match(pam_regex_string, amplicon_sequence[pam_start:pam_end]):
                guide_seq = amplicon_sequence[guide_start:guide_end]
                sum_base_edits = sum(all_sub_count_vector[base_edit_start:base_edit_end])
                # if a lot of edits are in the predicted base editor window, set base editor true
                # specifically, if at least min_pct_subs_in_base_editor_win % of substitutions happen in the predicted base editor window
                if sum_base_edits > min_pct_subs_in_base_editor_win * sum(all_sub_count_vector):
                    is_base_editor = True
                return guide_seq, is_base_editor

        # reverse direction
        pam_start = max_loc - 5 - offset
        pam_end = max_loc - 2 - offset
        guide_start = max_loc - 2 - offset
        guide_end = max_loc + 18 - offset
        base_edit_start = max_loc + 8 - offset
        base_edit_end = max_loc + 18 - offset
        if pam_start > 0 and guide_end < amp_len:
            if re.match(pam_regex_string, amplicon_sequence[pam_start:pam_end]):
                guide_seq = amplicon_sequence[guide_start:guide_end]
                sum_base_edits = sum(all_sub_count_vector[base_edit_start:base_edit_end])
                # if a lot of edits are in the predicted base editor window, set base editor true
                # specifically, if at least min_pct_subs_in_base_editor_win % of substitutions happen in the predicted base editor window
                if sum_base_edits > min_pct_subs_in_base_editor_win * sum(all_sub_count_vector):
                    is_base_editor = True
                return guide_seq, is_base_editor

    return (None, None)


######
# Fastq file manipulation
######

def force_merge_pairs(r1_filename, r2_filename, output_filename):
    """
    This can be useful in case paired end reads are too short to cover the amplicon.
    Note that this should be used with extreme caution because non-biological indels will appear at the site of read merging.
    R1------>     <-------R2
    becomes
    R1------><------R2

    input:
    r1_filename: path to fastq r1 (can be gzipped)
    r2_filename: path to fastq r2 (can be gzipped)
    output_filename: path to merged output filename

    returns:
    linecount: the number of lines of the resulting file
    """

    if r1_filename.endswith('.gz'):
        f1 = gzip.open(r1_filename, 'rt')
    else:
        f1 = open(r1_filename, 'r')
    if r2_filename.endswith('.gz'):
        f2 = gzip.open(r2_filename, 'rt')
    else:
        f2 = open(r2_filename, 'r')

    if output_filename.endswith('.gz'):
        f_out = gzip.open(output_filename, 'wt')
    else:
        f_out = open(output_filename, 'w')

    lineCount = 0
    id1 = f1.readline()
    while id1:
        lineCount += 1
        seq1 = f1.readline()
        seq1 = seq1.strip()
        plus1 = f1.readline()
        qual1 = f1.readline()
        qual1 = qual1.strip()

        id2 = f2.readline()
        seq2 = reverse_complement(f2.readline().strip()) + "\n"
        plus2 = f2.readline()
        qual2 = f2.readline()

        f_out.write(id1 + seq1 + seq2 + plus1 + qual1 + qual2)

        id1 = f1.readline()
    f1.close()
    f2.close()
    f_out.close()

    return (lineCount)


def split_interleaved_fastq(fastq_filename, output_filename_r1, output_filename_r2):
    """Split an interleaved fastq file into two files, one for each read pair.

    This assumes that the input fastq file is interleaved, i.e. that the reads are ordered as follows:
        R1
        R2
        R1
        R2
        ...

    And results in two files, one for each read pair:
        output_filename_r1
            R1
            R1
            ...
        output_filename_r2
            R2
            R2
            ...

    Parameters
    ----------
    fastq_filename : str
        Path to the input fastq file.
    output_filename_r1 : str
        Path to the output fastq file for r1.
    output_filename_r2 : str
        Path to the output fastq file for r2.

    Returns
    -------
    output_filename_r1 : str
        Path to the output fastq file for r1.
    output_filename_r2 : str
        Path to the output fastq file for r2.
    """
    if fastq_filename.endswith('.gz'):
        fastq_handle = gzip.open(fastq_filename, 'rt')
    else:
        fastq_handle = open(fastq_filename)

    try:
        fastq_splitted_outfile_r1 = gzip.open(output_filename_r1, 'wt')
        fastq_splitted_outfile_r2 = gzip.open(output_filename_r2, 'wt')
        [fastq_splitted_outfile_r1.write(line) if (i % 8 < 4) else fastq_splitted_outfile_r2.write(line) for i, line in enumerate(fastq_handle)]
    except:
        raise BadParameterException('Error in splitting read pairs from a single file')
    finally:
        fastq_handle.close()
        fastq_splitted_outfile_r1.close()
        fastq_splitted_outfile_r2.close()

    return output_filename_r1, output_filename_r2


######
# allele modification functions
######

def get_row_around_cut_asymmetrical(row,cut_point,plot_left,plot_right):
    cut_idx=row['ref_positions'].index(cut_point)
    return row['Aligned_Sequence'][cut_idx-plot_left+1:cut_idx+plot_right+1],row['Reference_Sequence'][cut_idx-plot_left+1:cut_idx+plot_right+1],row['Read_Status']=='UNMODIFIED',row['n_deleted'],row['n_inserted'],row['n_mutated'],row['#Reads'], row['%Reads']


def get_dataframe_around_cut_asymmetrical(df_alleles, cut_point,plot_left,plot_right,collapse_by_sequence=True):
    if df_alleles.shape[0] == 0:
        return df_alleles
    ref1 = df_alleles['Reference_Sequence'].iloc[0]
    ref1 = ref1.replace('-','')
    
    df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut_asymmetrical(row,cut_point,plot_left,plot_right),axis=1).values),
                    columns=['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated','#Reads','%Reads'])

    df_alleles_around_cut=df_alleles_around_cut.groupby(['Aligned_Sequence','Reference_Sequence','Unedited','n_deleted','n_inserted','n_mutated']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by=['#Reads', 'Aligned_Sequence', 'Reference_Sequence'], inplace=True, ascending=[False, True, True])
    df_alleles_around_cut['Unedited'] = df_alleles_around_cut['Unedited'] > 0
    return df_alleles_around_cut


def get_row_around_cut_debug(row, cut_point, offset):
    cut_idx = row['ref_positions'].index(cut_point)
    # don't check overflow -- it was checked when program started
    return row['Aligned_Sequence'][cut_idx - offset + 1:cut_idx + offset + 1], row['Reference_Sequence'][
                                                                               cut_idx - offset + 1:cut_idx + offset + 1], \
           row['Read_Status'] == 'UNMODIFIED', row['n_deleted'], row['n_inserted'], row['n_mutated'], row['#Reads'], \
           row['%Reads'], row['Aligned_Sequence'], row['Reference_Sequence']


def get_dataframe_around_cut_debug(df_alleles, cut_point, offset):
    df_alleles_around_cut = pd.DataFrame(
        list(df_alleles.apply(lambda row: get_row_around_cut_debug(row, cut_point, offset), axis=1).values),
        columns=['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated', '#Reads',
                 '%Reads', 'oSeq', 'oRef'])
    df_alleles_around_cut = df_alleles_around_cut.groupby(
        ['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated', 'oSeq',
         'oRef']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by=['#Reads', 'Aligned_Sequence', 'Reference_Sequence'], inplace=True, ascending=[False, True, True])
    df_alleles_around_cut['Unedited'] = df_alleles_around_cut['Unedited'] > 0
    return df_alleles_around_cut


def get_amplicon_info_for_guides(ref_seq, guides, guide_mismatches, guide_names, quantification_window_centers,
                                 quantification_window_sizes, quantification_window_coordinates, exclude_bp_from_left,
                                 exclude_bp_from_right, plot_window_size, guide_plot_cut_points,
                                 discard_guide_positions_overhanging_amplicon_edge=False):
    """
    gets cut site and other info for a reference sequence and a given list of guides

    input:
    ref_seq : reference sequence
    guides : a list of guide sequences
    guide_mismatches : a list of positions where a guide may have mismatches (for flexiguides)
    guide_names : a list of names for each guide
    quantification_window_centers : a list of positions where quantification is centered for each guide
    quantification_window_sizes : a list of lengths of quantification windows extending from quantification_window_center for each guide
    quantification_window_coordinates: if given, these override quantification_window_center and quantification_window_size for setting quantification window. These are specific for this amplicon.
    exclude_bp_from_left : these bp are excluded from the quantification window
    exclude_bp_from_right : these bp are excluded from the quantification window
    plot_window_size : length of window extending from quantification_window_center to plot
    guide_plot_cut_points : whether or not to add cut point to plot (prime editing flaps don't have cut points)
    discard_guide_positions_overhanging_amplicon_edge : if True, for guides that align to multiple positions, guide positions will be discarded if plotting around those regions would included bp that extend beyond the end of the amplicon.

    returns:
    this_sgRNA_sequences : list of sgRNAs that are in this amplicon
    this_sgRNA_intervals : indices of each guide
    this_sgRNA_cut_points : cut points for each guide (defined by quantification_window_center)
    this_sgRNA_plot_cut_points : whether or not a cut point is plotted
    this_sgRNA_plot_idxs : list of indices to be plotted for each sgRNA
    this_sgRNA_mismatches: list of mismatches between the guide and the amplicon
    this_sgRNA_names : list of names for each sgRNA (to disambiguate in case a sequence aligns to multiple positions)
    this_sgRNA_include_idxs : list of indices to be included in quantification per guide
    this_include_idxs : list of indices to be included in quantification
    this_exclude_idxs : list of indices to be excluded from quantification
    """
    ref_seq_length = len(ref_seq)

    this_sgRNA_sequences = []
    this_sgRNA_intervals = []
    this_sgRNA_cut_points = []
    this_sgRNA_plot_cut_points = []
    this_sgRNA_plot_idxs = []
    this_sgRNA_mismatches = []
    this_sgRNA_names = []
    this_sgRNA_include_idxs = []
    this_include_idxs = []
    this_exclude_idxs = []

    if exclude_bp_from_left:
        this_exclude_idxs += range(exclude_bp_from_left)

    if exclude_bp_from_right:
        this_exclude_idxs += range(ref_seq_length)[-exclude_bp_from_right:]

    window_around_cut = max(1, plot_window_size)

    seen_cut_points = {}  # keep track of cut points in case 2 gudes cut at same position (so they can get different names)
    seen_guide_names = {}  # keep track of guide names (so we don't assign a guide the same name as another guide)
    for guide_idx, current_guide_seq in enumerate(guides):
        if current_guide_seq == '':
            continue
        offset_fw = quantification_window_centers[guide_idx] + len(current_guide_seq) - 1
        offset_rc = (-quantification_window_centers[guide_idx]) - 1

        # .. run once with findall to get number of matches
        fw_regex = r'(?=(' + re.escape(current_guide_seq) + r'))'
        fw_matches = re.findall(fw_regex, ref_seq, flags=re.IGNORECASE)
        rv_regex = r'(?=(' + re.escape(reverse_complement(current_guide_seq)) + r'))'
        rv_matches = re.findall(rv_regex, ref_seq, flags=re.IGNORECASE)
        match_count = len(fw_matches) + len(rv_matches)

        # and now create the iter which will keep track of the locations of matches
        # (you can't get the length of an iter, and the findall only gives an array of matched strings and doesn't give locations of matches)
        fw_matches = re.finditer(fw_regex, ref_seq, flags=re.IGNORECASE)
        rv_matches = re.finditer(rv_regex, ref_seq, flags=re.IGNORECASE)

        # for every match, append:
        # this_sgRNA_cut_points, this_sgRNA_intervals,this_sgRNA_mismatches,this_sgRNA_names,this_sgRNA_sequences,this_sgRNA_include_idxs,include_idxs
        for m in fw_matches:
            cut_p = m.start() + offset_fw
            if discard_guide_positions_overhanging_amplicon_edge and (
                    (cut_p - window_around_cut + 1 < 0) or (cut_p + window_around_cut > ref_seq_length - 1)):
                continue
            this_sgRNA_cut_points.append(cut_p)
            this_sgRNA_plot_cut_points.append(guide_plot_cut_points[guide_idx])
            this_sgRNA_intervals.append((m.start(), m.start() + len(current_guide_seq) - 1))
            this_sgRNA_mismatches.append(guide_mismatches[guide_idx])

            if quantification_window_sizes[guide_idx] > 0:
                st = max(0, cut_p - quantification_window_sizes[guide_idx] + 1)
                en = min(ref_seq_length - 1, cut_p + quantification_window_sizes[guide_idx] + 1)
                this_include_idxs.extend(range(st, en))
                this_sgRNA_include_idxs.append(list(range(st, en)))
            else:
                this_sgRNA_include_idxs.append([])

            this_sgRNA_name = guide_names[guide_idx]
            if match_count == 1:
                this_sgRNA_names.append(this_sgRNA_name)
            else:
                if this_sgRNA_name == "":
                    this_sgRNA_name = current_guide_seq.upper()
                this_potential_name = this_sgRNA_name + "_" + str(m.start())
                curr_guide_idx = 1
                while this_potential_name in seen_guide_names:
                    curr_guide_idx += 1
                    this_potential_name = this_sgRNA_name + "_" + str(m.start()) + "_" + curr_guide_idx
                this_sgRNA_names.append(this_potential_name)
            this_sgRNA_sequences.append(current_guide_seq.upper())

        for m in rv_matches:
            cut_p = m.start() + offset_rc  # cut position
            if discard_guide_positions_overhanging_amplicon_edge and (
                    (cut_p - window_around_cut + 1 < 0) or (cut_p + window_around_cut > ref_seq_length - 1)):
                continue
            this_sgRNA_cut_points.append(cut_p)
            this_sgRNA_plot_cut_points.append(guide_plot_cut_points[guide_idx])
            this_sgRNA_intervals.append((m.start(), m.start() + len(current_guide_seq) - 1))
            this_sgRNA_mismatches.append([len(current_guide_seq) - (x + 1) for x in guide_mismatches[guide_idx]])

            if quantification_window_sizes[guide_idx] > 0:
                st = max(0, cut_p - quantification_window_sizes[guide_idx] + 1)
                en = min(ref_seq_length - 1, cut_p + quantification_window_sizes[guide_idx] + 1)
                this_include_idxs.extend(range(st, en))
                this_sgRNA_include_idxs.append(list(range(st, en)))
            else:
                this_sgRNA_include_idxs.append([])

            this_sgRNA_name = guide_names[guide_idx]
            if match_count == 1:
                this_sgRNA_names.append(this_sgRNA_name)
            else:
                if this_sgRNA_name == "":
                    this_sgRNA_name = current_guide_seq.upper()
                this_potential_name = this_sgRNA_name + "_" + str(m.start())
                curr_guide_idx = 1
                while this_potential_name in seen_guide_names:
                    curr_guide_idx += 1
                    this_potential_name = this_sgRNA_name + "_" + str(m.start()) + "_" + curr_guide_idx
                this_sgRNA_names.append(this_potential_name)
            this_sgRNA_sequences.append(current_guide_seq.upper())

    # create mask of positions in which to include/exclude indels for the quantification window
    # first, if exact coordinates have been given, set those
    given_include_idxs = []
    if quantification_window_coordinates is not None and quantification_window_coordinates != "0":
        coordinate_include_idxs = []
        theseCoords = str(quantification_window_coordinates).split("_")
        for coord in theseCoords:
            coordRE = re.match(r'^(\d+)-(\d+)$', coord)
            if coordRE:
                start = int(coordRE.group(1))
                end = int(coordRE.group(2)) + 1
                if end > ref_seq_length:
                    raise NTException("End coordinate " + str(end) + " for '" + str(coord) + "' in '" + str(
                        theseCoords) + "' is longer than the sequence length (" + str(ref_seq_length) + ")")
                coordinate_include_idxs.extend(range(start, end))
            else:
                raise NTException("Cannot parse analysis window coordinate '" + str(coord) + "' in '" + str(
                    theseCoords) + "'. Coordinates must be given in the form start-end e.g. 5-10 . Please check the --analysis_window_coordinate parameter.")
        given_include_idxs = coordinate_include_idxs
        this_include_idxs = coordinate_include_idxs
        this_sgRNA_include_idxs = [[]] * len(this_sgRNA_include_idxs)
    # otherwise, take quantification centers + windows calculated for each guide above
    elif this_sgRNA_cut_points and len(this_include_idxs) > 0:
        given_include_idxs = this_include_idxs
    else:
        this_include_idxs = range(ref_seq_length)
        this_sgRNA_include_idxs = [[]] * len(this_sgRNA_include_idxs)

    # flatten the arrays to avoid errors with old numpy library
    this_include_idxs = np.ravel(this_include_idxs)
    this_exclude_idxs = np.ravel(this_exclude_idxs)
    given_include_idxs = np.ravel(given_include_idxs)
    pre_exclude_include_idxs = this_include_idxs.copy()

    this_include_idxs = set(np.setdiff1d(this_include_idxs, this_exclude_idxs))
    if len(np.setdiff1d(given_include_idxs, list(this_include_idxs))) > 0:
        raise BadParameterException(
            'The quantification window has been partially exluded by the --exclude_bp_from_left or --exclude_bp_from_right parameters. Given: ' + str(
                given_include_idxs) + ' Pre: ' + str(pre_exclude_include_idxs) + ' Post: ' + str(this_include_idxs))
    if len(this_include_idxs) == 0:
        if len(pre_exclude_include_idxs) > 0:
            raise BadParameterException(
                'The quantification window around the sgRNA is excluded. Please decrease the exclude_bp_from_right and exclude_bp_from_left parameters.')
        else:
            raise BadParameterException(
                'The entire sequence has been excluded. Please enter a longer amplicon, or decrease the exclude_bp_from_right and exclude_bp_from_left parameters')

    if this_sgRNA_cut_points and plot_window_size > 0:
        for cut_p in this_sgRNA_cut_points:
            if cut_p - window_around_cut + 1 < 0:
                logging.warning('Offset around cut would extend to the left of the amplicon. Window will be truncated.')
                window_around_cut = cut_p + 1
            if cut_p + window_around_cut > ref_seq_length - 1:
                logging.warning('Offset around cut would be greater than reference sequence length. Window will be truncated.')
                window_around_cut = ref_seq_length - cut_p - 1
            st = max(0, cut_p - window_around_cut + 1)
            en = min(ref_seq_length - 1, cut_p + window_around_cut + 1)
            this_sgRNA_plot_idxs.append(sorted(list(range(st, en))))
    else:
        this_sgRNA_plot_idxs.append(range(ref_seq_length))

    this_include_idxs = np.sort(list(this_include_idxs))
    this_exclude_idxs = np.sort(list(this_exclude_idxs))

    return this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_sgRNA_names, this_sgRNA_include_idxs, this_include_idxs, this_exclude_idxs


def set_guide_array(vals, guides, property_name):
    """
    creates an int array of vals of the same length as guide, filling in missing values with the first item in vals
    vals: comma-separated string of values
    guides: list of guides
    property_name: property name, useful for creating warning to user

    returns: list of int vals
    """
    # if only one is given, set it for all guides
    vals_array = str(vals).split(",")
    ret_array = [int(vals_array[0])] * len(guides)
    # len(vals_array) is always one -- don't freak out if it's longer than 0 guides
    if len(vals_array) > 1 and len(vals_array) > len(guides):
        raise BadParameterException("More %s values were given than guides. Guides: %d %ss: %d" % (
        property_name, len(guides), property_name, len(vals_array)))

    if len(guides) == 0:
        return []

    for idx, val in enumerate(vals_array):
        if val != '':
            ret_array[idx] = int(val)
    return ret_array


def get_relative_coordinates(to_sequence, from_sequence):
    """Given an alignment, get the relative coordinates of the second sequence to the first.

    For example, from_sequence[i] matches to to_sequence[inds[i]]. A `-1`
    indicates a gap at the beginning of `to_sequence`.

    Parameters
    ----------
    to_sequence : str
        The alignment of the first sequence (where the coordinates are relative to)
    from_sequence : str
        The alignment of the second sequence

    Returns
    -------
    s1inds_gap_left : list of int
        The relative coordinates of the second sequence to the first, where gaps
        in the first sequence are filled with the left value.
    s1inds_gap_right : list of int
        The relative coordinates of the second sequence to the first, where gaps
        in the first sequence are filled with the right value.
    """
    s1inds_gap_left = []
    s1inds_gap_right = []
    s1idx_left = -1
    s1idx_right = 0
    s2idx = -1
    for ix in range(len(to_sequence)):
        if to_sequence[ix] != "-":
            s1idx_left += 1
        if from_sequence[ix] != "-":
            s2idx += 1
            s1inds_gap_left.append(s1idx_left)
            s1inds_gap_right.append(s1idx_right)
        if to_sequence[ix] != "-":
            s1idx_right += 1

    return s1inds_gap_left, s1inds_gap_right


def get_alignment_coordinates(to_sequence, from_sequence, aln_matrix, needleman_wunsch_gap_open,
                              needleman_wunsch_gap_extend):
    """
    gets the coordinates to align from from_sequence to to_sequence
    such that from_sequence[i] matches to to_sequence[inds[i]]
    params:
    to_sequence : sequence to align to
    from_sequence: sequence to align from
    aln_matrix: matrix specifying alignment substitution scores in the NCBI format
    needleman_wunsch_gap_open: needleman wunsch gap open penalty
    needleman_wunsch_gap_extend: needleman wunsch gap extend penalty

    returns:
    inds_l : if there's a gap in to_sequence, these values are filled with the left value
    inds_r : if there's a gap in to_sequence, these values are filled with the right value (after the gap)
    """
    this_gap_incentive = np.zeros(len(from_sequence) + 1, dtype=int)
    fws1, fws2, fwscore = CRISPResso2Align.global_align(to_sequence, from_sequence, matrix=aln_matrix,
                                                        gap_open=needleman_wunsch_gap_open,
                                                        gap_extend=needleman_wunsch_gap_extend,
                                                        gap_incentive=this_gap_incentive)
    return get_relative_coordinates(fws1, fws2)


def get_mismatches(seq_1, seq_2, aln_matrix, needleman_wunsch_gap_open, needleman_wunsch_gap_extend):
    """
    for a specific sequence seq_1, gets the location of mismatches as an array of length(seq_1) with reference to seq_2
    Only searches in the positive direction
    Only positions of seq_1 covered in seq_2 are included (e.g. if the base in seq_1 is deleted in seq_2 it won't be included)
    params:
    seq_1: sequence to compare
    seq_2: sequence to compare against
    aln_matrix: matrix specifying alignment substitution scores in the NCBI format. Can be read by: aln_matrix = CRISPResso2Align.read_matrix(aln_matrix_loc)
    needleman_wunsch_gap_open: needleman wunsch gap open penalty
    needleman_wunsch_gap_extend: needleman wunsch gap extend penalty

    returns:
    mismatches_arr: positions in seq_1 that mismatch seq_2
    """
    # when we set up the gap_flanking score, this should be set to the gap open penalty -- we'd like a good end-to-end alignment here
    coords_l, coords_r = get_alignment_coordinates(seq_2, seq_1, aln_matrix, needleman_wunsch_gap_open,
                                                   needleman_wunsch_gap_extend)
    mismatch_coords = []
    for i in range(len(seq_1)):
        if coords_l[i] < 0 or coords_l[i] >= len(seq_2) or seq_1[i] != seq_2[coords_l[i]] or coords_r[i] >= len(
                seq_2) or seq_1[i] != seq_2[coords_r[i]]:
            mismatch_coords.append(i)
    return mismatch_coords


def get_best_aln_pos_and_mismatches(guide_seq, within_amp_seq, aln_matrix, needleman_wunsch_gap_open,
                                    needleman_wunsch_gap_extend=0):
    """
    for a specific guide, gets the location, sequence, and mismatches for that guide at that location
    Only searches in the positive direction
    params:
    guide_seq: short guide sequence to look for
    within_amp_seq: longer amplicon to search within
    aln_matrix: matrix specifying alignment substitution scores in the NCBI format
    needleman_wunsch_gap_open: gap open penalty for needleman wunsch
    needleman_wunsch_gap_extend: gap extend pentaly for needleman wunsch (set as 0 by default to allow for longer stretches of gaps)

    returns:
    best_aln_seq: sequence in within_amp_seq of best hit
    best_aln_score: score of best alignment
    best_aln_mismatches: mismatches between best_aln_seq and guide_seq
    best_aln_start: start location of best match
    best_aln_end: end location of best match (1pb after last base)
    s1: aligned s1 (guide_seq)
    s2: aligned s2 (within_amp_seq)
    """
    ref_incentive = np.zeros(len(within_amp_seq) + 1, dtype=int)
    s1, s2, fw_score = CRISPResso2Align.global_align(guide_seq, within_amp_seq, matrix=aln_matrix,
                                                     gap_incentive=ref_incentive, gap_open=needleman_wunsch_gap_open,
                                                     gap_extend=0, )
    range_start_dashes = 0
    m = re.search(r'^-+', s1)
    if (m):
        range_start_dashes = m.span()[1]
    range_end_dashes = len(within_amp_seq)
    m = re.search(r'-+$', s1)
    if (m):
        range_end_dashes = m.span()[0]
    guide_seq_in_amp = s2[range_start_dashes:range_end_dashes].replace("-", "")

    best_aln_seq = guide_seq_in_amp
    best_aln_mismatches = get_mismatches(guide_seq, guide_seq_in_amp, aln_matrix, needleman_wunsch_gap_open,
                                         needleman_wunsch_gap_extend)
    best_aln_start = range_start_dashes
    best_aln_end = range_end_dashes
    return (best_aln_seq, fw_score, best_aln_mismatches, best_aln_start, best_aln_end, s1, s2)


def get_sgRNA_mismatch_vals(seq1, seq2, start_loc, end_loc, coords_l, coords_r, rev_coords_l, rev_coords_r):
    """
    given coordinates between one sequence and the other, given a start and end position, finds the locations of mismatches between the sequences between the start and end positions
    the sgRNA (as defined between start_loc and end_loc) is a substring of one of the samples, but the entire seqs are aligned to avoid any funny partial alignments
    seq1 : first sequence to compare
    seq2 : second sequence to compare
    start_loc : location to start comparison
    end_loc : location to stop comparison
    coords_l, coords_r : cooresponding indices between seq1 -> seq2
    rev_coords_l, rev_coords_r : cooresponding indices between seq2 -> seq1

    returns:
    mismatches between the two strings (for use as sgRNA_mismatches)
    """
    this_mismatches = []
    seen_inds = {}  # detect deletions in other string
    last_mismatch_val = rev_coords_l[coords_l[start_loc]]  # detect deletions in this string
    for i in range(coords_l[start_loc], coords_r[end_loc]):
        if seq1[i] != seq2[rev_coords_l[i]] or seq1[i] == "-":
            this_mismatches.append(i - coords_l[start_loc])
        this_last_mismatch_val = rev_coords_l[i]
        if this_last_mismatch_val in seen_inds:
            this_mismatches.append(i - coords_l[start_loc])
        elif this_last_mismatch_val != last_mismatch_val:  # if we skipped an index (deletion in other sequence)
            this_mismatches.append(i - coords_l[start_loc] - 0.5)
        seen_inds[this_last_mismatch_val] = 1
        last_mismatch_val = this_last_mismatch_val + 1
    return list(set(this_mismatches))


def get_quant_window_ranges_from_include_idxs(include_idxs):
    """Given a list of indexes, return the ranges that those indexes include.

    Parameters
    ----------
    include_idxs: list
       A list of indexes included in the quantification window, for example
       `[20, 21, 22, 35, 36, 37]`.

    Returns
    -------
    list
       A list of tuples representing the ranges included in the quantification
       window, for example `[(20, 22), (35, 37)]`. If there is a single index, it
       will be reported as `[(20, 20)]`.
    """
    quant_ranges = []
    if include_idxs is None or len(include_idxs) == 0:
        return quant_ranges
    start_idx = include_idxs[0]
    last_idx = include_idxs[0]
    for idx in include_idxs[1:]:
        if idx == last_idx + 1:
            last_idx = idx
        else:
            quant_ranges.append((start_idx, last_idx))
            start_idx = idx
            last_idx = idx
    quant_ranges.append((start_idx, last_idx))
    return quant_ranges


######
# terminal functions
######
def get_crispresso_logo():
    return (r'''
     _
    '  )
    .-'
   (____
C)|     \
  \     /
   \___/
''')


def get_crispresso_header(description, header_str):
    """
    Creates the CRISPResso header string with the header_str between two crispresso mugs
    """
    term_width = 80
    try:
        term_width = os.get_terminal_size().columns
    except:
        pass

    logo = get_crispresso_logo()
    logo_lines = logo.splitlines()
    max_logo_width = max([len(x) for x in logo_lines])

    output_line = ""
    if header_str is not None:
        header_str = header_str.rstrip()

        header_lines = header_str.splitlines()
        while len(header_lines) < len(logo_lines):
            header_lines = [""] + header_lines
        while len(header_lines) > len(logo_lines):
            logo_lines = [""] + logo_lines

        max_header_width = max([len(x) for x in header_lines])

        pad_space = int((term_width - (max_logo_width * 2) - max_header_width) / 4)
        pad_string = " " * pad_space

        for i in range(len(logo_lines))[::-1]:
            output_line = (logo_lines[i].ljust(max_logo_width) + pad_string + header_lines[i].ljust(
                max_header_width) + pad_string + logo_lines[i].ljust(max_logo_width)).center(
                term_width) + "\n" + output_line

    else:
        pad_space = int((term_width - max_logo_width) / 2)
        pad_string = " " * pad_space
        for i in range(len(logo_lines))[::-1]:
            output_line = (pad_string + logo_lines[i].ljust(max_logo_width) + pad_string).center(
                term_width) + "\n" + output_line

    output_line += '\n' + ('[CRISPResso version ' + __version__ + ']').center(term_width) + '\n' + (
        '[Note that as of version 2.3.0 FLASh and Trimmomatic have been replaced by fastp for read merging and trimming. Accordingly, the --flash_command and --trimmomatic_command parameters have been replaced with --fastp_command. Also, --trimmomatic_options_string has been replaced with --fastp_options_string.\n\nAlso in version 2.3.2, when running CRISPRessoPooled in mixed-mode (amplicon file and genome are provided) the default behavior will be as if the --demultiplex_only_at_amplicons parameter is provided. This change means that reads and amplicons do not need to align to the exact locations.]').center(
        term_width) + "\n" + ('[For support contact k.clement@utah.edu or support@edilytics.com]').center(term_width) + "\n"

    description_str = ""
    for str in description:
        str = str.strip()
        description_str += str.center(term_width) + "\n"

    return "\n" + description_str + output_line


def get_crispresso_footer():
    logo = get_crispresso_logo()
    logo_lines = logo.splitlines()

    max_logo_width = max([len(x) for x in logo_lines])
    pad_space = int((80 - max_logo_width) / 2)
    pad_string = " " * pad_space

    output_line = ""
    for i in range(len(logo_lines))[::-1]:
        output_line = pad_string + logo_lines[i].ljust(max_logo_width) + pad_string + "\n" + output_line

    return output_line

def format_cl_text(text, max_chars=None, spaces_to_tab=4):
    """
    Formats text for command line output
    params:
        text: text to format
        max_chars: maximum number of characters per line
        spaces_to_tab: number of spaces to add at the beginning of wrapped lines
    """
    if max_chars is None:
        try:
            max_chars = os.get_terminal_size().columns
        except:
            max_chars = 80

    broken_lines = []
    for line in text.split('\n'):
        if len(line) > max_chars:
            broken_line = ''
            while len(line) > max_chars:
                # Find the last space before max_chars characters
                last_space_index = line.rfind(' ', spaces_to_tab, max_chars)
                if last_space_index == -1:
                    # If no space found, break the line at max_chars characters
                    last_space_index = max_chars
                broken_line += line[:last_space_index] + '\n'
                line = ' ' * spaces_to_tab + line[last_space_index:].strip()
            broken_line += line
            broken_lines.append(broken_line)
        else:
            broken_lines.append(line)
    return '\n'.join(broken_lines)


def zip_results(results_folder):
    path_values = os.path.split(results_folder)
    output_folder = path_values[0]
    folder_id = path_values[1]
    if output_folder == "":
        cmd_to_zip = 'zip -m -r {0} {1}'.format(
            folder_id + ".zip", folder_id
        )
    else:
        cmd_to_zip = 'cd {0} && zip -m -r {1} {2} .'.format(
            output_folder, folder_id + ".zip", folder_id
        )
    sb.call(cmd_to_zip, shell=True)
    return


def is_C2Pro_installed():
    try:
        spec = importlib.util.find_spec("CRISPRessoPro")
        if spec is None:
            return False
        else:
            return True
    except:
        return False


def check_custom_config(args):
    """Check if the config_file argument was provided. If so load the configurations from the file, otherwise load default configurations.

    Parameters:
    -------------
    args : dict
        All arguments passed into the crispresso run.

    Returns:
    -------------
    config : dict
        A dict with a 'colors' key that contains hex color values for different report items as well as a 'guardrails' key that contains the guardrail values.
    """
    config =  {
        "colors": {
            'Substitution': '#0000FF',
            'Insertion': '#008000',
            'Deletion': '#FF0000',
            'A': '#7FC97F',
            'T': '#BEAED4',
            'C': '#FDC086',
            'G': '#FFFF99',
            'N': '#C8C8C8',
            '-': '#1E1E1E',
        },
        "guardrails": {
            'min_total_reads': 10000,
            'aligned_cutoff': 0.9,
            'alternate_alignment': 0.3,
            'min_ratio_of_mods_in_to_out': 0.01,
            'modifications_at_ends': 0.01,
            'outside_window_max_sub_rate': 0.002,
            'max_rate_of_subs': 0.3,
            'guide_len': 19,
            'amplicon_len': 50,
            'amplicon_to_read_length': 1.5
        }
    }

    logger = logging.getLogger(getmodule(stack()[1][0]).__name__)
    if not is_C2Pro_installed():
        return config

    if args.config_file:
        try:
            with open(args.config_file, "r") as json_file:
                custom_config = json.load(json_file)

            if 'guardrails' in custom_config.keys():
                for key in config['guardrails']:
                    if key not in custom_config['guardrails']:
                        logger.warn(f"Value for {key} not provided, defaulting.")
                        custom_config['guardrails'][key] = config['guardrails'][key]
                for key in custom_config['guardrails']:
                    if key not in config['guardrails']:
                        logger.warn(f"Key {key} is not a recognized guardrail parameter, skipping.")
            else:
                logger.warn("Json file does not contain the guardrails key. Defaulting all values.")
                custom_config['guardrails'] = config['guardrails']

            if 'colors' in custom_config.keys():
                for key in config['colors']:
                    if key not in custom_config['colors']:
                        logger.warn(f"Value for {key} not provided, defaulting")
                        custom_config['colors'][key] = config['colors'][key]
            else:
                logger.warn("Json file does not contain the colors key. Defaulting all values.")
                custom_config['colors'] = config['colors']

            return custom_config
        except Exception:
            if args.config_file:
                logger.warn("Cannot read config file '%s', defaulting config parameters." % args.config_file)
            else:
                logger.warn("No config file provided, defaulting config parameters.")
    return config


def safety_check(crispresso2_info, aln_stats, guardrails):
    """Check the results of analysis for potential issues and warns the user.

    Parameters
    ----------
    crispresso2_info : dict
        Dictionary of values describing the analysis
    aln_stats : dict
        Dictionary of alignment statistics and modification frequency and location
    guardrails : dict
        Contains the following:
            min_total_reads : int
                Cutoff value for total reads aligned
            alignedCutoff : float
                Check value for the percentage of reads aligned vs not aligned
            alternateAlignment : float
                Percentage of variance from expected reads to trigger guardrail
            minRatioOfModsInToOut : float
                Float representing the acceptable ratio of modifications in the window to out
            modificationsAtEnds : float
                The ratio of reads with modifications at the 0 or -1 spot
            outsideWindowMaxSubRate : float
                Ratio of subs allowed outside of the window
            maxRateOfSubs : float
                Allowed rate of subs accross the entire read
            guide_len : int
                Minimum guide length
            amplicon_len : int
                Minimum amplicon length
            ampliconToReadLen : float
                Comparison value between amplicons and reads
    """
    logger = logging.getLogger(getmodule(stack()[1][0]).__name__)
    messageHandler = GuardrailMessageHandler(logger)

    # Get amplicon and guide sequences and lengths
    amplicons = {}
    guide_groups = set()
    for name in crispresso2_info['results']['ref_names']:
        amplicons[name] = crispresso2_info['results']['refs'][name]['sequence_length']
        guide_groups.update(crispresso2_info['results']['refs'][name]['sgRNA_sequences'])
    unique_guides = {guide: len(guide) for guide in guide_groups}

    totalReadsGuardrail = TotalReadsGuardrail(messageHandler, guardrails['min_total_reads'])
    totalReadsGuardrail.safety(aln_stats['N_TOT_READS'])

    overallReadsAlignedGuard = OverallReadsAlignedGuardrail(messageHandler, guardrails['aligned_cutoff'])
    overallReadsAlignedGuard.safety(aln_stats['N_TOT_READS'], (aln_stats['N_CACHED_ALN'] + aln_stats['N_COMPUTED_ALN']))

    disproportionateReadsAlignedGuardrail = DisproportionateReadsAlignedGuardrail(messageHandler, guardrails['aligned_cutoff'])
    disproportionateReadsAlignedGuardrail.safety(aln_stats['N_TOT_READS'], crispresso2_info['results']['alignment_stats']['counts_total'])

    lowRatioOfModsInWindowToOut = LowRatioOfModsInWindowToOutGuardrail(messageHandler, guardrails['min_ratio_of_mods_in_to_out'])
    lowRatioOfModsInWindowToOut.safety(aln_stats['N_MODS_IN_WINDOW'], aln_stats['N_MODS_OUTSIDE_WINDOW'])

    highRateOfModificationAtEndsGuardrail = HighRateOfModificationAtEndsGuardrail(messageHandler, guardrails['modifications_at_ends'])
    highRateOfModificationAtEndsGuardrail.safety((aln_stats['N_CACHED_ALN'] + aln_stats['N_COMPUTED_ALN']), aln_stats['N_READS_IRREGULAR_ENDS'])

    highRateOfSubstitutionsOutsideWindowGuardrail = HighRateOfSubstitutionsOutsideWindowGuardrail(messageHandler, guardrails['outside_window_max_sub_rate'])
    highRateOfSubstitutionsOutsideWindowGuardrail.safety(aln_stats['N_GLOBAL_SUBS'], aln_stats['N_SUBS_OUTSIDE_WINDOW'])

    highRateOfSubstitutions = HighRateOfSubstitutionsGuardrail(messageHandler, guardrails['max_rate_of_subs'])
    highRateOfSubstitutions.safety(aln_stats['N_MODS_IN_WINDOW'], aln_stats['N_MODS_OUTSIDE_WINDOW'], aln_stats['N_GLOBAL_SUBS'])

    shortAmpliconSequence = ShortSequenceGuardrail(messageHandler, guardrails['amplicon_len'], 'amplicon')
    shortAmpliconSequence.safety(amplicons)

    shortGuideSequence = ShortSequenceGuardrail(messageHandler, guardrails['guide_len'], 'guide')
    shortGuideSequence.safety(unique_guides)

    longAmpliconShortReadsGuardrail = LongAmpliconShortReadsGuardrail(messageHandler, guardrails['amplicon_to_read_length'])
    longAmpliconShortReadsGuardrail.safety(amplicons, aln_stats['READ_LENGTH'])

    crispresso2_info['results']['guardrails'] = messageHandler.get_messages()

    return messageHandler.get_html_messages()


class GuardrailMessageHandler:
    """Class to handle message storage and display for guardrails"""
    def __init__(self, logger):
        """Create the message handler with an empty message array to collect html divs for the report

        Parameters:
        ------------
        logger : logger
            logger object used to display messages
        """
        self.logger = logger
        self.html_messages = []
        self.messages = {}

    def display_warning(self, guardrail, message):
        """Send the message to the logger to be displayed

        Parameters:
        -----------
        message : string
            Related guardrail message
        """
        self.logger.warning(message)
        self.messages[guardrail] = message

    def report_warning(self, message):
        """Create and store the html message to display on the report

        Parameters:
        -----------
        message : string
            Related guardrail message
        """
        html_warning = '<div class="alert alert-danger"><strong>Guardrail Warning!</strong>{0}</div>'.format(message)
        self.html_messages.append(html_warning)

    def get_messages(self):
        """Return the messages accumulated by the message handler"""
        return self.messages

    def get_html_messages(self):
        """Return the html messages accumulated by the message handler"""
        return self.html_messages


class TotalReadsGuardrail:
    """Guardrail class: check that the number of reads are above a minimum"""
    def __init__(self, messageHandler, minimum):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        minimum : int
            The comparison integer to determine if there are sufficent reads
        """
        self.messageHandler = messageHandler
        self.minimum = minimum
        self.message = " Low number of total reads: <{}.".format(minimum)

    def safety(self, total_reads):
        """Safety check, if total is below minimum send warnings

        Parameters:
        -----------
        total_reads : int
            The total reads, unaligned and aligned
        """
        if total_reads < self.minimum:
            self.message = self.message + " Total reads: {}.".format(total_reads)
            self.messageHandler.display_warning('TotalReadsGuardrail', self.message)
            self.messageHandler.report_warning(self.message)


class OverallReadsAlignedGuardrail:
    """Guardrail class: check if enough reads are aligned"""
    def __init__(self, messageHandler, cutoff):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : float
            The float representation of percentage of minimum reads to be aligned
        """
        self.messageHandler = messageHandler
        self.message = " <={val}% of reads were aligned.".format(val=(cutoff * 100))
        self.cutoff = cutoff

    def safety(self, total_reads, n_read_aligned):
        """Safety check, if total_reads divided by n_reads_aligned is lower than the cutoff

        Parameters:
        -----------
        total_reads : int
            Total reads, unaligned and aligned
        n_read_aligned : int
            Total aligned reads
        """
        if total_reads == 0:
            return
        if (n_read_aligned/total_reads) <= self.cutoff:
            self.message = self.message + " Total reads: {}, Aligned reads: {}.".format(total_reads, n_read_aligned)
            self.messageHandler.display_warning('OverallReadsAlignedGuardrail', self.message)
            self.messageHandler.report_warning(self.message)


class DisproportionateReadsAlignedGuardrail:
    """Guardrail class: check if the distribution of reads is roughly even across amplicons"""
    def __init__(self, messageHandler, cutoff):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : float
            The float representation of the accepted percentage deviation from the expected distribution
        """
        self.messageHandler = messageHandler
        self.message = " Disproportionate percentages of reads were aligned to amplicon: "
        self.cutoff = cutoff

    def safety(self, n_read_aligned, reads_aln_amplicon):
        """Safety check, if total_reads divided by n_reads_aligned is higher or lower than the total_reads divided by the number of amplicons

        Parameters:
        -----------
        total_reads : int
            Total reads, unaligned and aligned
        reads_aln_amplicon : dict
            A dictionary with the names of amplicons as the key and the number of reads aligned as the value
        """
        if len(reads_aln_amplicon.keys()) <= 1:
            return
        expected_per_amplicon = n_read_aligned / len(reads_aln_amplicon.keys())
        for amplicon, aligned in reads_aln_amplicon.items():
            if aligned <= (expected_per_amplicon * self.cutoff) or aligned >= (expected_per_amplicon * (1 - self.cutoff)):
                amplicon_message = self.message + amplicon + ", Percent of aligned reads aligned to this amplicon: {}%.".format(round((aligned/n_read_aligned) * 100, 2))
                self.messageHandler.display_warning('DisproportionateReadsAlignedGuardrail', amplicon_message)
                self.messageHandler.report_warning(amplicon_message)


class LowRatioOfModsInWindowToOutGuardrail:
    """Guardrail class: check the ratio of modifications in the quantification window to out of it"""
    def __init__(self, messageHandler, cutoff):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : float
            The float representation of the maximum percentage of modifications outside of the quantification window
        """
        self.messageHandler = messageHandler
        self.message = " <={}% of modifications were inside of the quantification window.".format(cutoff * 100)
        self.cutoff = cutoff

    def safety(self, mods_in_window, mods_outside_window):
        """Safety check, if the modifications in the window are below a ratio of the total

        Parameters:
        -----------
        mods_in_window : int
            The number of mods in the quantification window
        mods_outside_window : int
            The number of mods outside of the quantification window
        """
        total_mods = mods_in_window + mods_outside_window
        if total_mods == 0:
            return
        if ((mods_in_window / total_mods) <= self.cutoff):
            self.message = self.message + " Total modifications: {}, Modifications in window: {}, Modifications outside window: {}.".format(total_mods, mods_in_window, mods_outside_window)
            self.messageHandler.display_warning('LowRatioOfModsInWindowToOutGuardrail', self.message)
            self.messageHandler.report_warning(self.message)


class HighRateOfModificationAtEndsGuardrail:
    """Guardrail class: check the ratio of modifications in the quantification window to out of it"""
    def __init__(self, messageHandler, percentage_start_end):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        percentage_start_end : float
            The float representation of the maximum percentage reads that have modifications on either end
        """
        self.messageHandler = messageHandler
        self.message = " >={}% of reads have modifications at the start or end.".format(round(percentage_start_end * 100, 2))
        self.percent = percentage_start_end

    def safety(self, total_reads, irregular_reads):
        """Safety check, comparison between the number of irregular reads to total reads

        Parameters:
        -----------
        total_reads : int
            The number of mods in the quantification window
        irregular_reads : int
            The number of mods outside of the quantification window
        """
        if total_reads == 0:
            return
        if (irregular_reads / total_reads) >= self.percent:
            self.message = self.message + " Total reads: {}, Irregular reads: {}.".format(total_reads, irregular_reads)
            self.messageHandler.display_warning('HighRateOfModificationAtEndsGuardrail', self.message)
            self.messageHandler.report_warning(self.message)


class HighRateOfSubstitutionsOutsideWindowGuardrail:
    """Guardrail class: check the ratio of global substitutions to substitutions outside of the quantification window"""
    def __init__(self, messageHandler, cutoff):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : float
            The float representation of how many of the total substitutions can be outside of the quantification window
        """
        self.messageHandler = messageHandler
        self.message = " >={}% of substitutions were outside of the quantification window.".format(cutoff * 100)
        self.cutoff = cutoff

    def safety(self, global_subs, subs_outside_window):
        """Safety check, comparison between the number of global subs to subs outside of the quantification window

        Parameters:
        -----------
        global_subs : int
            The number of mods in the quantification window
        subs_outside_window : int
            The number of mods outside of the quantification window
        """
        if global_subs == 0:
            return
        if ((subs_outside_window / global_subs) >= self.cutoff):
            self.message = self.message + " Total substitutions: {}, Substitutions outside window: {}.".format(global_subs, subs_outside_window)
            self.messageHandler.display_warning('HighRateOfSubstitutionsOutsideWindowGuardrail', self.message)
            self.messageHandler.report_warning(self.message)


class HighRateOfSubstitutionsGuardrail:
    """Guardrail class: check the ratio of global substitutions to total modifications"""
    def __init__(self, messageHandler, cutoff):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : float
            The float representation of how many of the total modifications can be subsitutions
        """
        self.messageHandler = messageHandler
        self.message = " >={}% of modifications were substitutions. This could potentially indicate poor sequencing quality.".format(cutoff * 100)
        self.cutoff = cutoff

    def safety(self, mods_in_window, mods_outside_window, global_subs):
        """Safety check, comparison between subsitutions and total modifications

        Parameters:
        -----------
        mods_in_window : int
            Modifications inside of the quantification window
        mods_outside_window : int
            Modifications outside of the quantification window
        global_subs : int
            Total subsitutions across all reads
        """
        total_mods = mods_in_window + mods_outside_window
        if total_mods == 0:
            return
        if ((global_subs / total_mods) >= self.cutoff):
            self.message = self.message + " Total modifications: {}, Substitutions: {}.".format(int(total_mods), global_subs)
            self.messageHandler.display_warning('HighRateOfSubstitutionsGuardrail', self.message)
            self.messageHandler.report_warning(self.message)


class ShortSequenceGuardrail:
    """Guardrail class: Check to make sure that sequences (amplicons and guides) are above a certain length"""
    def __init__(self, messageHandler, cutoff, sequence_type):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : int
            Integer value to measure the length of the sequence
        sequence_type : string
            Amplicon or guide
        """
        self.messageHandler = messageHandler
        self.cutoff = cutoff
        self.message = f" {string.capwords(sequence_type)} length <{cutoff}: "

    def safety(self, sequences):
        """Safety check, comparison between sequence lengths and minimum lengths

        Parameters:
        -----------
        sequences : dict
            Dictionary with the name of the sequence as the key and the length of the sequence as the value
        """
        for name, length in sequences.items():
            if length < self.cutoff:
                sequence_message = self.message + name + ", Length: {}.".format(length)
                self.messageHandler.display_warning('ShortSequenceGuardrail', sequence_message)
                self.messageHandler.report_warning(sequence_message)


class LongAmpliconShortReadsGuardrail:
    """Guardrail class: Check to make sure that the reads are close in size to amplicons"""
    def __init__(self, messageHandler, cutoff):
        """Assign variables and create guardrail message

        Parameters:
        -----------
        messageHandler : GuardrailMessagehandler
            Guardrail message handler to create and display warnings
        cutoff : float
            The value multiplied by the read length to make sure the amplicon isn't too much longer than the reads.

        """
        self.messageHandler = messageHandler
        self.cutoff = cutoff
        self.message = " Amplicon length is greater than {}x the length of the reads: ".format(cutoff)

    def safety(self, amplicons, read_len):
        """Safety check, comparison between amplicon length and read length

        Parameters:
        -----------
        amplicons : dict
            Dictionary with the name of the amplicon as the key and the length of the sequence as the value
        read_len : int
            Average length of reads
        """
        for name, length in amplicons.items():
            if length > (read_len * self.cutoff):
                sequence_message = self.message + name + ", Amplicon length: {}, Read length: {}.".format(length, read_len)
                self.messageHandler.display_warning('LongAmpliconShortReadsGuardrail', sequence_message)
                self.messageHandler.report_warning(sequence_message)
