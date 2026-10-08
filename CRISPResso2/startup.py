"""Pandas-free helpers used before CORE is imported.

CLI startup prints the header and can serve ``--help`` / ``--version`` from
this module so a cold start does not wait on numpy/pandas/matplotlib.
"""

import argparse
import importlib.metadata
import json
import os
import sys
import textwrap
from pathlib import Path


class _StartupState:
    header_printed = False


def read_version():
    """Return the installed CRISPResso2 version string."""
    return importlib.metadata.version('CRISPResso2')


__version__ = read_version()


def header_already_printed():
    """Return whether the startup header has already been written to stdout."""
    return _StartupState.header_printed


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
    """Creates the CRISPResso header string with the header_str between two crispresso mugs
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


def print_startup_header(description, header_str=None):
    """Write the ASCII header to stdout immediately and mark it as printed."""
    sys.stdout.write(get_crispresso_header(description, header_str))
    sys.stdout.flush()
    _StartupState.header_printed = True


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

    json_path = Path(__file__).parent / 'args.json'

    with open(json_path, 'r') as json_file:
        args_dict = json.load(json_file)
        args_dict = args_dict["CRISPResso_args"]
    type_mapper = {
        "str": str,
        "int": int,
        "float": float,
    }

    def _add_args_from_dict(args_dict, tool, parser, type_mapper):
        """Add arguments from a parsed JSON args dict to the parser."""
        for key, value in args_dict.items():
            tools = value.get('tools', [])
            if tool in tools:
                action = value.get('action')
                required = value.get('required', False)
                default = value.get('default')
                type_value = value.get('type', 'str')
                arg_help = value.get('help', '') if value.get('help') != "SUPPRESS" else argparse.SUPPRESS

                if action:
                    parser.add_argument(*value['keys'], help=arg_help, action=action)
                elif required and default is None:
                    parser.add_argument(*value['keys'], help=arg_help, type=type_mapper[type_value], required=True)
                elif required:
                    parser.add_argument(*value['keys'], help=arg_help, default=default, type=type_mapper[type_value], required=True)
                else:
                    kwargs = {'help': arg_help, 'type': type_mapper[type_value]}
                    if default is not None: kwargs['default'] = default
                    parser.add_argument(*value['keys'], **kwargs)

    _add_args_from_dict(args_dict, tool, parser, type_mapper)

    # Load CRISPRessoPro args if installed. Prefer find_spec so --help does not
    # import Pro (and its heavy deps) just to read args.json.
    try:
        import importlib.util
        spec = importlib.util.find_spec('CRISPRessoPro')
        origin = spec.origin if spec is not None else None
        if origin:
            pro_args_path = Path(origin).parent / 'args.json'
            if pro_args_path.exists():
                with open(pro_args_path, 'r') as f:
                    pro_args_dict = json.load(f)
                pro_args_dict = pro_args_dict.get("CRISPRessoPro_args", {})
                _add_args_from_dict(pro_args_dict, tool, parser, type_mapper)
    except (ImportError, ModuleNotFoundError, ValueError):
        pass

    return parser
