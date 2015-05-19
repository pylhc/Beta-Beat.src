import __init__  # @UnusedImport
import os
import sys
import re
import madxrunner
import optparse
from Utilities import iotools

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
LIB = "lib"


def _parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-o", "--output",
                    help="If defined, it will write the processed MADX script into the file.",
                    metavar="OUTPUT", dest="output")
    parser.add_option("-l", "--log",
                    help="File where to write the MADX output.",
                    metavar="LOG", dest="log")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit(-1)
    file_to_load = args[0]
    output_file = options.output
    log_file = options.log
    return file_to_load, output_file, log_file


def resolve_and_run_file(input_file, output_file=None, log_file=None):
    _check_files(input_file, output_file, log_file)
    main_file_content = iotools.read_all_lines_in_textfile(input_file)
    full_madx_script = resolve(main_file_content, output_file, log_file)
    run(full_madx_script, log_file)


def resolve_and_run_string(input_string, output_file=None, log_file=None):
    _check_files(None, output_file, log_file)
    full_madx_script = resolve(input_string, output_file, log_file)
    run(full_madx_script, log_file)


def resolve(input_string, output_file=None, log_file=None):
    """
    Resolves the !@requires annotations of the input_string, and returns
    the resulting script. If a output_file is given, the
    result is also written into the file. The MADX log will be
    written to log_file if given and to sys.stdout if not.
    """
    required_calls = "option, -echo;\n" + \
                     _resolve_required_macros(input_string) + \
                     "option, echo;\n\n"
    full_madx_script = required_calls + input_string
    if not output_file is None:
        with open(output_file, "w") as output_stream:
            output_stream.write(full_madx_script)
    return full_madx_script


def run(full_madx_script, log_file=None):
    """
    Runs MADX with the given, already resolved, MADX script.
    If log_file is given the MADX log will be written there.
    """
    if log_file is None:
        log_stream = sys.stdout
    else:
        log_stream = open(log_file, "w")
    madxrunner.runForInputString(full_madx_script, stdout=log_stream, stderr=log_stream)


def _resolve_required_macros(file_content):
    """
    Recursively searches for "!@require lib" annotations in the input script,
    adding the required macros library calls to the script header.
    """
    call_commands = ""
    for line in file_content.split("\n"):
        match = re.search("^!@require\s+([^\s]+).*$", line)
        if not match is None:
            required_macros = match.group(1)
            required_macros_file = os.path.abspath(os.path.join(CURRENT_PATH, LIB, required_macros))
            if not os.path.exists(required_macros_file):
                print >> sys.stderr, "Trying to import non existing macros library \"" + required_macros + "\". Aborting."
                sys.exit(-1)
            call_commands += _resolve_required_macros(iotools.read_all_lines_in_textfile(required_macros_file))
            call_commands += "call, file = \"" + required_macros_file + "\";\n"
    return call_commands


def _check_files(input_file, output_file, log_file):
    if not input_file is None:
        if not os.path.exists(input_file) or not os.path.isfile(input_file):
            print >> sys.stderr, "Cannot read input file: " + input_file + ". Aborting."
            sys.exit(-1)
    if not output_file is None:
        try:
            open(output_file, "a").close()
        except:
            print >> sys.stderr, "Cannot write output file: " + output_file + ". Aborting."
            sys.exit(-1)
    if not log_file is None:
        try:
            open(log_file, "a").close()
        except:
            print >> sys.stderr, "Cannot write log file: " + log_file + ". Aborting."
            sys.exit(-1)


if __name__ == "__main__":
    _file_to_load, _output_file, _log_file = _parse_args()
    resolve_and_run_file(_file_to_load, _output_file, _log_file)
