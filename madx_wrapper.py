"""
:module: madx.madx_wrapper

Runs MADX with a file or a string as an input, it processes @required macros.
If defined, writes the processed MADX script and logging output into files.
TODO: write tests
TODO: Possibly more anotation from MADX...
"""
from os.path import abspath, join, dirname
import sys
import re
import subprocess
import optparse

LIB = abspath(join(dirname(__file__), "madx", "lib"))
if "darwin" in sys.platform:
    MADX_PATH = abspath(join(dirname(__file__), "madx", "bin", "madx-macosx64-intel"))
elif "win" in sys.platform:
    MADX_PATH = abspath(join(dirname(__file__), "madx", "bin", "madx-win64-gnu.exe"))
else:
    MADX_PATH = abspath(join(dirname(__file__), "madx", "bin", "madx-linux64-gnu"))


def _parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file", metavar="FILE", dest="file",
                      help="The file with the annotated MADX input to run.")
    parser.add_option("-o", "--output", metavar="OUTPUT", dest="output",
                      help="If defined, it will write the processed MADX script into the file.")
    parser.add_option("-l", "--log", metavar="LOG", dest="log",
                      help="File where to write the MADX log output.")
    parser.add_option("-m", "--madx", metavar="MADX", dest="madx_path",
                      help="Path to the MAD-X executable to use", default=MADX_PATH)
    (options, args) = parser.parse_args()
    if len(args) > 1 or ((options.file is None) == (len(args) == 0)):
        raise IOError("No input found: either pass the file as first parameter or use --file (-f)")
    if len(args) == 1:
        return args[0], options.output, options.log, options.madx_path
    return options.file, options.output, options.log, options.madx_path


def resolve_and_run_file(input_file, output_file=None, log_file=None, madx_path=MADX_PATH):
    """Runs MADX in a subprocess.

    Attributes:
        input_file: MADX input file
        output_file: If given writes resolved MADX script.
        log_file: If given writes MADX logging output.
        madx_path: Path to MADX executable
    """
    input_string = _read_input_file(input_file)
    return resolve_and_run_string(input_string, output_file=output_file, log_file=log_file,
                                  madx_path=madx_path)


def resolve_and_run_string(input_string, output_file=None, log_file=None, madx_path=MADX_PATH):
    """Runs MADX in a subprocess.

    Attributes:
        input_string: MADX input string
        output_file: If given writes resolved MADX script.
        log_file: If given writes MADX logging output.
        madx_path: Path to MADX executable
    """
    _check_log_and_output_files(output_file, log_file)
    full_madx_script = _resolve(input_string, output_file)
    return _run(full_madx_script, log_file, madx_path)


def _resolve(input_string, output_file=None):
    """Resolves the !@requires annotations of the input_string, and returns the resulting script."""
    macro_calls = "option, -echo;\n" + _resolve_required_macros(input_string) + "option, echo;\n\n"
    full_madx_script = macro_calls + input_string
    if output_file is not None:
        with open(output_file, "w") as output:
            output.write(full_madx_script)
    return full_madx_script


def _run(full_madx_script, log_file=None, madx_path=MADX_PATH):
    if log_file is None:
        ret_value = _run_for_input_string(full_madx_script, madx_path)
    else:
        with open(log_file, "w") as output:
            ret_value = _run_for_input_string(full_madx_script, madx_path, logs=output)
    return ret_value


def _run_for_input_string(input_string, madx, logs=sys.stdout):
    process = subprocess.Popen(madx, shell=False, stdin=subprocess.PIPE, stdout=logs, stderr=logs)
    process.communicate(input_string.encode("utf-8"))
    return process.wait()


def _resolve_required_macros(file_content):
    """
    Recursively searches for "!@require lib" MADX annotations in the input script,
    adding the required macros library calls to the script header.
    """
    call_commands = ""
    for line in file_content.split("\n"):
        match = re.search("^!@require\s+([^\s]+).*$", line)
        if match is not None:
            required_macros = _add_macro_lib_ending(match.group(1))
            required_macros_file = abspath(join(LIB, required_macros))
            call_commands += _resolve_required_macros(_read_input_file(required_macros_file))
            call_commands += "call, file = \"" + required_macros_file + "\";\n"
    return call_commands


def _add_macro_lib_ending(macro_lib_name):
    macro_lib_name = macro_lib_name.strip()
    if macro_lib_name.endswith(".macros.madx"):
        return macro_lib_name
    else:
        return macro_lib_name + ".macros.madx"


def _read_input_file(input_file):
    with open(input_file) as text_file:
        return text_file.read()


def _check_log_and_output_files(output_file, log_file):
    if output_file is not None:
        open(output_file, "a").close()
    if log_file is not None:
        open(log_file, "a").close()


if __name__ == "__main__":
    resolve_and_run_file(*_parse_args())
