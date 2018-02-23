"""
Module used to run either a MADX macro contained in one lib or a full template.

To run a macro use as:
madx_commander.py macro 'do_track_single_particle(0.00035, 0.0002, 500, "track")' -m "lhc_runI, tracking"

To run a template use as:
madx_commander.py template 'lhc_best_knowledge_madx(RUN, PATH, NUM_BEAM, MATCHER, ENERGY, QMX, QMY, DPP, STOP, QX, QY, QDX, QDY)'
See madx_templates_runner.py.

"""
import sys
import optparse
import madx_wrapper
from madx_templates_runner import MadxTemplates


def _parse_args():
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-m", "--macros",
                    help="List of comma-separated macros libraries to load. The '.macros.madx' ending is assumed",
                    metavar="MACROS", dest="macros")
    parser.add_option("-o", "--output",
                    help="If defined, it will write the processed MADX script into the file.",
                    metavar="OUTPUT", dest="output")
    parser.add_option("-l", "--log",
                    help="File where to write the MADX output.",
                    metavar="LOG", dest="log")
    options, args = parser.parse_args()
    if len(args) != 2:
        print parser.print_help()
        sys.exit(-1)
    command = args[0]
    call = args[1]
    output_file = options.output
    log_file = options.log
    if command == "macro":
        macros_list = _solve_macros_libs(options.macros)
        run_madx_command(call, macros_list, output_file, log_file)
    elif command == "template":
        if not options.macros is None:
            print >> sys.stderr, "The option -m isn't compatible with template mode."
            print parser.print_help()
            sys.exit(-1)
        _run_madx_template(call.strip("'").strip('"').strip())


def _solve_macros_libs(macros):
    macros_list = []
    raw_macros_list = macros.strip('"').strip("'").split(",")
    for macro_lib in raw_macros_list:
        macro_lib = madx_wrapper._add_macro_lib_ending(macro_lib)
    return macros_list


def _run_madx_template(template_call):
    madxTemplates = MadxTemplates()  # @UnusedVariable Used inside the string
    exec "madxTemplates." + template_call + ";"


def run_madx_command(command, macros_list, output_file, log_file):
    madx_script = create_madx_script(command, macros_list)
    madx_wrapper.resolve_and_run_string(madx_script, output_file, log_file)


def create_madx_script(command, macros_list):
    madx_script = _generate_requires(macros_list)
    madx_script += _generate_macro_call(command)
    return madx_script


def _generate_requires(macros_list):
    requires = ""
    for macro_lib in macros_list:
        requires += "!@require " + macro_lib + "\n"
    return requires


def _generate_macro_call(command):
    return "exec, " + command.strip('"').strip("'") + ";\n"


if __name__ == "__main__":
    raise DeprecationWarning("This will be removed soon, please consider using madx_wrapper.py")
    _parse_args()
