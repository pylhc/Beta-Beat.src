from optparse import OptionParser
import os
import sys

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


def _parse_args():
    parser = OptionParser()
    parser.add_option("-m", "--model",
                    help="Model twiss file to use.",
                    metavar="MODEL", dest="model_twiss")
    parser.add_option("-p", "--path",
                    help="Path to the modifiers file.",
                    metavar="MODIFIERS", dest="path")
    parser.add_option("--tunex",
                    help="Horizontal tune.",
                    metavar="TUNEX", dest="tunex")
    parser.add_option("--tuney",
                    help="Vertical tune.",
                    metavar="TUNEY", dest="tuney")
    parser.add_option("-e", "--energy",
                    help="The energy of the beam. It must be: 0.45TeV, 3.5TeV, 4TeV or 6.5TeV.",
                    metavar="ENERGY", dest="energy")
    options, _ = parser.parse_args()
    if not os.path.isfile(options.model_twiss):
        print >> sys.stderr, "Cannot find the specified twiss.dat file."
        sys.exit(-1)
    modifiers_path = options.path
    if modifiers_path is None:
        modifiers_path = os.path.dirname(options.model_twiss)
    if options.energy not in ["0.45TeV", "3.5TeV", "4TeV", "6.5TeV"]:
        print >> sys.stderr, "No valid energy specified, it must be one of: 0.45TeV, 3.5TeV, 4TeV or 6.5TeV"
        sys.exit(-1)
    return options.model_twiss, modifiers_path, options.tunex, options.tuney, options.energy


def process_masks(model_twiss, modifiers_path, tunex, tuney, energy):
    model_path = os.path.dirname(model_twiss)

    with open(os.path.join(CURRENT_PATH, "job.systematic.mask"), "r") as error_mask_lines:
        with open(os.path.join(model_path, "job.systematic.mask"), "w") as error_mask_output:
            for line in error_mask_lines:
                new_line = line
                new_line = new_line.replace("%PATH", modifiers_path)
                new_line = new_line.replace("%QMX", tunex)
                new_line = new_line.replace("%QMY", tuney)
                new_line = new_line.replace("%ENERGY", energy)

                error_mask_output.write(new_line)


if __name__ == "__main__":
    _model_twiss, _modifiers_path, _tunex, _tuney, _energy = _parse_args()
    process_masks(_model_twiss, _modifiers_path, _tunex, _tuney, _energy)
