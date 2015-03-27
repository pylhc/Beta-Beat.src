import __init__  # @UnusedImport
from optparse import OptionParser
from Python_Classes4MAD import metaclass
import os
import sys
import json
import re


def _parse_args():
    parser = OptionParser()
    parser.add_option("-m", "--model",
                    help="Model where to read dispersion data",
                    metavar="MODEL", dest="model")
    parser.add_option("-f", "--files",
                    help="List of comma-separated sdds files from analysis",
                    metavar="FILES", dest="files")
    parser.add_option("-o", "--output",
                    help="File where to output the results",
                    metavar="OUT", dest="output_file")
    options, _ = parser.parse_args()

    return options.model, options.files, options.output_file


def compute_dp_over_p(model, files, output_file):
    model_twiss, lin_x_twiss_list = _parse_twiss_files(model, files)
    file_to_dpp_dict = {}
    for lin_x_twiss in lin_x_twiss_list:
        dp_over_p = _single_file_dp_over_p(lin_x_twiss, model_twiss)
        file_to_dpp_dict[os.path.basename(lin_x_twiss.filename)] = dp_over_p
    _json_output(file_to_dpp_dict, output_file)


def _single_file_dp_over_p(lin_x_twiss, model_twiss):
    numer = 0
    denom = 0
    for bpm_name in lin_x_twiss.NAME:
        if bpm_name in model_twiss.NAME:
            if _is_arc_bpm(bpm_name):
                dispersion = model_twiss.DX[model_twiss.indx[bpm_name]] * 1e3  # We need it in mm
                closed_orbit = lin_x_twiss.CO[lin_x_twiss.indx[bpm_name]]
                numer += dispersion * closed_orbit
                denom += dispersion ** 2
    dp_over_p = numer / denom
    return dp_over_p


def _is_arc_bpm(bpm_name):
    match = re.search("BPM.*\.([0-9]+)[RL].\..*", bpm_name, re.IGNORECASE)
    if match:
        return int(match.group(1)) > 14  # The arc bpms are from BPM.14... and up
    else:
        return False


def _parse_twiss_files(model, files):
    model_twiss = metaclass.twiss(model)
    lin_x_twiss_list = []
    file_list = [file_name.strip() for file_name in files.strip("\"").split(",")]
    for file_name in file_list:
        if not file_name.endswith("_linx"):
            file_name += "_linx"
        if not os.path.isfile(file_name):
            print >> sys.stderr, "Cannot find file_name: " + file_name + ". Aborting."
            sys.exit()
        else:
            lin_x_twiss_list.append(metaclass.twiss(file_name))
    return model_twiss, lin_x_twiss_list


def _json_output(file_to_dpp_dict, output_file):
    json_string = json.dumps(file_to_dpp_dict)
    print json_string
    with open(output_file, "w") as output_file_data:
        output_file_data.write(json_string)


if __name__ == "__main__":
    _model, _files, _output_file = _parse_args()
    compute_dp_over_p(_model, _files, _output_file)
