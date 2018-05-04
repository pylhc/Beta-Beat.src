import os
import sys
import json
import datetime
import math
from optparse import OptionParser


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
AFS_ANACONDA_PATH = "/afs/cern.ch/work/o/omc/anaconda/bin/python"


def _parse_args():
    parser = OptionParser()
    (_, args) = parser.parse_args()
    if len(args) != 1:
        print >> sys.stderr, "Input only the path to the file to convert."
        sys.exit(-1)
    return args[0]


def main(input_file_path):
    esrf_data = loadmat(input_file_path)

    all_x = esrf_data["allx"]
    num_turns_x, num_bpms_x, num_sets_x = all_x.shape
    all_y = esrf_data["allz"]
    num_turns_y, num_bpms_y, num_sets_y = all_y.shape

    if num_turns_x != num_turns_y:
        raise ValueError("Number of turns in x and y do not match")

    if num_sets_x != num_sets_y:
        raise ValueError("Number of measurements in x and y do not match")

    bpm_data = json.load(open(os.path.join(CURRENT_PATH, "bpm_names.json"), "r"))

    previous_file_tbt_data = None
    this_file_tbt_data = None
    for set_number in range(num_sets_x):
        set_string = str(set_number + 1) if len(str(set_number + 1)) > 1 else "0" + str(set_number + 1)
        output_file_name = os.path.basename(input_file_path).replace(".mat", "." + set_string + '.sdds')
        output_file_path = os.path.join(os.path.dirname(input_file_path), output_file_name)
        previous_file_tbt_data = this_file_tbt_data
        this_file_tbt_data = {}
        this_file_tbt_data["0"] = {}
        this_file_tbt_data["1"] = {}
        with open(output_file_path, "w") as output_file:
            output_file.write("#SDDSASCIIFORMAT v1 \n")
            output_file.write("#Beam: ESRFonly1 \n")
            output_file.write("#Created: " + str(datetime.datetime.now()) + " by ESRF mat to sdds converter\n")
            output_file.write("#numberofturns : " + str(num_turns_x) + "\n",)
            output_file.write("#numberofmonitors: " + str(num_bpms_x + num_bpms_y) + " \n")
            output_file.write("#NTURNS: " + str(num_turns_x) + "\n")

            for bpm_number in range(num_bpms_x):
                _write_sdds_data_line(bpm_number, set_number, all_x, "0", bpm_data, output_file,
                                      this_file_tbt_data, previous_file_tbt_data)

            for bpm_number in range(num_bpms_y):
                _write_sdds_data_line(bpm_number, set_number, all_y, "1", bpm_data, output_file,
                                      this_file_tbt_data, previous_file_tbt_data)


def _write_sdds_data_line(bpm_number, set_number, all_data, plane_string, bpm_data, output_file, this_file_tbt_data, previous_file_tbt_data):
    bpm_name = bpm_data["names"][bpm_number].replace("PICK", "BPM")
    bpm_s = bpm_data["s"][bpm_number]
    tbt_data = map(float, map(lambda x: "{0:f}".format(x), all_data[:, bpm_number, set_number]))
    if _check_tbt_data(bpm_name, plane_string, tbt_data, previous_file_tbt_data):
        output_file.write(plane_string + " " + bpm_name + " " + str(bpm_s) + " ")
        this_file_tbt_data[plane_string][bpm_name] = tbt_data
        output_file.write(" ".join(map(str, tbt_data)))
        output_file.write("\n")


def _check_tbt_data(bpm_name, plane_string, tbt_data, previous_file_tbt_data):
    if any(map(math.isnan, tbt_data)):
        print "Ignoring BPM", bpm_name, "found NaNs for plane", plane_string
        return False
    elif previous_file_tbt_data is None:
        return True
    else:
        if bpm_name in previous_file_tbt_data[plane_string]:
            if previous_file_tbt_data[plane_string][bpm_name] == tbt_data:
                print "Ignoring BPM", bpm_name, "found repeated values from previous file for plane", plane_string
                return False
    return True


def _launch_with_afs_anaconda(input_file_path):
    command = AFS_ANACONDA_PATH + " " + os.path.abspath(__file__) + " " + input_file_path
    print "Cannot import scipy, trying using afs version, command: " + command
    os.system(command)


if __name__ == "__main__":
    _input_file_path = _parse_args()
    try:
        from scipy.io import loadmat
        main(_input_file_path)
    except ImportError:
        _launch_with_afs_anaconda(_input_file_path)
