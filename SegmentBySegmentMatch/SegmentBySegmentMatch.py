import os
import shutil
import optparse
from Utilities import iotools
import subprocess
from Python_Classes4MAD.metaclass import twiss
from Python_Classes4MAD import madxrunner


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
SBS_DATA_B1_PATH = "/user/slops/data/LHC_DATA/OP_DATA/Betabeat/10-7-2012/LHCB1/Results/0-4-46"
SBS_DATA_B2_PATH = "/user/slops/data/LHC_DATA/OP_DATA/Betabeat/10-7-2012/LHCB2/Results/0-5-39"


def parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--iteraction_point",
                    help="Which interaction point: 1, 2, 3...",
                    metavar="IP", default="8", dest="ip")
    (options, args) = parser.parse_args()
    return options, args


def main(options):
    ip = options.ip

    match_temporary_path = os.path.join(CURRENT_PATH, "match")
    beam1_temporary_path = os.path.join(match_temporary_path, "Beam1")
    beam2_temporary_path = os.path.join(match_temporary_path, "Beam2")

    iotools.create_dirs(os.path.join(beam1_temporary_path, "sbs"))
    iotools.create_dirs(os.path.join(beam2_temporary_path, "sbs"))

    print "Copying necessary files into temporary folder..."
    iotools.copy_item(os.path.join(CURRENT_PATH, "genconstraints.py"), match_temporary_path)
    iotools.copy_item(os.path.join(CURRENT_PATH, "genphases.py"), match_temporary_path)
    iotools.copy_item(os.path.join(CURRENT_PATH, "genvariables.py"), match_temporary_path)

    _copy_beam1_temp_files(ip, beam1_temporary_path)
    _apply_replace_to_beam1_files(beam1_temporary_path, os.path.join(beam1_temporary_path, "sbs"), ip)

    _copy_beam2_temp_files(ip, beam2_temporary_path)
    _apply_replace_to_beam2_files(beam2_temporary_path, os.path.join(beam2_temporary_path, "sbs"), ip)

    print "Getting matching range..."
    ((range_beam1_start_s, range_beam1_start_name),
    (range_beam1_end_s, range_beam1_end_name)) = _get_match_bpm_range(os.path.join(beam1_temporary_path, "sbs",
                                                                           "sbsphasext_IP" + ip + ".out"))
    range_beam1 = range_beam1_start_name + " / " + range_beam1_end_name

    ((range_beam2_start_s, range_beam2_start_name),
    (range_beam2_end_s, range_beam2_end_name)) = _get_match_bpm_range(os.path.join(beam2_temporary_path, "sbs",
                                                                           "sbsphasext_IP" + ip + ".out"))
    range_beam2 = range_beam2_start_name + " / " + range_beam2_end_name
    print "Matching range for Beam 1:", range_beam1_start_name, range_beam1_end_name
    print "Matching range for Beam 2:", range_beam2_start_name, range_beam2_end_name

    print "Running MADX..."
    madx_script_path = os.path.join(match_temporary_path, "matchIP" + ip + ".madx")
    _prepare_script_and_run_madx(ip, range_beam1, range_beam2, madx_script_path, match_temporary_path)

    print "Running GNUPlot..."
    _prepare_and_run_gnuplot(ip, match_temporary_path, range_beam1_start_s, range_beam1_end_s, range_beam2_start_s, range_beam2_end_s)

    print "Cleaning temporary folder..."
    _clean_up_temporary_dir(match_temporary_path)

    print "Done"


def _copy_beam1_temp_files(ip, beam1_temporary_path):
    _copy_files_with_extension(SBS_DATA_B1_PATH, beam1_temporary_path, ".out")
    _copy_files_with_extension(os.path.join(SBS_DATA_B1_PATH, "sbs"),
                               os.path.join(beam1_temporary_path, "sbs"),
                               ".madx")
    _copy_files_which_contains(os.path.join(SBS_DATA_B1_PATH, "sbs"),
                               os.path.join(beam1_temporary_path, "sbs"),
                               "IP" + str(ip))
    _copy_files_with_extension(os.path.join(SBS_DATA_B1_PATH, "sbs"),
                               os.path.join(beam1_temporary_path, "sbs"),
                               ".py")


def _apply_replace_to_beam1_files(beam1_temporary_path, beam1_temporary_sbs_path, ip):
    strings_to_replace_madx = [("//", "/"),
                               (SBS_DATA_B1_PATH, beam1_temporary_path),
                               ("stop;", "return;"),
                               ("install,", "!install,")]
    strings_to_replace_ip = [("! Write here some correction", "return;")]
    _replace_in_files_with_extension(beam1_temporary_sbs_path, strings_to_replace_madx, ".madx")
    _replace_in_files_with_extension(beam1_temporary_sbs_path, strings_to_replace_ip, "t_IP" + ip + ".madx")


def _copy_beam2_temp_files(ip, beam2_temporary_path):
    _copy_files_with_extension(os.path.join(SBS_DATA_B2_PATH), beam2_temporary_path, ".out")
    _copy_files_with_extension(os.path.join(SBS_DATA_B2_PATH, "sbs"),
                               os.path.join(beam2_temporary_path, "sbs"),
                               ".madx")
    _copy_files_which_contains(os.path.join(SBS_DATA_B2_PATH, "sbs"),
                               os.path.join(beam2_temporary_path, "sbs"), "IP" + str(ip))
    _copy_files_with_extension(os.path.join(SBS_DATA_B2_PATH, "sbs"),
                               os.path.join(beam2_temporary_path, "sbs"),
                               ".py")


def _apply_replace_to_beam2_files(beam2_temporary_path, beam2_temporary_sbs_path, ip):
    strings_to_replace_madx = [("//", "/"),
                               (SBS_DATA_B2_PATH, beam2_temporary_path),
                               ("stop;", "return;"),
                               ("install,", "!install,"),
                               ("label=b0", "label=b02"),
                               ("beta0=b0", "beta0=b02"),
                               ("label=b1", "label=b2"),
                               ("beta0=b1", "beta0=b2")]
    strings_to_replace_ip = [("! Write here some correction", "return;")]
    strings_to_replace_madx_slash = [(".*back propagation.*", ("return;"))]
    _replace_in_files_with_extension(beam2_temporary_sbs_path, strings_to_replace_madx, ".madx")
    _replace_in_files_with_extension(beam2_temporary_sbs_path, strings_to_replace_ip, "t_IP" + ip + ".madx")
    _replace_in_files_with_extension(beam2_temporary_sbs_path, strings_to_replace_madx_slash, "t_IP" + ip + ".madx", "/")


def _get_match_bpm_range(file_path):
    twiss_data = twiss(file_path)
    bpms_with_distances_list = zip(twiss_data.MODEL_S, twiss_data.NAME)
    bpms_with_distances_list.sort()
    return bpms_with_distances_list[0], bpms_with_distances_list[-1]


def _prepare_script_and_run_madx(ip, range_beam1, range_beam2, madx_script_path, match_temporary_path):
    iotools.copy_item("matchIP.madx", madx_script_path)
    _replace_in_file(madx_script_path, [("__IPNO__", str(ip)), ("__RANGEB1__", range_beam1), ("__RANGEB2__", range_beam2)])
    madx_binary_path = madxrunner.get_sys_dependent_path_to_mad_x()
    call_command = madx_binary_path + " < " + madx_script_path
    process = subprocess.Popen(call_command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               shell=True,
                               cwd=match_temporary_path)
    process.communicate()


def _prepare_and_run_gnuplot(ip, match_temporary_path, range_beam1_start_s, range_beam1_end_s, range_beam2_start_s, range_beam2_end_s):
    beam1_plot_path = os.path.join(match_temporary_path, "IP" + ip + "B1.gplot")
    beam2_plot_path = os.path.join(match_temporary_path, "IP" + ip + "B2.gplot")
    beam1_plot_replacements = [("__IPNO__", str(ip)), ("__BEAMNO__", "1"),
        ("__FILENAME__", "../IP" + str(ip) + "B1.eps"),
        ("__srangestart__", str(range_beam1_start_s)),
        ("__srangeend__", str(range_beam1_end_s))]
    beam2_plot_replacements = [("__IPNO__", str(ip)),
        ("__BEAMNO__", "2"),
        ("__FILENAME__", "../IP" + str(ip) + "B2.eps"),
        ("__srangestart__", str(range_beam2_start_s)),
        ("__srangeend__", str(range_beam2_end_s))]
    iotools.copy_item("templ.gplot", beam1_plot_path)
    iotools.copy_item("templ.gplot", beam2_plot_path)
    if str(ip) == "2":
        iotools.copy_item("templIP2B1.gplot", beam1_plot_path)
        beam1_qx, beam1_qy = _get_q_value(match_temporary_path, 1)
        beam1_plot_replacements.append(("__QX__", beam1_qx))
        beam1_plot_replacements.append(("__QY__", beam1_qy))
    elif str(ip) == "8":
        iotools.copy_item("templIP8B2.gplot", beam2_plot_path)
        beam2_qx, beam2_qy = _get_q_value(match_temporary_path, 2)
        beam2_plot_replacements.append(("__QX__", beam2_qx))
        beam2_plot_replacements.append(("__QY__", beam2_qy))
    _replace_in_file(beam1_plot_path, beam1_plot_replacements)
    _replace_in_file(beam2_plot_path, beam2_plot_replacements)
    proccess_beam1 = subprocess.Popen("gnuplot " + beam1_plot_path, shell=True, cwd=match_temporary_path)
    proccess_beam2 = subprocess.Popen("gnuplot " + beam2_plot_path, shell=True, cwd=match_temporary_path)
    proccess_beam1.communicate()
    proccess_beam2.communicate()


def _get_q_value(match_temporary_path, beam_num):
    with open(os.path.join(match_temporary_path, "Beam" + str(beam_num), "getphasex.out")) as file_content:
        for line in file_content:
            line_pieces = line.split(" ")
            if line_pieces[1] == "Q1":
                qx = line_pieces[3].strip()
            elif line_pieces[1] == "Q2":
                qy = line_pieces[3].strip()
    return qx, qy


def _clean_up_temporary_dir(match_temporary_path):
    os.unlink(os.path.join(match_temporary_path, "ats"))
    os.unlink(os.path.join(match_temporary_path, "db"))
    os.unlink(os.path.join(match_temporary_path, "db5"))
    os.unlink(os.path.join(match_temporary_path, "ds"))
    os.unlink(os.path.join(match_temporary_path, "lt"))
    iotools.delete_content_of_dir(match_temporary_path)


def _copy_files_with_extension(src, dest, ext):
    _copy_files_with_filter(src, dest, lambda file_name: file_name.endswith(ext))


def _copy_files_which_contains(src, dest, substring):
    _copy_files_with_filter(src, dest, lambda file_name: substring in file_name)


def _copy_files_with_filter(src, dest, filter_function):
    src_files = _get_filtered_file_list(src, filter_function)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        shutil.copy(full_file_name, dest)


def _replace_in_files_with_extension(src, replace_pairs, ext, sed_sep="#"):
    src_files = _get_filtered_file_list(src, lambda file_name: file_name.endswith(ext))
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        _replace_in_file(full_file_name, replace_pairs, sed_sep)


def _replace_in_file(full_file_name, replace_pairs, sed_sep="#"):  # TODO: Use a python function instead of sed
    sed_command = "sed -i "
    for pattern, replace in replace_pairs:
        full_command = sed_command + "'s" + sed_sep + pattern + sed_sep + replace + sed_sep + "g' " + full_file_name
        subprocess.call(full_command, shell=True)


def _get_filtered_file_list(src, filter_function):
    filtered_file_list = []
    original_file_list = os.listdir(src)
    for file_name in original_file_list:
        if os.path.isfile(os.path.join(src, file_name)) and filter_function(file_name):
            filtered_file_list.append(file_name)
    return filtered_file_list


if __name__ == "__main__":
    (options, args) = parse_args()
    main(options)
