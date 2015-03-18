import __init__  # @UnusedImport used for appending paths
import os
import shutil
import optparse
from Utilities import iotools
import subprocess
from Python_Classes4MAD.metaclass import twiss
from Python_Classes4MAD import madxrunner
import json
import numpy as np
from SegmentBySegment import SegmentBySegment
import Utilities


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
ALL_LISTS_BEAM1_PATH = os.path.join(CURRENT_PATH, '..', 'MODEL', 'LHCB', 'fullresponse', 'LHCB1', 'AllLists.json')
ALL_LISTS_BEAM2_PATH = os.path.join(CURRENT_PATH, '..', 'MODEL', 'LHCB', 'fullresponse', 'LHCB2', 'AllLists.json')

ERROR_CONSTRAINT_FACTOR = 1000
MAX_WEIGHT = 4


def parse_args():
    parser = optparse.OptionParser()
    parser.add_option("--run",
                    help="LHC run 1 or 2",
                    metavar="RUN", default="1", dest="lhc_run")
    parser.add_option("--ip",
                    help="Which interaction point: 1, 2, 3...",
                    metavar="IP", default="1", dest="ip")
    parser.add_option("--beam1",
                    help="Path to the measurement files for beam 1",
                    metavar="BEAM1", dest="b1")
    parser.add_option("--beam2",
                    help="Path to the measurement files for beam 2",
                    metavar="BEAM2", dest="b2")
    parser.add_option("-t", "--temp",
                    help="Path to the a temporary folder",
                    metavar="TEMP", default="", dest="temp")
    parser.add_option("--exclude",
                    help="Variables to exclude",
                    metavar="EXCLUDE", default="", dest="exclude")
    parser.add_option("-r", "--useerrors",
                    help="Use errors in constraint generation",
                    metavar="USE_ERRORS", action="store_true", default=False, dest="use_errors")
    (options, args) = parser.parse_args()
    return options, args


def main(options, args):
    if len(args) > 0:
        command = args[0]
    else:
        command = "match"
    ip = options.ip
    temporary_path = options.temp
    match_temporary_path = os.path.join(temporary_path, "match")
    if command == "variables":
        exclude_vars_string = options.exclude
        generate_variables(ip, match_temporary_path, exclude_vars_string)
    elif command == "constraints":
        sbs_data_b1_path = options.b1
        sbs_data_b2_path = options.b2
        exclude_constr_string = options.exclude
        generate_constraints(ip, options.use_errors, sbs_data_b1_path, sbs_data_b2_path, match_temporary_path, exclude_constr_string)
    elif command == "clean":
        clean_up_temporary_dir(match_temporary_path)
    else:
        sbs_data_b1_path = options.b1
        sbs_data_b2_path = options.b2
        temporary_path = options.temp
        match(options.lhc_run, ip, sbs_data_b1_path, sbs_data_b2_path, match_temporary_path, options.use_errors)


def match(lhc_run, ip, sbs_data_b1_path, sbs_data_b2_path, match_temporary_path, use_errors):

    print "+++ Starting Segment by Segment Match +++"

    beam1_temporary_path = os.path.join(match_temporary_path, "Beam1")
    beam2_temporary_path = os.path.join(match_temporary_path, "Beam2")

    iotools.create_dirs(os.path.join(beam1_temporary_path, "sbs"))
    iotools.create_dirs(os.path.join(beam2_temporary_path, "sbs"))

    _check_and_run_genvariables(ip, match_temporary_path)
    _check_and_run_genconstraints(ip, sbs_data_b1_path, sbs_data_b2_path, match_temporary_path, use_errors)
    run_genphases(ip, match_temporary_path, sbs_data_b1_path, sbs_data_b2_path)

    print "Copying files into temporary folder..."
    iotools.copy_item(os.path.join(CURRENT_PATH, "dumpB1.gplot"), match_temporary_path)
    iotools.copy_item(os.path.join(CURRENT_PATH, "dumpB2.gplot"), match_temporary_path)

    _copy_beam1_temp_files(ip, sbs_data_b1_path, beam1_temporary_path)
    _copy_beam2_temp_files(ip, sbs_data_b2_path, beam2_temporary_path)

    print "Getting matching range..."
    ((range_beam1_start_s, range_beam1_start_name),
    (range_beam1_end_s, range_beam1_end_name)) = _get_match_bpm_range(os.path.join(beam1_temporary_path, "sbs",
                                                                           "sbsphasext_IP" + ip + ".out"))
    ((range_beam2_start_s, range_beam2_start_name),
    (range_beam2_end_s, range_beam2_end_name)) = _get_match_bpm_range(os.path.join(beam2_temporary_path, "sbs",
                                                                           "sbsphasext_IP" + ip + ".out"))
    print "Matching range for Beam 1:", range_beam1_start_name, range_beam1_end_name
    print "Matching range for Beam 2:", range_beam2_start_name, range_beam2_end_name

    print "Running MADX..."
    label = "IP" + str(ip)
    _prepare_script_and_run_madx(lhc_run, label, beam1_temporary_path, beam2_temporary_path, match_temporary_path,
                                 range_beam1_start_name, range_beam1_end_name,
                                 range_beam2_start_name, range_beam2_end_name)

    print "Writting sbs files from MADX results..."
    _write_sbs_data(ip, beam1_temporary_path, beam2_temporary_path, range_beam1_start_name, range_beam2_start_name)

    print "Building changeparameters_match.madx..."
    _build_changeparameters_file(match_temporary_path)

    print "Running GNUPlot..."
    _prepare_and_run_gnuplot(ip, match_temporary_path,
                             range_beam1_start_s, range_beam1_end_s, range_beam2_start_s, range_beam2_end_s)

    print "+++  Ended Segment by Segment Match  +++"
    return 0


def _check_and_run_genvariables(ip, match_temporary_path):
    for file_name in ["applycorrection.seqx", "dvariables.seqx", "genchangpars.seqx",
                      "svariables.seqx", "variablesb1.seqx", "variablesb2.seqx", "variablesc.seqx"]:
        full_file_path = os.path.join(match_temporary_path, file_name)
        if not os.path.exists(full_file_path):  # TODO: Here the variables should be recreated if they are for a different IP
            print "File " + file_name + " not found, generating new variables files..."
            generate_variables(ip, match_temporary_path)
            break


def _check_and_run_genconstraints(ip, sbs_data_b1_path, sbs_data_b2_path, match_temporary_path, use_errors):
    for file_name in ["constraintsb1.seqx", "constraintsb2.seqx"]:
        full_file_path = os.path.join(match_temporary_path, file_name)
        if not os.path.exists(full_file_path):  # TODO: Here the constraints should be recreated if they are for a different IP
            print "File " + file_name + " not found, generating new constraints files..."
            generate_constraints(ip, use_errors, sbs_data_b1_path, sbs_data_b2_path, match_temporary_path)
            break


def generate_variables(ip, variables_path=os.path.join(CURRENT_PATH, "match"), exclude_string=""):
    if not os.path.exists(variables_path):
        iotools.create_dirs(variables_path)

    exclude_list = _parse_exclude_string(exclude_string)

    variables_beam1 = json.load(file(ALL_LISTS_BEAM1_PATH, 'r'))['getListsByIR'][1]
    variables_common, variables_beam2 = json.load(file(ALL_LISTS_BEAM2_PATH, 'r'))['getListsByIR']

    ip_string = str(ip)

    apply_correction_file = open(os.path.join(variables_path, "applycorrection.seqx"), 'w')
    variables_common_file = open(os.path.join(variables_path, "variablesc.seqx"), 'w')
    variables_beam1_file = open(os.path.join(variables_path, "variablesb1.seqx"), 'w')
    variables_beam2_file = open(os.path.join(variables_path, "variablesb2.seqx"), 'w')
    variables_s_file = open(os.path.join(variables_path, "svariables.seqx"), 'w')
    variables_d_file = open(os.path.join(variables_path, "dvariables.seqx"), 'w')

    param_change_generator_file = open(os.path.join(variables_path, "genchangpars.seqx"), 'w')
    param_change_generator_file.write('select,flag=save, clear;')

    variables = variables_beam1[ip_string]
    param_change_generator_file.write('!B1\n')
    _vars_to_files(apply_correction_file, variables_beam1_file, variables_s_file,
                   variables_d_file, param_change_generator_file, variables, exclude_list)

    variables = variables_beam2[ip_string]
    param_change_generator_file.write('\n!B2\n')
    _vars_to_files(apply_correction_file, variables_beam2_file, variables_s_file,
                   variables_d_file, param_change_generator_file, variables, exclude_list)

    variables = variables_common[ip_string]
    param_change_generator_file.write('\n!B1 and B2\n')
    _vars_to_files(apply_correction_file, variables_common_file, variables_s_file,
                   variables_d_file, param_change_generator_file, variables, exclude_list)

    variables_common_file.close()
    variables_beam1_file.close()
    variables_beam2_file.close()
    variables_s_file.close()
    variables_d_file.close()

    param_change_generator_file.write('\n save, file=\"changeparameters.madx\";\n')
    param_change_generator_file.close()


def _vars_to_files(apply_correction_file, variables_file, variables_s_file, variables_d_file, param_change_generator_file, variables, excluded):
    for variable in variables:
        if variable not in excluded:
            variables_file.write('   vary, name=d' + variable + ', step:=1e-4;\n')
        variables_s_file.write(' ' + variable + '_0 = ' + variable + ';\n')
        variables_d_file.write(' ' + variable + ' := ' + variable + '_0 + d' + variable + ';\n')
        param_change_generator_file.write('select,flag=save,pattern=\"d' + variable + '\";\n')
        apply_correction_file.write(variable + ' = ' + variable + '_0 + d' + variable + ';\n')


def generate_constraints(ip, use_errors, sbs_data_b1_path, sbs_data_b2_path, constraints_path=os.path.join(CURRENT_PATH, "match"), exclude_string=""):
    full_data_beam1 = twiss(os.path.join(sbs_data_b1_path, 'getphasex.out'))
    x_tune_beam1 = full_data_beam1.Q1
    y_tune_beam1 = full_data_beam1.Q2

    full_data_beam2 = twiss(os.path.join(sbs_data_b2_path, 'getphasex.out'))
    x_tune_beam2 = full_data_beam2.Q1
    y_tune_beam2 = full_data_beam2.Q2

    sbs_x_data_beam1 = twiss(os.path.join(sbs_data_b1_path, 'sbs', 'sbsphasext_IP' + ip + '.out'))
    sbs_y_data_beam1 = twiss(os.path.join(sbs_data_b1_path, 'sbs', 'sbsphaseyt_IP' + ip + '.out'))
    sbs_x_data_beam2 = twiss(os.path.join(sbs_data_b2_path, 'sbs', 'sbsphasext_IP' + ip + '.out'))
    sbs_y_data_beam2 = twiss(os.path.join(sbs_data_b2_path, 'sbs', 'sbsphaseyt_IP' + ip + '.out'))

    constr_file_beam1 = open(os.path.join(constraints_path, 'constraintsb1.seqx'), 'w')
    constr_file_beam2 = open(os.path.join(constraints_path, 'constraintsb2.seqx'), 'w')

    if exclude_string.strip() == "":
        exclude_list_x = ""
        exclude_list_y = ""
    else:
        exclude_both_planes = exclude_string.split(";")
        exclude_list_x = _parse_exclude_string(exclude_both_planes[0])
        exclude_list_y = _parse_exclude_string(exclude_both_planes[1])

    _write_constraints_file(sbs_x_data_beam1, constr_file_beam1, ip, 1, "x", x_tune_beam1, exclude_list_x, use_errors)
    _write_constraints_file(sbs_y_data_beam1, constr_file_beam1, ip, 1, "y", y_tune_beam1, exclude_list_y, use_errors)
    _write_constraints_file(sbs_x_data_beam2, constr_file_beam2, ip, 2, "x", x_tune_beam2, exclude_list_x, use_errors)
    _write_constraints_file(sbs_y_data_beam2, constr_file_beam2, ip, 2, "y", y_tune_beam2, exclude_list_y, use_errors)


def _write_constraints_file(sbs_data, constr_file, ip, beam, plane, tune, exclude_list, use_errors):
    if plane == "x":
        constr_file.write('\n!!!! BEAM ' + str(beam) + ' H !!!!!\n\n')
    else:
        constr_file.write('\n!!!! BEAM ' + str(beam) + ' V !!!!!\n\n')

    for index in range(0, len(sbs_data.NAME)):
        name = sbs_data.NAME[index]
        if name not in exclude_list:
            phase = sbs_data.PROPPHASEX[index] if plane == "x" else sbs_data.PROPPHASEY[index]
            error = sbs_data.ERRPROPPHASEX[index] if plane == "x" else sbs_data.ERRPROPPHASEY[index]
            s = sbs_data.S[index]

            if abs(phase) > 0.25 or error == 0.:
                weight = 1e-6
            elif not use_errors:
                weight = 1
            else:
                weight = 1 / ((ERROR_CONSTRAINT_FACTOR * error) ** 2)
                if weight > MAX_WEIGHT:
                    weight = MAX_WEIGHT

            constr_file.write('   constraint, weight = ' + str(weight) + ' , ')
            constr_file.write('expr =  dmu' + plane + name + ' = ' + str(phase) + '; ')

            constr_file.write('!   S = ' + str(s))
            constr_file.write(';\n')


def _parse_exclude_string(exclude_string):
    if not exclude_string == "":
        exclude_list = [var_name.strip() for var_name in exclude_string.strip('"').split(",")]
    else:
        exclude_list = []
    return exclude_list


# TODO: Refactor this function
def run_genphases(ip, match_temporary_path, sbs_data_b1_path, sbs_data_b2_path):

    sbs_x_data_beam1 = twiss(os.path.join(sbs_data_b1_path, 'sbs', 'sbsphasext_IP' + ip + '.out'))
    sbs_y_data_beam1 = twiss(os.path.join(sbs_data_b1_path, 'sbs', 'sbsphaseyt_IP' + ip + '.out'))
    sbs_x_data_beam2 = twiss(os.path.join(sbs_data_b2_path, 'sbs', 'sbsphasext_IP' + ip + '.out'))
    sbs_y_data_beam2 = twiss(os.path.join(sbs_data_b2_path, 'sbs', 'sbsphaseyt_IP' + ip + '.out'))

    phases_file = open(os.path.join(match_temporary_path, "phases.seqx"), 'w')
    phases0_beam1_file = open(os.path.join(match_temporary_path, "phases0b1.seqx"), 'w')
    phases0_beam2_file = open(os.path.join(match_temporary_path, "phases0b2.seqx"), 'w')

    beam1_s_list = sbs_x_data_beam1.MODEL_S
    sorted_index_beam1 = np.argsort(beam1_s_list)

    beam2_s_list = sbs_x_data_beam2.MODEL_S
    sorted_index_beam2 = np.argsort(beam2_s_list)

    phases_file.write('\n!!!! BEAM 1 H !!!!!\n\n')
    phases0_beam1_file.write('\n!!!! BEAM 1 H !!!!!\n\n')

    for name in sbs_x_data_beam1.NAME:
        phases0_beam1_file.write('mux0' + name + ' = ')
        phases0_beam1_file.write('table(twiss, ' + name + ', mux) - ')
        phases0_beam1_file.write('table(twiss, ' + sbs_x_data_beam1.NAME[sorted_index_beam1[0]] + ', mux);\n')

        phases_file.write('mux' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', mux) - ')
        phases_file.write('table(twiss, ' + sbs_x_data_beam1.NAME[sorted_index_beam1[0]] + ', mux);\n')

    phases_file.write('\n')

    for name in sbs_x_data_beam1.NAME:
        phases_file.write('dmux' + name + ' := ')
        phases_file.write('mux' + name + ' - ' 'mux0' + name + ';\n')

    phases_file.write('\n!!!! BEAM 2 H !!!!!\n\n')

    phases0_beam2_file.write('\n!!!! BEAM 2 H !!!!!\n\n')

    for idx in sorted_index_beam2:
        name = sbs_x_data_beam2.NAME[idx]

        phases0_beam2_file.write('mux0' + name + ' = ')
        phases0_beam2_file.write('table(twiss, ' + name + ', mux) - ')
        phases0_beam2_file.write('table(twiss, ' + sbs_x_data_beam2.NAME[sorted_index_beam2[0]] + ', mux);\n')

        phases_file.write('mux' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', mux) - ')
        phases_file.write('table(twiss, ' + sbs_x_data_beam2.NAME[sorted_index_beam2[0]] + ', mux);\n')

    phases_file.write('\n')

    for name in sbs_x_data_beam2.NAME:
        phases_file.write('dmux' + name + ' := ')
        phases_file.write('mux' + name + ' - ' 'mux0' + name + ';\n')

    beam1_s_list = sbs_y_data_beam1.MODEL_S
    sorted_index_beam1 = np.argsort(beam1_s_list)

    beam2_s_list = sbs_y_data_beam2.MODEL_S
    sorted_index_beam2 = np.argsort(beam2_s_list)

    phases_file.write('\n!!!! BEAM 1 V !!!!!\n\n')

    phases0_beam1_file.write('\n!!!! BEAM 1 V !!!!!\n\n')

    for name in sbs_y_data_beam1.NAME:

        phases0_beam1_file.write('muy0' + name + ' = ')
        phases0_beam1_file.write('table(twiss, ' + name + ', muy) - ')
        phases0_beam1_file.write('table(twiss, ' + sbs_y_data_beam1.NAME[sorted_index_beam1[0]] + ', muy);\n')

        phases_file.write('muy' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', muy) - ')
        phases_file.write('table(twiss, ' + sbs_y_data_beam1.NAME[sorted_index_beam1[0]] + ', muy);\n')

    phases_file.write('\n')

    for name in sbs_y_data_beam1.NAME:
        phases_file.write('dmuy' + name + ' := ')
        phases_file.write('muy' + name + ' - ' 'muy0' + name + ';\n')

    phases_file.write('\n!!!! BEAM 2 V !!!!!\n\n')

    phases0_beam2_file.write('\n!!!! BEAM 2 V !!!!!\n\n')

    for name in sbs_y_data_beam2.NAME:
        phases0_beam2_file.write('muy0' + name + ' = ')
        phases0_beam2_file.write('table(twiss, ' + name + ', muy) - ')
        phases0_beam2_file.write('table(twiss, ' + sbs_y_data_beam2.NAME[sorted_index_beam2[0]] + ', muy);\n')

        phases_file.write('muy' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', muy) - ')
        phases_file.write('table(twiss, ' + sbs_y_data_beam2.NAME[sorted_index_beam2[0]] + ', muy);\n')

    phases_file.write('\n')

    for name in sbs_y_data_beam2.NAME:
        phases_file.write('dmuy' + name + ' := ')
        phases_file.write('muy' + name + ' - ' 'muy0' + name + ';\n')

    phases0_beam1_file.close()
    phases0_beam2_file.close()
    phases_file.close()


def _copy_beam1_temp_files(ip, sbs_data_b1_path, beam1_temporary_path):
    _copy_files_with_extension(sbs_data_b1_path, beam1_temporary_path, ".out")
    _copy_files_with_extension(os.path.join(sbs_data_b1_path, "sbs"),
                               os.path.join(beam1_temporary_path, "sbs"),
                               ".madx")
    _copy_files_which_contains(os.path.join(sbs_data_b1_path, "sbs"),
                               os.path.join(beam1_temporary_path, "sbs"),
                               "IP" + str(ip))
    _copy_files_with_extension(os.path.join(sbs_data_b1_path, "sbs"),
                               os.path.join(beam1_temporary_path, "sbs"),
                               ".py")


def _copy_beam2_temp_files(ip, sbs_data_b2_path, beam2_temporary_path):
    _copy_files_with_extension(os.path.join(sbs_data_b2_path), beam2_temporary_path, ".out")
    _copy_files_with_extension(os.path.join(sbs_data_b2_path, "sbs"),
                               os.path.join(beam2_temporary_path, "sbs"),
                               ".madx")
    _copy_files_which_contains(os.path.join(sbs_data_b2_path, "sbs"),
                               os.path.join(beam2_temporary_path, "sbs"), "IP" + str(ip))
    _copy_files_with_extension(os.path.join(sbs_data_b2_path, "sbs"),
                               os.path.join(beam2_temporary_path, "sbs"),
                               ".py")


def _get_match_bpm_range(file_path):
    twiss_data = twiss(file_path)
    bpms_with_distances_list = zip(twiss_data.S, twiss_data.NAME)
    bpms_with_distances_list.sort()
    return bpms_with_distances_list[0], bpms_with_distances_list[-1]


def _prepare_script_and_run_madx(lhc_run, label, beam1_temporary_path, beam2_temporary_path, match_temporary_path,
                                 b1_range_start, b1_range_end, b2_range_start, b2_range_end):
    bb_path = os.path.abspath(os.path.join(CURRENT_PATH, ".."))
    sbs_path = os.path.join(bb_path, "SegmentBySegment")
    dict_for_replacing = dict(
        LHC_RUN=lhc_run,
        PATHB1=os.path.join(beam1_temporary_path, "sbs"),
        PATHB2=os.path.join(beam2_temporary_path, "sbs"),
        MATCH=match_temporary_path,
        LABEL=label,
        BBPATH=bb_path,
        SBSPATH=sbs_path,
        STARTFROMB1=b1_range_start,
        ENDATB1=b1_range_end,
        STARTFROMB2=b2_range_start,
        ENDATB2=b2_range_end
        )

    mask_file = os.path.join(CURRENT_PATH, "job.match.madx")
    madx_script_path = os.path.join(match_temporary_path, "job.match" + label + ".madx")

    Utilities.iotools.replace_keywords_in_textfile(mask_file, dict_for_replacing, madx_script_path)

    madxrunner.runForInputFile(madx_script_path, stdout=open(os.path.join(match_temporary_path, "match_madx_out.log"), "w"))


def _write_sbs_data(ip, beam1_temporary_path, beam2_temporary_path, range_beam1_start_name, range_beam2_start_name):
    save_path_b1 = os.path.join(beam1_temporary_path, "sbs")
    save_path_b2 = os.path.join(beam2_temporary_path, "sbs")
    input_data_b1 = SegmentBySegment._InputData(beam1_temporary_path)
    input_data_b2 = SegmentBySegment._InputData(beam2_temporary_path)
    prop_models_b1 = SegmentBySegment._PropagatedModels(save_path_b1, "IP" + str(ip))
    prop_models_b2 = SegmentBySegment._PropagatedModels(save_path_b2, "IP" + str(ip))

    SegmentBySegment.getAndWriteData("IP" + ip, input_data_b1, None, prop_models_b1, save_path_b1, False, False, False, "LHCB1", None)
    SegmentBySegment.getAndWriteData("IP" + ip, input_data_b2, None, prop_models_b2, save_path_b2, False, False, False, "LHCB2", None)


def _prepare_and_run_gnuplot(ip, match_temporary_path, range_beam1_start_s, range_beam1_end_s, range_beam2_start_s, range_beam2_end_s):
    beam1_plot_path = os.path.join(match_temporary_path, "IP" + ip + "B1.gplot")
    beam2_plot_path = os.path.join(match_temporary_path, "IP" + ip + "B2.gplot")
    beam1_plot_replacements = dict(
                                   IPNO=str(ip),
                                   BEAMNO="1",
                                   FILENAME="IP" + str(ip) + "B1.eps",
                                   SRANGESTART=str(range_beam1_start_s),
                                   SRANGEEND=str(range_beam1_end_s)
                                   )

    beam2_plot_replacements = dict(
                                   IPNO=str(ip),
                                   BEAMNO="2",
                                   FILENAME="IP" + str(ip) + "B2.eps",
                                   SRANGESTART=str(range_beam2_start_s),
                                   SRANGEEND=str(range_beam2_end_s)
                                   )
    beam1_plot_template = os.path.join(CURRENT_PATH, "templ.gplot")
    beam2_plot_template = os.path.join(CURRENT_PATH, "templ.gplot")

    Utilities.iotools.replace_keywords_in_textfile(beam1_plot_template, beam1_plot_replacements, beam1_plot_path)
    Utilities.iotools.replace_keywords_in_textfile(beam2_plot_template, beam2_plot_replacements, beam2_plot_path)

    proccess_beam1 = subprocess.Popen("gnuplot " + beam1_plot_path, shell=True, cwd=match_temporary_path)
    proccess_beam2 = subprocess.Popen("gnuplot " + beam2_plot_path, shell=True, cwd=match_temporary_path)

    proccess_beam1.communicate()
    proccess_beam2.communicate()


def clean_up_temporary_dir(match_temporary_path):
    os.unlink(os.path.join(match_temporary_path, "ats"))
    os.unlink(os.path.join(match_temporary_path, "db"))
    os.unlink(os.path.join(match_temporary_path, "db5"))
    os.unlink(os.path.join(match_temporary_path, "ds"))
    os.unlink(os.path.join(match_temporary_path, "lt"))
    iotools.delete_content_of_dir(match_temporary_path)


def _build_changeparameters_file(match_temporary_path):
    original_changeparameters_file = open(os.path.join(match_temporary_path, "changeparameters.madx"), "r")
    changeparameters_match_file = open(os.path.join(match_temporary_path, "changeparameters_match.madx"), "w")
    for original_line in original_changeparameters_file.readlines():
        parts = original_line.split("=")
        variable_name = parts[0].replace("d", "", 1).strip()
        variable_value = -float(parts[1].replace(";", "").strip())
        if variable_value < 0.0:
            print >> changeparameters_match_file, variable_name, " = ", variable_name, " - ", abs(variable_value), ";"
        else:
            print >> changeparameters_match_file, variable_name, " = ", variable_name, " + ", abs(variable_value), ";"
    print >> changeparameters_match_file, "return;"


def _copy_files_with_extension(src, dest, ext):
    _copy_files_with_filter(src, dest, lambda file_name: file_name.endswith(ext))


def _copy_files_which_contains(src, dest, substring):
    _copy_files_with_filter(src, dest, lambda file_name: substring in file_name)


def _copy_files_with_filter(src, dest, filter_function):
    src_files = _get_filtered_file_list(src, filter_function)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        shutil.copy(full_file_name, dest)


def _get_filtered_file_list(src, filter_function):
    filtered_file_list = []
    original_file_list = os.listdir(src)
    for file_name in original_file_list:
        if os.path.isfile(os.path.join(src, file_name)) and filter_function(file_name):
            filtered_file_list.append(file_name)
    return filtered_file_list


if __name__ == "__main__":
    (options, args) = parse_args()
    main(options, args)
