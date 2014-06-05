import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import os
import optparse
import metaclass  # @UnresolvedImport
import math


def parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-l", "--label",
                    help="Label: IP1 IP2 IP3 IP4...",
                    metavar="LABEL", default="", dest="label")
    parser.add_option("-p", "--path",
                    help="Path to model files and output",
                    metavar="PATH", default="./", dest="path")
    parser.add_option("-f", "--fast",
                    help="1 for very fast calculation only writing initvals file",
                    metavar="FAST", default="0", dest="fast")
    parser.add_option("-e", "--exp",
                    help="path to experimental files, only used if fast!=1",
                    metavar="EXP", default="./", dest="exp")
    parser.add_option("-s", "--start",
                    help="Start BPM",
                    metavar="START", default="./", dest="start")
    parser.add_option("-m", "--method",
                    help="Method",
                    metavar="ME", default="_free", dest="ME")
    parser.add_option("-w", "--w",  # Path to Chromaticity functions
                        help="Path to  chromaticity functions, by default this is skiped",
                        metavar="wpath", default="0", dest="wpath")
    options = parser.parse_args()[0]
    main(options.path, options.exp, options.label, options.fast, options.start, options.ME, options.wpath)


def main(output_path, experimental_path, label, fast, start_bpm, method, w_path):
    model_twiss = _try_to_load_twiss(output_path, "twiss_" + label + ".dat")
    if not model_twiss:
        print >> sys.stderr, "Can't continue without model..."
        return

    method = method
    if method == "driven":
        method = ""

    # This calculates the coupling matrix with the f values
    model_twiss.Cmatrix()
    elements_names = model_twiss.NAME

    write_initial_values_file(output_path, model_twiss, elements_names)

    if fast != "1":
        model_play_twiss = _try_to_load_twiss(output_path, "twiss_" + label + "_cor.dat")
        if not model_play_twiss:
            print >> sys.stderr, "Can't continue without model play..."
            return
        model_play_twiss.Cmatrix()

        write_coupling_files(output_path, experimental_path, label, method, model_twiss, model_play_twiss)
        write_phase_files(output_path, experimental_path, label, method, start_bpm, model_twiss, model_play_twiss)
        write_w_function_files(output_path, w_path, label, model_twiss, model_play_twiss)
        write_dispersion_files(output_path, experimental_path, label, model_twiss, model_play_twiss)


def write_initial_values_file(path, model_twiss, elements_names):
    initial_values_file_path = os.path.join(path, "initvals.dat")
    initial_values_file = open(initial_values_file_path, "w")
    f1001 = model_twiss.f1001[model_twiss.indx[elements_names[0]]]
    f1010 = model_twiss.f1010[model_twiss.indx[elements_names[0]]]
    print >> initial_values_file, "f1001r=", f1001.real, ";"
    print >> initial_values_file, "f1001i=", f1001.imag, ";"
    print >> initial_values_file, "f1010r=", f1010.real, ";"
    print >> initial_values_file, "f1010i=", f1010.imag, ";"
    initial_values_file.close()


def write_coupling_files(output_path, experimental_path, label, method, model_twiss, model_play_twiss):
    sbs_coupling_file_path = os.path.join(output_path, "sbscouple_" + label + ".out")
    sbs_coupling_file = open(sbs_coupling_file_path, "w")

    exp_coupling_file = _try_to_load_twiss(experimental_path, "getcouple" + method + ".out")
    if not exp_coupling_file:
        print >> sys.stderr, "Won't write coupling file"
        return

    print >> sbs_coupling_file, _get_headers_for("coupling")
    print >> sbs_coupling_file, _get_data_type_header_for("coupling")

    elements_names = model_twiss.NAME
    for element_name in elements_names:
        if element_name in exp_coupling_file.indx:

            element_index = exp_coupling_file.indx[element_name]

            f1001 = model_twiss.f1001[model_twiss.indx[element_name]]
            f1010 = model_twiss.f1010[model_twiss.indx[element_name]]
            f1001p = model_play_twiss.f1001[model_play_twiss.indx[element_name]]
            f1010p = model_play_twiss.f1010[model_play_twiss.indx[element_name]]

            print >> sbs_coupling_file, element_name, exp_coupling_file.S[element_index], abs(f1001), f1001.real,\
                    f1001.imag, abs(f1010), f1010.real, f1010.imag, exp_coupling_file.F1001W[element_index],\
                    exp_coupling_file.FWSTD1[element_index], exp_coupling_file.F1001R[element_index],\
                    exp_coupling_file.F1001I[element_index], exp_coupling_file.F1010W[element_index],\
                    exp_coupling_file.FWSTD2[element_index], exp_coupling_file.F1010R[element_index],\
                    exp_coupling_file.F1010I[element_index], abs(f1001p), f1001p.real, f1001p.imag, abs(f1010p),\
                    f1010p.real, f1010p.imag, model_twiss.S[model_twiss.indx[element_name]]
        else:
            print "Element " + element_name + " not found in experimental coupling file"
    sbs_coupling_file.close()


def write_phase_files(output_path, experimental_path, label, method, start_bpm, model_twiss, model_play_twiss):
    exp_phase_twiss_x = _try_to_load_twiss(experimental_path, "getphasetotx" + method + ".out")
    exp_phase_twiss_y = _try_to_load_twiss(experimental_path, "getphasetoty" + method + ".out")
    if not exp_phase_twiss_x or not exp_phase_twiss_y:
        print >> sys.stderr, "Won't write phase files"
        return

    valid_bpms_x = _get_bpms_in_experiment_and_model(exp_phase_twiss_x, model_twiss)
    valid_bpms_y = _get_bpms_in_experiment_and_model(exp_phase_twiss_y, model_twiss)

    if len(valid_bpms_x) == 0 or len(valid_bpms_y) == 0:
        print >> sys.stderr, "No common bmps in experiment and model, won't write phases file"
        return

    if start_bpm not in valid_bpms_x or start_bpm not in valid_bpms_y:
        print >> sys.stderr, "Selected start bpms not in experimental bpms list, won't write phases file"
        return

    sbs_phase_x_file = open(os.path.join(output_path, 'sbsphasext_' + label + '.out'), 'w')

    print >> sbs_phase_x_file, _get_headers_for("phase_x")
    print >> sbs_phase_x_file, _get_data_type_header_for("phase_x")

    for bpm_name in valid_bpms_x:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_exp = exp_phase_twiss_x.indx[bpm_name]

        phase_increase_model = (model_twiss.MUX[bpm_index_model] - model_twiss.MUX[model_twiss.indx[start_bpm]]) % 1

        phase_increase_exp = (exp_phase_twiss_x.PHASEX[bpm_index_exp] -
                              exp_phase_twiss_x.PHASEX[exp_phase_twiss_x.indx[start_bpm]]) % 1

        exp_model_phase_difference = (phase_increase_exp - phase_increase_model) % 1
        if exp_model_phase_difference > 0.5:
            exp_model_phase_difference = exp_model_phase_difference - 1
        play_model_phase_difference = model_play_twiss.MUX[bpm_index_model_play] - model_twiss.MUX[bpm_index_model]

        print >> sbs_phase_x_file,\
                bpm_name, exp_phase_twiss_x.S[bpm_index_exp], phase_increase_exp, exp_model_phase_difference,\
                exp_phase_twiss_x.STDPHX[bpm_index_exp], play_model_phase_difference, model_twiss.S[bpm_index_model]

    sbs_phase_x_file.close()

    sbs_phase_y_file = open(os.path.join(output_path, 'sbsphaseyt_' + label + '.out'), 'w')

    print >> sbs_phase_y_file, _get_headers_for("phase_y")
    print >> sbs_phase_y_file, _get_data_type_header_for("phase_y")

    for bpm_name in valid_bpms_y:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_exp = exp_phase_twiss_y.indx[bpm_name]

        phase_increase_model = (model_twiss.MUY[bpm_index_model] - model_twiss.MUY[model_twiss.indx[start_bpm]]) % 1

        phase_increase_exp = (exp_phase_twiss_y.PHASEY[bpm_index_exp] -
                              exp_phase_twiss_y.PHASEY[exp_phase_twiss_y.indx[start_bpm]]) % 1

        exp_model_phase_difference = (phase_increase_exp - phase_increase_model) % 1
        if exp_model_phase_difference > 0.5:
            exp_model_phase_difference = exp_model_phase_difference - 1
        play_model_phase_difference = model_play_twiss.MUY[bpm_index_model_play] - model_twiss.MUY[bpm_index_model]

        print >> sbs_phase_y_file,\
                bpm_name, exp_phase_twiss_y.S[bpm_index_exp], phase_increase_exp, exp_model_phase_difference,\
                exp_phase_twiss_y.STDPHY[bpm_index_exp], play_model_phase_difference, model_twiss.S[bpm_index_model]


def write_w_function_files(output_path, w_files_path, label, model_twiss, model_play_twiss):
    w_twiss_x = _try_to_load_twiss(w_files_path, "wx.out")
    w_twiss_y = _try_to_load_twiss(w_files_path, "wy.out")
    if not w_twiss_x or not w_twiss_y:
        print "Won't write w function files"
        return

    sbs_w_file_x = open(os.path.join(output_path, "sbsWx_" + label + ".out"), "w")

    print >> sbs_w_file_x, _get_headers_for("w_function_x")
    print >> sbs_w_file_x, _get_data_type_header_for("w_function_x")

    valid_bmps_x = _get_bpms_in_experiment_and_model(w_twiss_x, model_twiss)
    for bpm_name in valid_bmps_x:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_w = w_twiss_x.indx[bpm_name]
        print >> sbs_w_file_x,\
                bpm_name, w_twiss_x.S[bpm_index_w], w_twiss_x.WX[bpm_index_w], w_twiss_x.WXERR[bpm_index_w],\
                model_twiss.WX[bpm_index_model], model_play_twiss.WX[bpm_index_model_play],\
                w_twiss_x.PHIX[bpm_index_w], w_twiss_x.PHIXERR[bpm_index_w], model_twiss.PHIX[bpm_index_model],\
                model_play_twiss.PHIX[bpm_index_model_play], model_twiss.S[bpm_index_model]
    sbs_w_file_x.close()

    sbs_w_file_y = open(os.path.join(output_path, "sbsWy_" + label + ".out"), "w")

    print >> sbs_w_file_y, _get_headers_for("w_function_y")
    print >> sbs_w_file_y, _get_data_type_header_for("w_function_y")

    valid_bmps_y = _get_bpms_in_experiment_and_model(w_twiss_y, model_twiss)
    for bpm_name in valid_bmps_y:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_w = w_twiss_y.indx[bpm_name]
        print >> sbs_w_file_y,\
                bpm_name, w_twiss_y.S[bpm_index_w], w_twiss_y.WY[bpm_index_w], w_twiss_y.WYERR[bpm_index_w],\
                model_twiss.WY[bpm_index_model], model_play_twiss.WY[bpm_index_model_play],\
                w_twiss_y.PHIY[bpm_index_w], w_twiss_y.PHIYERR[bpm_index_w], model_twiss.PHIY[bpm_index_model],\
                model_play_twiss.PHIY[bpm_index_model_play], model_twiss.S[bpm_index_model]
    sbs_w_file_y.close()


def write_dispersion_files(output_path, experimental_path, label, model_twiss, model_play_twiss):
    dispersion_twiss_x = _try_to_load_twiss(experimental_path, "getDx.out")
    dispersion_twiss_y = _try_to_load_twiss(experimental_path, "getDy.out")
    normalized_dispersion_twiss_x = _try_to_load_twiss(experimental_path, "getNDx.out")

    if not dispersion_twiss_x or not dispersion_twiss_y or not normalized_dispersion_twiss_x:
        print "Won't write dispersion files"
        return

    sbs_dispersion_x_file = open(os.path.join(output_path, 'sbsDx_' + label + '.out'), 'w')

    print >> sbs_dispersion_x_file, _get_headers_for("dx")
    print >> sbs_dispersion_x_file, _get_data_type_header_for("dx")

    valid_bpms_x = _get_bpms_in_experiment_and_model(dispersion_twiss_x, model_twiss)

    for bpm_name in valid_bpms_x:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_dispersion = dispersion_twiss_x.indx[bpm_name]
        print >> sbs_dispersion_x_file, bpm_name, dispersion_twiss_x.S[bpm_index_dispersion],\
                dispersion_twiss_x.DX[bpm_index_dispersion], dispersion_twiss_x.STDDX[bpm_index_dispersion],\
                model_twiss.DX[bpm_index_model], model_play_twiss.DX[bpm_index_model_play],\
                model_twiss.S[bpm_index_model]
    sbs_dispersion_x_file.close()

    sbs_dispersion_y_file = open(os.path.join(output_path, 'sbsDy_' + label + '.out'), 'w')

    print >> sbs_dispersion_y_file, _get_headers_for("dy")
    print >> sbs_dispersion_y_file, _get_data_type_header_for("dy")

    valid_bpms_y = _get_bpms_in_experiment_and_model(dispersion_twiss_y, model_twiss)

    for bpm_name in valid_bpms_y:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_dispersion = dispersion_twiss_y.indx[bpm_name]
        print >> sbs_dispersion_y_file, bpm_name, dispersion_twiss_y.S[bpm_index_dispersion],\
                dispersion_twiss_y.DY[bpm_index_dispersion], dispersion_twiss_y.STDDY[bpm_index_dispersion],\
                model_twiss.DY[bpm_index_model], model_play_twiss.DY[bpm_index_model_play],\
                model_twiss.S[bpm_index_model]
    sbs_dispersion_x_file.close()

    sbs_norm_dispersion_x_file = open(os.path.join(output_path, 'sbsNDx_' + label + '.out'), 'w')

    print >> sbs_norm_dispersion_x_file, _get_headers_for("ndx")
    print >> sbs_norm_dispersion_x_file, _get_data_type_header_for("ndx")

    valid_bpms_x = _get_bpms_in_experiment_and_model(normalized_dispersion_twiss_x, model_twiss)

    for bpm_name in valid_bpms_x:
        bpm_index_model = model_twiss.indx[bpm_name]
        bpm_index_model_play = model_play_twiss.indx[bpm_name]
        bpm_index_dispersion = normalized_dispersion_twiss_x.indx[bpm_name]
        print >> sbs_dispersion_x_file, bpm_name, normalized_dispersion_twiss_x.S[bpm_index_dispersion],\
                normalized_dispersion_twiss_x.NDX[bpm_index_dispersion],\
                normalized_dispersion_twiss_x.STDNDX[bpm_index_dispersion],\
                model_twiss.DX[bpm_index_model] / math.sqrt(model_twiss.BETX[bpm_index_model]),\
                model_play_twiss.DX[bpm_index_model_play] / math.sqrt(model_play_twiss.BETX[bpm_index_model_play]),\
                model_twiss.S[bpm_index_model]

    sbs_dispersion_x_file.close()


def _try_to_load_twiss(twiss_path, twiss_file_name):
    twiss_full_path = os.path.join(twiss_path, twiss_file_name)
    if not os.path.isfile(twiss_full_path):
        print "Twiss file " + twiss_file_name + " doesn't exist"
        return False
    try:
        stderr = sys.stderr
        sys.stderr = sys.stdout  # Temporary redirect of standard error output to avoid annoying messages from twiss
        twiss_file = metaclass.twiss(twiss_full_path)
        sys.stderr = stderr
    except:
        print "Can't read " + twiss_file_name + " twiss file"
        return False

    return twiss_file


def _get_bpms_in_experiment_and_model(experimental_twiss, model_twiss):
        common_bpms_list = []
        for bpm_name in model_twiss.NAME:
            if bpm_name.upper() in experimental_twiss.indx:
                common_bpms_list.append((experimental_twiss.S[experimental_twiss.indx[bpm_name]], bpm_name))
            else:
                print bpm_name, " in model but not in experiment"
        common_bpms_list.sort()
        return zip(*common_bpms_list)[1]


def _get_headers_for(file_type):
    return {"coupling": '* NAME   S   f1001 f1001re  f1001im    f1010   f1010re   f1010im  f1001_EXP f1001err_EXP  f1001re_EXP  f1001im_EXP    f1010_EXP  f1010err_EXP   f1010re_EXP   f1010im_EXP     f1001_PLAY   f1001re_PLAY  f1001im_PLAY    f1010_PLAY   f1010re_PLAY   f1010im_PLAY  S_MODEL',
            "phase_x": "* NAME  S PHASEX  PHASEXT    ERRORX PHASE_PLAY MODEL_S ",
            "phase_y": "* NAME  S PHASEY  PHASEYT    ERRORY PHASE_PLAY MODEL_S ",
            "w_function_x": "* NAME  S  WX   WXERR WX_MDL WX_PLAY  PHIX PHIXERR PHIX_MDL PHIX_PLA    MODEL_S ",
            "w_function_y": "* NAME  S  WY   WYERR WY_MDL WY_PLAY  PHIY PHIYERR PHIY_MDL PHIY_PLAY      MODEL_S ",
            "dx": "* NAME  S  DX   STDDX DX_MDL DX_PLAY MODEL_S ",
            "dy": "* NAME  S  DY   STDDY DY_MDL DY_PLAY MODEL_S ",
            "ndx": "* NAME  S  NDX   STDNDX NDX_MDL NDX_PLAY MODEL_S "
    }[file_type]


def _get_data_type_header_for(file_type):
    return {"coupling": '$  %s    %le   %le    %le  %le   %le    %le    %le  %le   %le    %le  %le   %le    %le    %le   %le  %le   %le    %le    %le   %le  %le %le',
            "phase_x": "$ %s    %le  %le  %le        %le     %le       %le  ",
            "phase_y": "$ %s   %le  %le   %le        %le    %le     %le ",
            "w_function_x": "$ %s   %le  %le   %le  %le    %le     %le %le  %le    %le     %le ",
            "w_function_y": "$ %s   %le  %le   %le  %le    %le     %le  %le %le %le %le ",
            "dx": "$ %s   %le  %le   %le  %le    %le     %le ",
            "dy": "$ %s   %le  %le   %le  %le    %le     %le ",
            "ndx": "$ %s   %le  %le   %le   %le    %le     %le "
    }[file_type]

if __name__ == "__main__":
    parse_args()
