'''
Created on 11/09/09

@author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)
'''


from __future__ import print_function
import os
import sys
import argparse

from math import sqrt
import json
import numpy

from os.path import abspath, join, dirname, pardir
new_path = abspath(join(dirname(__file__), pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

import utils.iotools
from model import manager, creator
from model.accelerators.lhc import Lhc
from model.accelerators.psbooster import Psbooster
from model.accelerators.accelerator import Element
from Python_Classes4MAD.metaclass import twiss
import madx_wrapper

from SegmentBySegment.sbs_writers import (
    sbs_beta_writer,
    sbs_beta_beating_writer,
    sbs_phase_writer,
    sbs_dispersion_writer,
    sbs_coupling_writer,
    sbs_chromatic_writer,
    sbs_special_element_writer,
)

import logging


#===================================================================================================
# parse_args()-function
#===================================================================================================
def _parse_args(args=None):
    '''
    Parses the arguments, checks for valid input and returns tupel
    It needs also the input needed to define the accelerator:
    --accel=<accel_name>
    and all the rest of the parameters needed to define given accelerator.
    e.g. for LHC runII 2017 beam 1
    --accel lhc --lhcmode lhc_runII_2017 --beam 1
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--path",  # assumes that output is same as input
                        help="Path to measurement files",
                        default="./", dest="path")
    parser.add_argument("-s", "--start",
                        help="give start,endbpm,name (multiple allowed) eg: start1,end1,name1,start2,end2,name2,...",
                        metavar="SEGF", default="./", dest="segf")
    parser.add_argument("-t", "--twiss",
                        help="basic twiss file, the modifiers.madx is assumed to be in the same directory",
                        metavar="TWISS", default="./", dest="twiss")
    parser.add_argument("-p", "--save",
                        help="Output path",
                        metavar="SAVE", default="./", dest="save")
    parser.add_argument("-m", "--mad",  # assumes that output is same as input
                        help="mad link",
                        default="", dest="mad")
    parser.add_argument("-b", "--bbsource",  # assumes that output is same as input
                        help="beta beat source",
                        metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="bb")
    parser.add_argument("-x", "--take",
                        help="If present, it will run MADX on previously "
                             "created scripts, without remaking them.",
                        dest="madpass", action="store_true")
    parser.add_argument("-c", "--cuts",
                        help="cut on error of beta in percentage",
                        metavar="cuts", default="10", dest="cuts")
    parser.add_argument("-w", "--w",  # Path to Chromaticity functions
                        help="Path to  chromaticity functions, by default this is skiped",
                        metavar="wpath", default="0", dest="wpath")
    options, accel_args = parser.parse_known_args(args)

    accel_cls = manager.get_accel_class(accel_args)
    return accel_cls, options


#==============================================================================
# main()-function
#==============================================================================


def main(accel_cls, options):
    '''
    :Parameters:
        'options': Values
            Values instance with all options from OptionParser
    :Return: int
        0 if execution was successful otherwise !=0
    '''
    if (sys.flags.debug):
        logging.basicConfig(level=logging.INFO)
        
    print("+++ Starting Segment by Segment +++")
    print("Using accelerator class: " + accel_cls.__name__)
    measurement_path = options.path
    w_path = options.wpath
    if w_path == "0":
        w_path = measurement_path
    input_data = _InputData(measurement_path, w_path)

    save_path = options.save + os.path.sep
    utils.iotools.create_dirs(save_path)

    elements_data = options.segf.split(',')
    error_cut = float(options.cuts)

    twiss_file = options.twiss
    print("Input model twiss file", twiss_file)
    input_model = _try_to_load_twiss(twiss_file)
    if input_model is None:
        raise IOError("Cannot read input model, aborting.")

    twiss_directory = os.path.dirname(twiss_file)

    elements_names, start_bpms, end_bpms = structure_elements_info(elements_data)

    summaries = _Summaries(save_path)

    for element_name in elements_names:

        print("Started processing", element_name)

        start_bpm_name, end_bpm_name, is_element = get_good_bpms(input_data, error_cut, input_model, start_bpms, end_bpms, element_name)
        if issubclass(accel_cls, Psbooster):
            (start_bpm_horizontal_data,
            start_bpm_vertical_data,
            end_bpm_horizontal_data,
            end_bpm_vertical_data) = gather_data_amplitude(input_data, start_bpm_name, end_bpm_name)
        else:
            (start_bpm_horizontal_data,
            start_bpm_vertical_data,
            end_bpm_horizontal_data,
            end_bpm_vertical_data) = gather_data(input_data, start_bpm_name, end_bpm_name)
        element_has_dispersion, start_bpm_dispersion, end_bpm_dispersion = _get_dispersion_parameters(input_data, start_bpm_name, end_bpm_name)

        element_has_coupling, f_ini, f_end = _get_coupling_parameters(input_data, start_bpm_name, end_bpm_name)

        element_has_chrom, chrom_ini, chrom_end = _get_chrom_parameters(input_data, start_bpm_name, end_bpm_name)

        accel_instance = accel_cls.get_segment(
            element_name,
            start_bpm_name,
            end_bpm_name,
            os.path.join(save_path, "modifiers.madx"),
            twiss_file
        )

        if not options.madpass:
            _run4mad(save_path,
                     accel_instance,
                     start_bpm_horizontal_data,
                     start_bpm_vertical_data,
                     end_bpm_horizontal_data,
                     end_bpm_vertical_data,
                     start_bpm_dispersion,
                     end_bpm_dispersion,
                     f_ini,
                     f_end,
                     chrom_ini,
                     chrom_end,
                     options.path,
                     twiss_directory,
                     input_data.couple_method,
                     options.bb,
                     options.mad)

        else:
            print("Just rerunning mad")
            mad_file_path, log_file_path = _get_files_for_mad(save_path,
                                                              element_name)
            madx_wrapper.resolve_and_run_file(mad_file_path,
                                              log_file=log_file_path)

        propagated_models = _PropagatedModels(save_path, element_name)

        kmod_data_file_x, kmod_data_file_y = _get_kmod_files()

        getAndWriteData(element_name,
                        input_data,
                        input_model,
                        propagated_models,
                        save_path,
                        is_element,
                        element_has_dispersion,
                        element_has_coupling,
                        element_has_chrom,
                        accel_instance,
                        summaries,
                        kmod_data_file_x,
                        kmod_data_file_y)
        print("Everything done for", element_name, "\n")

    summaries.write_summaries_to_files()

    print("+++  Ended Segment by Segment   +++")

# END main() ---------------------------------------------------------------------------------------


def structure_elements_info(elements_data):
    start_bpms = {}
    end_bpms = {}
    elements_names = []

    step = 0
    for position in range(len(elements_data)):
        position = position + step
        if position == (len(elements_data)):
            break

        if "BPM" in elements_data[position]:
            segment_name = elements_data[position + 2]
            start_bpms[segment_name] = elements_data[position]
            end_bpms[segment_name] = elements_data[position + 1]
            elements_names.append(segment_name)
            step = 2 + step

        else:
            step = 0 + step
            elements_names.append(elements_data[position])

    return elements_names, start_bpms, end_bpms


def get_good_bpms(input_data, errorcut, twiss_data, start_bpms, end_bpms, element_name):
    #  checking if end_bpms element is BPM or instrument (IP,collimators)
    if element_name in start_bpms and element_name in end_bpms:
        is_element = False
        segment = [start_bpms[element_name], end_bpms[element_name]]
        start_bpm_name, end_bpm_name = _filter_and_find(input_data.beta_x,
                                                        input_data.beta_y,
                                                        "null", segment,
                                                        twiss_data, errorcut)
    elif element_name in start_bpms or element_name in end_bpms:
        raise SegmentBySegmentError(
            "Something strange ....Did you properly define the input ?\n"
            "Like: BPM1,BPM2,ARC12,IP1\n"
            "Structure must be for Segment => BPML,BPMR,NAME\n"
            "For Instrument just name\n"
        )
    else:
        is_element = True
        start_bpm_name, end_bpm_name = _filter_and_find(input_data.beta_x,
                                                        input_data.beta_y,
                                                        element_name,
                                                        [],
                                                        twiss_data,
                                                        errorcut)
    return start_bpm_name, end_bpm_name, is_element


def gather_data(input_data, startbpm, endbpm):
    start_bpm_horizontal_data = [input_data.beta_x.BETX[input_data.beta_x.indx[startbpm]],
                                 input_data.beta_x.ALFX[input_data.beta_x.indx[startbpm]]]
    start_bpm_vertical_data = [input_data.beta_y.BETY[input_data.beta_y.indx[startbpm]],
                               input_data.beta_y.ALFY[input_data.beta_y.indx[startbpm]]]
    end_bpm_horizontal_data = [input_data.beta_x.BETX[input_data.beta_x.indx[endbpm]],
                               input_data.beta_x.ALFX[input_data.beta_x.indx[endbpm]]]
    end_bpm_vertical_data = [input_data.beta_y.BETY[input_data.beta_y.indx[endbpm]],
                             input_data.beta_y.ALFY[input_data.beta_y.indx[endbpm]]]
    return start_bpm_horizontal_data, start_bpm_vertical_data, end_bpm_horizontal_data, end_bpm_vertical_data


def gather_data_amplitude(input_data, startbpm, endbpm):
    start_bpm_horizontal_data = [input_data.amplitude_beta_x.BETX[input_data.amplitude_beta_x.indx[startbpm]],
                                 input_data.beta_x.ALFX[input_data.beta_x.indx[startbpm]]]
    start_bpm_vertical_data = [input_data.amplitude_beta_y.BETY[input_data.amplitude_beta_y.indx[startbpm]],
                               input_data.beta_y.ALFY[input_data.beta_y.indx[startbpm]]]
    end_bpm_horizontal_data = [input_data.amplitude_beta_x.BETX[input_data.amplitude_beta_x.indx[endbpm]],
                               input_data.beta_x.ALFX[input_data.beta_x.indx[endbpm]]]
    end_bpm_vertical_data = [input_data.amplitude_beta_y.BETY[input_data.amplitude_beta_y.indx[endbpm]],
                             input_data.beta_y.ALFY[input_data.beta_y.indx[endbpm]]]
    return start_bpm_horizontal_data, start_bpm_vertical_data, end_bpm_horizontal_data, end_bpm_vertical_data


def _get_dispersion_parameters(input_data, startbpm, endbpm):
    start_bpm_dispersion = [0, 0, 0, 0]
    end_bpm_dispersion = [0, 0, 0, 0]
    element_has_dispersion = False
    if input_data.has_dispersion:
        if startbpm not in input_data.dispersion_x.indx:
            print("Start BPM ", startbpm, " not found in horizontal dispersion measurement, will not compute dispersion.")
        elif startbpm not in input_data.dispersion_y.indx:
            print("Start BPM ", startbpm, " not found in vertical dispersion measurement, will not compute dispersion.")
        elif endbpm not in input_data.dispersion_x.indx:
            print("End BPM ", endbpm, " not found in horizontal dispersion measurement, will not compute dispersion.")
        elif endbpm not in input_data.dispersion_y.indx:
            print("End BPM ", endbpm, " not found in vertical dispersion measurement, will not compute dispersion.")
        else:
            dxx_start = input_data.dispersion_x.DX[input_data.dispersion_x.indx[startbpm]]
            dxp_start = input_data.dispersion_x.DPX[input_data.dispersion_x.indx[startbpm]]
            dyy_start = input_data.dispersion_y.DY[input_data.dispersion_y.indx[startbpm]]
            dyp_start = input_data.dispersion_y.DPY[input_data.dispersion_y.indx[startbpm]]
            start_bpm_dispersion = [dxx_start, dxp_start, dyy_start, dyp_start]

            dxx_end = input_data.dispersion_x.DX[input_data.dispersion_x.indx[endbpm]]
            dxp_end = input_data.dispersion_x.DPX[input_data.dispersion_x.indx[endbpm]]
            dyy_end = input_data.dispersion_y.DY[input_data.dispersion_y.indx[endbpm]]
            dyp_end = input_data.dispersion_y.DPY[input_data.dispersion_y.indx[endbpm]]
            end_bpm_dispersion = [dxx_end, dxp_end, dyy_end, dyp_end]
            element_has_dispersion = True
            print("Start and end BPMs found in dispersion measurement.")
    return element_has_dispersion, start_bpm_dispersion, end_bpm_dispersion


def _get_coupling_parameters(input_data, startbpm, endbpm):
    f_ini = {}
    f_end = {}
    f_ini["f1001r"] = 0.
    f_ini["f1001i"] = 0.
    f_ini["f1010r"] = 0.
    f_ini["f1010i"] = 0.
    f_ini["f1001std"] = 0.
    f_ini["f1010std"] = 0.
    f_end["f1001r"] = 0.
    f_end["f1001i"] = 0.
    f_end["f1010r"] = 0.
    f_end["f1010i"] = 0.
    f_end["f1001std"] = 0.
    f_end["f1010std"] = 0.
    element_has_coupling = False
    if input_data.has_coupling:
        if startbpm not in input_data.couple.indx:
            print("Start BPM ", startbpm, " not found in coupling measurement, will not compute coupling.")
        elif not endbpm in input_data.couple.indx:
            print("End BPM ", endbpm, " not found in coupling measurement, will not compute coupling.")
        else:
            f_ini["f1001r"] = input_data.couple.F1001R[input_data.couple.indx[startbpm]]
            f_ini["f1001i"] = input_data.couple.F1001I[input_data.couple.indx[startbpm]]
            f_ini["f1010r"] = input_data.couple.F1010R[input_data.couple.indx[startbpm]]
            f_ini["f1010i"] = input_data.couple.F1010I[input_data.couple.indx[startbpm]]
            f_ini["f1001std"] = input_data.couple.FWSTD1[input_data.couple.indx[startbpm]]
            f_ini["f1010std"] = input_data.couple.FWSTD2[input_data.couple.indx[startbpm]]
            f_end["f1001r"] = input_data.couple.F1001R[input_data.couple.indx[endbpm]]
            f_end["f1001i"] = input_data.couple.F1001I[input_data.couple.indx[endbpm]]
            f_end["f1010r"] = input_data.couple.F1010R[input_data.couple.indx[endbpm]]
            f_end["f1010i"] = input_data.couple.F1010I[input_data.couple.indx[endbpm]]
            f_end["f1001std"] = input_data.couple.FWSTD1[input_data.couple.indx[endbpm]]
            f_end["f1010std"] = input_data.couple.FWSTD2[input_data.couple.indx[endbpm]]
            element_has_coupling = True
            print("Start and end BPMs found in coupling measurement.")
    else:
        print("No coupling measurement")
    return element_has_coupling, f_ini, f_end


def _get_chrom_parameters(input_data, startbpm, endbpm):
    chrom_ini = {}
    chrom_end = {}
    chrom_ini["wx"] = 0.
    chrom_ini["wy"] = 0.
    chrom_ini["phi_x"] = 0.
    chrom_ini["phi_y"] = 0.
    chrom_end["wx"] = 0.
    chrom_end["wy"] = 0.
    chrom_end["phi_x"] = 0.
    chrom_end["phi_y"] = 0.
    element_has_chrom = False
    if input_data.has_chromatic:
        if not startbpm in input_data.wx.indx:
            print("Start BPM ", startbpm, " not found in chromatic measurement, will not compute chromatic.")
        elif not endbpm in input_data.wy.indx:
            print("End BPM ", endbpm, " not found in chromatic measurement, will not compute chromatic.")
        else:
            chrom_ini["wx"] = input_data.wx.WX[input_data.wx.indx[startbpm]]
            chrom_ini["wy"] = input_data.wy.WY[input_data.wy.indx[startbpm]]
            chrom_ini["phi_x"] = input_data.wx.PHIX[input_data.wx.indx[startbpm]]
            chrom_ini["phi_y"] = input_data.wy.PHIY[input_data.wy.indx[startbpm]]
            chrom_end["wx"] = input_data.wx.WX[input_data.wx.indx[endbpm]]
            chrom_end["wy"] = input_data.wy.WY[input_data.wy.indx[endbpm]]
            chrom_end["phi_x"] = input_data.wx.PHIX[input_data.wx.indx[endbpm]]
            chrom_end["phi_y"] = input_data.wy.PHIY[input_data.wy.indx[endbpm]]
            element_has_chrom = True
            print("Start and end BPMs found in chromatic measurement.")
    return element_has_chrom, chrom_ini, chrom_end


def _get_kmod_files():
    try:
        kmod_path_x = twiss(_join_output_with("getkmodbetax.out"))
        kmod_path_y = twiss(_join_output_with("getkmodbetay.out"))
    except IOError:
        kmod_path_x = None
        kmod_path_y = None
    return kmod_path_x, kmod_path_y


def _get_calibrated_betas(plane):
    calibration_data = None
    try:
        amplitude_beta = _get_twiss_for_one_of("getampbeta{}_free.out".format(plane),
                                               "getampbeta{}.out".format(plane))
    except IOError:
        return None
    try:
        calibration_data = twiss(_join_output_with("calibration_" + plane + ".out"))
    except IOError:
        return amplitude_beta
    calibrated_betas = CalibratedBetas(plane)
    index_counter = 0
    for bpm_name in calibration_data.NAME:
        if bpm_name in amplitude_beta.NAME:
            amp_index = amplitude_beta.indx[bpm_name]
            amp_beta = getattr(amplitude_beta, "BET" + plane.upper())[amp_index]
            calibrated_betas.NAME.append(bpm_name)
            calibrated_betas.S.append(amplitude_beta.S[amp_index])
            calibrated_betas.indx[bpm_name] = index_counter
            index_counter += 1
            getattr(calibrated_betas, "BET" + plane.upper()).append(amp_beta)
            try:
                err_amp_beta = getattr(amplitude_beta, "STDBET" + plane.upper())[amp_index]
            except AttributeError:
                err_amp_beta = getattr(amplitude_beta, "ERRBET" + plane.upper())[amp_index]
            getattr(calibrated_betas, "BET" + plane.upper() + "STD").append(err_amp_beta)
            mdl_amp_beta = getattr(amplitude_beta, "BET" + plane.upper() + "MDL")[amp_index]
            getattr(calibrated_betas, "BET" + plane.upper() + "MDL").append(mdl_amp_beta)
    return calibrated_betas


class CalibratedBetas(object):
    def __init__(self, plane):
        self.NAME = []
        self.S = []
        setattr(self, "BET" + plane.upper(), [])
        setattr(self, "BET" + plane.upper() + "STD", [])
        setattr(self, "BET" + plane.upper() + "MDL", [])
        self.indx = {}


#===================================================================================================
# helper-functions
#===================================================================================================


def _filter_and_find(beta_x_twiss, beta_y_twiss, element_name, segment_bpms_names, model, errorcut):
    '''
    Automatic BPM filter and BPM finder

    :Parameters:
        'beta_x_twiss': twiss
            twiss for horizontal beta function
        'beta_y_twiss': twiss
            twiss for vertcal beta function
        'element_name':
            instrument,collimator or IP . Specify as "null" if it as a segment
        'segment_bpms_names':
            leftbpm and rightbpm where you want to compute the segment_bpms_names
        'model': twiss
            twiss model
        'errorcut': float
            errorcut for beta in percentage (err/bet)
    :Return: list
        Left and right BPM in a list
    :Exception:
        If n BPMS is  < 3 => System exit
    '''
    # initial points
    translate = {}  # this dict associates the position in the ring with the element name and if it's valid or not
    locations_list = []

    if element_name != "null":
        is_segment = False
        if element_name in model.indx:
            element_s = model.S[model.indx[element_name]]
        else:
            raise SegmentBySegmentError(element_name + " not found in model.")
        locations_list.append(element_s)
        translate[element_s] = [True, element_name]
        print("You selected an element")
    else:
        is_segment = True
        left_bpm_name = segment_bpms_names[0]
        right_bpm_name = segment_bpms_names[1]
        if left_bpm_name in model.indx and right_bpm_name in model.indx:
            left_bpm_s = model.S[model.indx[left_bpm_name]]
            right_bpm_s = model.S[model.indx[right_bpm_name]]
        else:
            raise SegmentBySegmentError(
                left_bpm_name + " " + right_bpm_name + " not found in model."
            )
        locations_list.append(left_bpm_s)
        locations_list.append(right_bpm_s)

        if (left_bpm_name in beta_y_twiss.indx and
            left_bpm_name in beta_x_twiss.indx and
            right_bpm_name in beta_y_twiss.indx and
            right_bpm_name in beta_x_twiss.indx):
            translate[left_bpm_s] = [True, left_bpm_name]
            translate[right_bpm_s] = [True, right_bpm_name]
        else:
            translate[left_bpm_s] = [False, left_bpm_name]
            translate[right_bpm_s] = [False, right_bpm_name]

        print("You selected a segment")

    elements_names_in_model = model.NAME

    number_of_good_bpms = 0

    # filtering
    for current_element_name in elements_names_in_model:

        if "BPM" in current_element_name:
            current_element_s = model.S[model.indx[current_element_name]]

            if current_element_name in beta_y_twiss.indx and current_element_name in beta_x_twiss.indx:
                beta_y = beta_y_twiss.BETY[beta_y_twiss.indx[current_element_name]]
                beta_x = beta_x_twiss.BETX[beta_x_twiss.indx[current_element_name]]
                err_beta_x = beta_x_twiss.ERRBETX[beta_x_twiss.indx[current_element_name]]
                err_beta_y = beta_y_twiss.ERRBETY[beta_y_twiss.indx[current_element_name]]

                stdbetax_exists = True
                try:
                    beta_x_twiss.STDBETX[beta_x_twiss.indx[current_element_name]]
                except AttributeError:
                    stdbetax_exists = False
                if stdbetax_exists:
                    stdbetax = beta_x_twiss.STDBETX[beta_x_twiss.indx[current_element_name]]
                    total_err_x = sqrt(err_beta_x ** 2 + stdbetax ** 2)
                else:
                    total_err_x = err_beta_x
                stdbetay_exists = True
                try:
                    beta_y_twiss.STDBETY[beta_y_twiss.indx[current_element_name]]
                except AttributeError:
                    stdbetay_exists = False
                if stdbetay_exists:
                    stdbetay = beta_y_twiss.STDBETY[beta_y_twiss.indx[current_element_name]]
                    total_err_y = sqrt(err_beta_y ** 2 + stdbetay ** 2)
                else:
                    total_err_y = err_beta_y

                if (beta_x > 0 and
                        beta_y > 0 and
                        beta_x > total_err_x and
                        beta_y > total_err_y and
                        total_err_x > 0 and
                        total_err_y > 0 and
                        ((total_err_x / beta_x) * 100) < errorcut and
                        ((total_err_y / beta_y) * 100) < errorcut):

                    translate[current_element_s] = [True, current_element_name]
                    number_of_good_bpms = number_of_good_bpms + 1
                    locations_list.append(current_element_s)

                elif is_segment and (left_bpm_name in current_element_name or right_bpm_name in current_element_name):
                    translate[current_element_s] = [False, current_element_name]

    locations_list.sort()
    if number_of_good_bpms < 3:
        raise SegmentBySegmentError(
            "Not enough BPMs! Less than 3 BPMs remaining after filtering!"
        )

    # finding the BPMs
    if not is_segment:
        element_location_index = locations_list.index(element_s)

        if element_location_index == 0:
            selected_left_bpm = translate[locations_list[len(locations_list) - 2]][1]
            selected_right_bpm = translate[locations_list[1]][1]
        elif element_location_index == (len(locations_list) - 1):
            selected_left_bpm = translate[locations_list[element_location_index - 1]][1]
            selected_right_bpm = translate[locations_list[0]][1]
        else:
            selected_left_bpm = translate[locations_list[element_location_index - 1]][1]
            selected_right_bpm = translate[locations_list[element_location_index + 1]][1]

    else:
        left_bpm_is_good = translate[left_bpm_s][0]
        right_bpm_is_good = translate[right_bpm_s][0]

        if left_bpm_is_good:
            selected_left_bpm = translate[left_bpm_s][1]
        else:
            element_location_index = locations_list.index(left_bpm_s)
            if element_location_index == 0:
                selected_left_bpm = translate[locations_list[len(locations_list) - 2]][1]
            else:
                selected_left_bpm = translate[locations_list[element_location_index - 1]][1]

        if right_bpm_is_good:
            selected_right_bpm = translate[right_bpm_s][1]
        else:
            element_location_index = locations_list.index(right_bpm_s)
            if element_location_index == (len(locations_list) - 1):
                selected_right_bpm = translate[locations_list[0]][1]
            else:
                selected_right_bpm = translate[locations_list[element_location_index + 1]][1]

    print(selected_left_bpm, selected_right_bpm, "Will be used for the propagation")

    return [selected_left_bpm, selected_right_bpm]


def getAndWriteData(
    element_name, input_data, input_model, propagated_models,
    save_path, is_element,
    element_has_dispersion, element_has_coupling, element_has_chrom,
    accel_instance, summaries,
    kmod_data_x, kmod_data_y
):
    '''
    Function that returns the optics function at the given element

    :Parameters:
        # TODO: rewrite this
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''

    print("Start writing files for", element_name)

    if hasattr(summaries, 'beta'):
        beta_summary = summaries.beta
        disp_summary = summaries.dispersion
        coupling_summary = summaries.coupling
        chrom_summary = summaries.chrom
    else:
        beta_summary = None
        disp_summary = None
        coupling_summary = None
        chrom_summary = None

    (beta_x, err_beta_x, alfa_x, err_alfa_x,
     beta_y, err_beta_y, alfa_y, err_alfa_y) = sbs_beta_writer.write_beta(
        element_name, is_element,
        input_data.beta_x, input_data.beta_y,
        input_model,
        propagated_models,
        save_path, beta_summary
    )
    if not is_element:
        sbs_phase_writer.write_phase(
            element_name,
            input_data.total_phase_x, input_data.total_phase_y,
            input_data.beta_x, input_data.beta_y,
            propagated_models, save_path
        )
        sbs_beta_beating_writer.write_beta_beat(
            element_name,
            input_data.beta_x, input_data.beta_y,
            input_data.amplitude_beta_x, input_data.amplitude_beta_y,
            kmod_data_x, kmod_data_y,
            propagated_models, save_path)
    if element_has_dispersion:
        sbs_dispersion_writer.write_dispersion(element_name, is_element,
                                                           input_data.dispersion_x, input_data.dispersion_y, input_data.normalized_dispersion_x,
                                                           input_model,
                                                           propagated_models, save_path, disp_summary)
    if element_has_coupling:
        sbs_coupling_writer.write_coupling(element_name, is_element, input_data.couple, input_model, propagated_models, save_path, coupling_summary)

    if element_has_chrom:
        sbs_chromatic_writer.write_chromatic(element_name, is_element, input_data.wx, input_data.wy, input_model, propagated_models, save_path, chrom_summary)

    if "IP" in element_name and is_element:
        sbs_special_element_writer.write_ip(input_data.beta_x, input_data.beta_y,
                                                        beta_x, err_beta_x, alfa_x, err_alfa_x,
                                                        beta_y, err_beta_y, alfa_y, err_alfa_y,
                                                        input_model, input_data.phase_x, input_data.phase_y, element_name,
                                                        accel_instance, save_path)
    # TODO: This need to be fixed before using (find BPMs)
    elif False and "ADT" in element_name and is_element:  # if False to avoid going inside
        sbs_special_element_writer.write_transverse_damper(propagated_models, element_name, input_model,
                                                                       save_path, input_data.phase_x, input_data.phase_y,
                                                                       input_data.beta_x, input_data.beta_y,
                                                                       accel_instance)


def _run4mad(save_path,
             accel_instance,
             start_bpm_horizontal_data,
             start_bpm_vertical_data,
             end_bpm_horizontal_data,
             end_bpm_vertical_data,
             start_bpm_dispersion,
             end_bpm_dispersion,
             f_ini,
             f_end,
             chrom_ini,
             chrom_end,
             exppath,
             twiss_directory,
             coupling_method,
             bb_path,
             madx_exe_path):
    _copy_modifiers_and_corrections_locally(save_path, twiss_directory,
                                            accel_instance)

    betx_ini = start_bpm_horizontal_data[0]
    bety_ini = start_bpm_vertical_data[0]
    betx_end = end_bpm_horizontal_data[0]
    bety_end = end_bpm_vertical_data[0]

    if betx_ini < 0. or bety_ini < 0.:
        raise SegmentBySegmentError(
            "Negative betas in initial BPM of segment!"
        )
    if betx_end < 0. or bety_end < 0.:
        raise SegmentBySegmentError("Negative betas in last BPM of segment!")

    alfx_ini = start_bpm_horizontal_data[1]
    alfy_ini = start_bpm_vertical_data[1]
    alfx_end = -end_bpm_horizontal_data[1]
    alfy_end = -end_bpm_vertical_data[1]

    ini_r11, ini_r12, ini_r21, ini_r22 = _get_R_terms(
        betx_ini, bety_ini, alfx_ini, alfy_ini,
        f_ini["f1001r"], f_ini["f1001i"],
        f_ini["f1010r"], f_ini["f1010i"]
    )
    end_r11, end_r12, end_r21, end_r22 = _get_R_terms(
        betx_end, bety_end, alfx_end, alfy_end,
        f_end["f1001r"], f_end["f1001i"],
        f_end["f1010r"], f_end["f1010i"]
    )

    mad_file_path, log_file_path = _get_files_for_mad(
        save_path, accel_instance.label
    )

    accel_instance.start = Element(accel_instance.start.name.replace("-", "_"),
                                   accel_instance.start.s)
    accel_instance.end = Element(accel_instance.end.name.replace("-", "_"),
                                 accel_instance.end.s)

    measurement_dict = dict(
        betx_ini=betx_ini,
        bety_ini=bety_ini,
        alfx_ini=alfx_ini,
        alfy_ini=alfy_ini,
        dx_ini=start_bpm_dispersion[0],
        dy_ini=start_bpm_dispersion[2],
        dpx_ini=start_bpm_dispersion[1],
        dpy_ini=start_bpm_dispersion[3],
        wx_ini=chrom_ini["wx"],
        phix_ini=chrom_ini["phi_x"],
        wy_ini=chrom_ini["wy"],
        phiy_ini=chrom_ini["phi_y"],
        wx_end=chrom_end["wx"],
        phix_end=chrom_end["phi_x"],
        wy_end=chrom_end["wy"],
        phiy_end=chrom_end["phi_y"],
        ini_r11=ini_r11,
        ini_r12=ini_r12,
        ini_r21=ini_r21,
        ini_r22=ini_r22,
        end_r11=end_r11,
        end_r12=end_r12,
        end_r21=end_r21,
        end_r22=end_r22,
        betx_end=betx_end,
        bety_end=bety_end,
        alfx_end=alfx_end,
        alfy_end=alfy_end,
        dx_end=end_bpm_dispersion[0],
        dy_end=end_bpm_dispersion[2],
        dpx_end=-end_bpm_dispersion[1],
        dpy_end=-end_bpm_dispersion[3],
    )

    with open(os.path.join(save_path, "measurement_" + accel_instance.label + ".madx"), "w") as measurement_file:
        for name, value in measurement_dict.iteritems():
            measurement_file.write(name + " = " + str(value) + ";\n")

    _runmad(accel_instance, save_path, mad_file_path, log_file_path)


def _get_R_terms(betx, bety, alfx, alfy, f1001r, f1001i, f1010r, f1010i):

    ga11 = 1 / numpy.sqrt(betx)
    ga12 = 0
    ga21 = alfx / numpy.sqrt(betx)
    ga22 = numpy.sqrt(betx)
    Ga = numpy.reshape(numpy.array([ga11, ga12, ga21, ga22]), (2, 2))

    gb11 = 1 / numpy.sqrt(bety)
    gb12 = 0
    gb21 = alfy / numpy.sqrt(bety)
    gb22 = numpy.sqrt(bety)
    Gb = numpy.reshape(numpy.array([gb11, gb12, gb21, gb22]), (2, 2))

    J = numpy.reshape(numpy.array([0, 1, -1, 0]), (2, 2))

    absf1001 = numpy.sqrt(f1001r ** 2 + f1001i ** 2)
    absf1010 = numpy.sqrt(f1010r ** 2 + f1010i ** 2)

    gamma2 = 1. / (1. + 4. * (absf1001 ** 2 - absf1010 ** 2))
    c11 = (f1001i + f1010i)
    c22 = (f1001i - f1010i)
    c12 = -(f1010r - f1001r)
    c21 = -(f1010r + f1001r)
    Cbar = numpy.reshape(2 * numpy.sqrt(gamma2) *
                         numpy.array([c11, c12, c21, c22]), (2, 2))

    C = numpy.dot(numpy.linalg.inv(Ga), numpy.dot(Cbar, Gb))
    jCj = numpy.dot(J, numpy.dot(C, -J))
    c = numpy.linalg.det(C)
    r = -c / (c - 1)
    R = numpy.transpose(numpy.sqrt(1 + r) * jCj)
    return numpy.ravel(R)


def _get_files_for_mad(save_path, element_name):
    mad_file_name = 't_' + str(element_name) + '.madx'
    log_file_name = element_name + "_mad.log"
    mad_file_path = os.path.join(save_path, mad_file_name)
    log_file_path = os.path.join(save_path, log_file_name)
    return mad_file_path, log_file_path


def _copy_modifiers_and_corrections_locally(save_path, twiss_directory, accel_instance):
    modifiers_file_path = os.path.join(twiss_directory, 'modifiers.madx')
    if os.path.isfile(modifiers_file_path):
        utils.iotools.copy_item(modifiers_file_path, save_path)
    else:
        print("Cannot find modifiers.madx file, will create an empty file.")
        open(os.path.join(save_path, 'modifiers.madx'), "a").close()

    correction_file_comments = ""
    if (accel_instance.label.lower().startswith("ip") and
            issubclass(type(accel_instance), Lhc)):
        correction_file_comments = _get_corrections_file_comments_for_ip(accel_instance)

    output_corrections_file_path = os.path.join(save_path, "corrections_" + accel_instance.label + ".madx")
    if not os.path.isfile(output_corrections_file_path):
        corrections_file_path = os.path.join(twiss_directory, "corrections_" + accel_instance.label + ".madx")
        if os.path.isfile(corrections_file_path):
            utils.iotools.copy_item(corrections_file_path, save_path)
        else:
            print("Cannot find corrections file, will create an empty file.")
            with open(os.path.join(save_path, "corrections_" + accel_instance.label + ".madx"), "a") as corrections_file:
                corrections_file.write(correction_file_comments)
    else:
        print("corrections file found in output path.")


def _get_corrections_file_comments_for_ip(accel_instance):
    comments_string = ""
    for variable in accel_instance.get_segment_vars():
        comments_string += "! " + variable + "\n"
    return comments_string


def _runmad(accel_instance, path, madx_file_path, log_file_path):
    creator.create_model(accel_instance, "segment", path,
                         logfile=log_file_path, writeto=madx_file_path)
    print("MAD done, log file:", log_file_path)


def _try_to_load_twiss(file_path):
    syserr = sys.stderr
    if not os.path.isfile(file_path):
        return None
    try:
        sys.stderr = open(os.devnull, "w")
        twiss_data = twiss(file_path)
    except ValueError:
        twiss_data = None
    finally:
        sys.stderr = syserr
    return twiss_data


def _get_twiss_for_one_of(*file_names):
    for file_name in file_names:
        single_path = _join_output_with(file_name)
        if _exists(single_path):
            return twiss(single_path)
    raise IOError("None of the files exist:\n\t" + "\n\t".join(file_names))


def _exists(path_to_file):
    return os.path.isfile(path_to_file) and os.path.exists(path_to_file)


def _all_exists_in_output_path(*file_names):
    for file_name in file_names:
        if not _exists(_join_output_with(file_name)):
            return False
    return True


__output_path = None
__save_path = None


def _set_output_path(output_path):
    global __output_path
    __output_path = output_path


def _join_output_with(a_string):
    return os.path.join(__output_path, a_string)


def _set_save_path(save_path):
    global __save_path
    __save_path = save_path


def _get_twiss_for_file_in_save_path(file_name):
    return twiss(os.path.join(__save_path, file_name))


class _PropagatedModels(object):

    def __init__(self, save_path, element_name):
        self.__save_path = save_path

        self.corrected = self.__get_twiss_for_file('twiss_' + element_name + '_cor.dat')
        self.propagation = self.__get_twiss_for_file('twiss_' + element_name + '.dat')
        self.back_propagation = self.__get_twiss_for_file('twiss_' + element_name + '_back.dat')
        self.corrected_back_propagation = self.__get_twiss_for_file('twiss_' + element_name + '_cor_back.dat')

    def __get_twiss_for_file(self, file_name):
        twiss_file = _try_to_load_twiss(os.path.join(self.__save_path, file_name))
        if twiss_file is None:
            raise SegmentBySegmentError(
                "Cannot load " + file_name + ". See mad log."
            )
        return twiss_file


class _Summaries(object):

    def __init__(self, save_path):
        self.beta = sbs_beta_writer.get_beta_summary_file(save_path)
        self.coupling = sbs_coupling_writer.get_coupling_summary_file(save_path)
        self.dispersion = sbs_dispersion_writer.get_dispersion_summary_file(save_path)
        self.chrom = sbs_chromatic_writer.get_chrom_summary_file(save_path)

    def write_summaries_to_files(self):
        if not self.beta._TfsFileWriter__tfs_table.is_empty():
            self.beta.write_to_file()
        if not self.coupling._TfsFileWriter__tfs_table.is_empty():
            self.coupling.write_to_file()
        if not self.dispersion._TfsFileWriter__tfs_table.is_empty():
            self.dispersion.write_to_file()
        if not self.chrom._TfsFileWriter__tfs_table.is_empty():
            self.chrom.write_to_file()


class _InputData(object):

    def __init__(self, measurement_path, w_path=None):
        _set_output_path(measurement_path)
        self.beta_x = _get_twiss_for_one_of("getbetax_free.out", "getbetax.out")
        self.QXX = self.beta_x.Q1
        self.QYY = self.beta_x.Q2

        self.beta_y = _get_twiss_for_one_of("getbetay_free.out", "getbetay.out")

        self.amplitude_beta_x = _get_calibrated_betas("x")
        self.amplitude_beta_y = _get_calibrated_betas("y")

        self.phase_x = _get_twiss_for_one_of("getphasex_free.out", "getphasex.out")
        self.phase_y = _get_twiss_for_one_of("getphasey_free.out", "getphasey.out")

        self.total_phase_x = _get_twiss_for_one_of("getphasetotx_free.out", "getphasetotx.out")
        self.total_phase_y = _get_twiss_for_one_of("getphasetoty_free.out", "getphasetoty.out")

        self.has_dispersion = False
        self.has_coupling = False
        self.couple_method = "_free"

        ### check if dispersion exist
        if self.__try_to_load_dispersion_files():
            self.has_dispersion = True
            print("Dispersion files OK", self.dispersion_x.DX[0], self.dispersion_x.NAME[0])
        else:
            self.has_dispersion = False
            print("No dispersion files... will continue without taking into account dispersion")

        ### check if coupling exist
        if self.__try_to_load_coupling_files():
            self.has_coupling = True
            self.couple_method = "driven"
        else:
            self.has_coupling = False

        if _exists(_join_output_with("getcouple_free.out")):
            self.couple = twiss(_join_output_with("getcouple_free.out"))
            self.couple_method = "_free"
            self.has_coupling = True
            print("Free coupling found")
        if not self.has_coupling:
            print("No coupling file... will continue without taking into account coupling")

        ### check if chromatic exists
        if w_path is None:
            w_path = measurement_path
        self.has_chromatic = self.__try_to_load_chromatic_files(w_path)

    def __try_to_load_dispersion_files(self):
        if _all_exists_in_output_path("getDx.out", "getNDx.out", "getDy.out"):
            self.dispersion_x = self.__try_to_load_twiss_from_output("getDx.out")
            self.normalized_dispersion_x = self.__try_to_load_twiss_from_output("getNDx.out")
            self.dispersion_y = self.__try_to_load_twiss_from_output("getDy.out")
            return (not self.dispersion_x is None and
                    not self.normalized_dispersion_x is None and
                    not self.dispersion_y is None)
        else:
            return False

    def __try_to_load_coupling_files(self):
        if _all_exists_in_output_path("getcouple.out", "getcoupleterms.out"):
            self.couple = self.__try_to_load_twiss_from_output("getcouple.out")
            self.couple_terms = self.__try_to_load_twiss_from_output("getcoupleterms.out")
            return (not self.couple is None and
                    not self.couple_terms is None)
        else:
            return False

    def __try_to_load_chromatic_files(self, w_path):
        wx_twiss = _try_to_load_twiss(os.path.join(w_path, "chrombetax.out"))
        wy_twiss = _try_to_load_twiss(os.path.join(w_path, "chrombetay.out"))
        if wx_twiss is None or wy_twiss is None:
            print("No chromatic files... will continue without taking into account chromatic")
            return False
        else:
            self.wx = wx_twiss
            self.wy = wy_twiss
            return True

    def __try_to_load_twiss_from_output(self, file_name):
        try:
            twiss_data = twiss(_join_output_with(file_name))
        except ValueError:
            print("Imposible to read twiss file ", file_name, file=sys.stderr)
            return None
        return twiss_data


class SegmentBySegmentError(Exception):
    pass


#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":
    (_accel_cls, _options) = _parse_args()

    return_value = main(_accel_cls, _options)
