'''
Created on 11/09/09

@author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)

@version: 0.34

TODO: Description

Change history:
#  !=> SegementBySegment_0.0.py : - Construction of the program (11/09/09)
#  !=> SegementBySegment_0.1.py : - Adding the dispersion to the program (23/09/09)
#                                 - Make distuinsh between BPMs and instuments (24/09/09)
#                                 - Adding output tables (sbs...out) (21/10/09)
#                                 - Adding coupling propogation (21/10/09)
#                                 - changing output when segment is start and end (24/11/09)
#  !=> SegementBySegment_0.22.py :- Fixing some small bugs
#                                 - Added Instruments
#
#  !=> SegementBySegment_0.23.py :- Adding summary list for instruments
#                                 - Applying better cuts for beta
#  !=> SegementBySegment_0.24.py :- Including one/two beam matching (two beam for IP) (7/12/09)
#
#  !=> SegementBySegment_0.25.py :- Adding phase to the propagation (10/12/09)
#
#  !=> SegementBySegment_0.26.py :- General clean-up (19/01/10)
#                                 - Adding  Total phase to the propagation (19/01/10)
#                                 - Coupling propogation (20/01/10)
#
#  !=> SegementBySegment_0.27.py :- Fixing mistake of total phase (1/03/10)
#                                 - Include real and imaginary part in the output of coupling (1/03/10)
#
#  !=> SegementBySegment_0.28.py : -Coupling initial conditions added to MADX segment  24 March 2010
#                                  -New file getfterms.py required to convert MADX C matrix to observable f terms
#        4 April 2010              -Introducing the modifiers.madx
#                                   file, it should be in twisspath or an
#                                   empty file will be created.
#                                   The MAD mask file has also been
#                                   updated
#  !=> SegementBySegment_0.29.py : - Update for IP2  beam1 and IP8 beam 2... problem with S coordinate
#  !=>  SegementBySegment_0.30.py : - New filtering...
#  !=>  SegementBySegment_0.31.py : - New and easier filtering
#  !=>  SegementBySegment_0.31.py : - Adjusting filtering
#                                   - Adding C-matrix,angle,waist calculation,...
#                                   - Adding IP calculation from getllm... if BPM is not found you take another.
#   Added system save_path to the beginning to avoid having to define it in the shell
#
# 0.33 20120820 (tbach):
# - changed savepath = options.SAVE to savepath = options.SAVE + "/"
# - fixed errors and warnings from static code analysis (unused variables, not defined variables, variable names overriding bultin reserved key words, broken indentation)
# - fixed mixed indentation (hint: if you want to a diff which ignores the changed indentation, use sdiff -W -s
# - deleted code which is commented out
 - 0.34 vimaier 28th May 2013:
    Restructured code into parse_args(),main(),helper-functions and main invocation
    Adapted docstrings.
    Removed unused parameters from getIPfrompara
    Removed unused function getTwiss
    Previous version can be found in git commit eca47b2d06ab7eab911c830f30dc222a8bcd91a0

'''


import os
import sys
import optparse

from math import sqrt
import __init__  # @UnusedImport used for appending paths
import Utilities.iotools
import json
from Python_Classes4MAD.metaclass import twiss
from Python_Classes4MAD import madxrunner

import numpy
from numpy.linalg import inv, det
import sbs_writers.sbs_beta_writer
import sbs_writers.sbs_phase_writer
import sbs_writers.sbs_dispersion_writer
import sbs_writers.sbs_coupling_writer
import sbs_writers.sbs_chromatic_writer
import sbs_writers.sbs_special_element_writer


#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
    ''' Parses the arguments, checks for valid input and returns tupel '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                    help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                    metavar="ACCEL", default="LHCB1", dest="accel")
    parser.add_option("-f", "--path",  # assumes that output is same as input
                    help="Path to measurement files",
                    metavar="PATH", default="./", dest="path")
    parser.add_option("-i", "--path2",  # assumes that output is same as input
                    help="Path to second measurement files",
                    metavar="PATH2", default="./", dest="path2")
    parser.add_option("-s", "--start",
                    help="give start,endbpm,name (multiple allowed) eg: start1,end1,name1,start2,end2,name2,...",
                    metavar="SEGF", default="./", dest="segf")
    parser.add_option("-t", "--twiss",
                    help="basic twiss file, the modifiers.madx is assumed to be in the same direcory",
                    metavar="TWISS", default="./", dest="twiss")
    parser.add_option("-r", "--response",
                    help="switch calculate corrections 0 no, 1 yes",
                    metavar="COR", default="0", dest="cor")
    parser.add_option("-g", "--graphic",
                    help="choice between graphic or table output 0=graphic and 1=table output",
                    metavar="GRA", default="1", dest="gra")
    parser.add_option("-p", "--save",
                    help="Output path",
                    metavar="SAVE", default="./", dest="save")
    parser.add_option("-m", "--mad",  # assumes that output is same as input
                    help="mad link",
              metavar="mad", default="", dest="mad")
    parser.add_option("-b", "--bbsource",  # assumes that output is same as input
                    help="beta beat source",
                    metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="bb")
    parser.add_option("-x", "--take",  # take or create mad input, default should be 0 for creating
                    help="take or create madx 0/1",
                    metavar="mad", default="0", dest="madpass")
    parser.add_option("-c", "--cuts",
                    help="cut on error of beta in percentage",
                    metavar="cuts", default="10", dest="cuts")
    parser.add_option("-w", "--w",  # Path to Chromaticity functions
                        help="Path to  chromaticity functions, by default this is skiped",
                        metavar="wpath", default="0", dest="wpath")

    (options, args) = parser.parse_args()

    return options, args


#===================================================================================================
# main()-function
#===================================================================================================


def main(options):
    '''
    :Parameters:
        'options': Values
            Values instance with all options from OptionParser
    :Return: int
        0 if execution was successful otherwise !=0
    '''

    print "+++ Starting Segment by Segment +++"
    measurement_path = options.path
    w_path = options.wpath
    if w_path == "0":
        w_path = measurement_path
    input_data = _InputData(measurement_path, w_path)

    save_path = options.save + os.path.sep
    Utilities.iotools.create_dirs(save_path)

    elements_data = options.segf.split(',')
    error_cut = float(options.cuts)
    selected_accelerator = options.accel
    selected_accelerator_no_suffix = options.accel.replace("_II", "II")

    twiss_file = options.twiss
    if twiss_file == "./":
        twiss_file = os.path.join(options.bb, "MODEL", selected_accelerator_no_suffix, "nominal.opt", "twiss_elements.dat")

    print "Input model twiss file", twiss_file
    input_model = _try_to_load_twiss(twiss_file)
    if input_model is None:
        print >> sys.stderr, "Cannot read input model, aborting."
        sys.exit(-1)

    twiss_directory = os.path.dirname(twiss_file) + os.path.sep

    elements_names, start_bpms, end_bpms = structure_elements_info(elements_data)

    summaries = _Summaries(save_path)

    for element_name in elements_names:

        print "Started processing", element_name

        start_bpm_name, end_bpm_name, is_element = get_good_bpms(input_data, error_cut, input_model, start_bpms, end_bpms, element_name)

        (start_bpm_horizontal_data,
         start_bpm_vertical_data,
         end_bpm_horizontal_data,
         end_bpm_vertical_data) = gather_data(input_data, start_bpm_name, end_bpm_name)

        if input_data.has_dispersion:
            start_bpm_dispersion, end_bpm_dispersion = get_dispersion_info(input_data, start_bpm_name, end_bpm_name)
        else:
            start_bpm_dispersion = [0, 0, 0, 0, 0, 0]
            end_bpm_dispersion = [0, 0, 0, 0, 0, 0]
            print "No dispersion"

        element_has_coupling, f_ini, f_end = _get_coupling_parameters(input_data, start_bpm_name, end_bpm_name)

        element_has_chrom, chrom_ini, chrom_end = _get_chrom_parameters(input_data, start_bpm_name, end_bpm_name)

        if str(options.madpass) == "0":
            _run4mad(save_path,
                    start_bpm_horizontal_data,
                    start_bpm_vertical_data,
                    end_bpm_horizontal_data,
                    end_bpm_vertical_data,
                    start_bpm_dispersion,
                    end_bpm_dispersion,
                    start_bpm_name,
                    end_bpm_name,
                    element_name,
                    f_ini,
                    f_end,
                    chrom_ini,
                    chrom_end,
                    options.path + "/",
                    twiss_directory,
                    input_data.couple_method,
                    selected_accelerator,
                    options.bb,
                    options.mad)

        else:
            print "Just rerunning mad"
            mad_file_path, log_file_path = _get_files_for_mad(save_path, element_name)
            _runmad(mad_file_path, log_file_path, options.mad)

        reversetable(save_path, element_name)

        propagated_models = _PropagatedModels(save_path, element_name)

        getAndWriteData(element_name,
                        input_data,
                        input_model,
                        propagated_models,
                        save_path,
                        is_element,
                        element_has_coupling,
                        element_has_chrom,
                        selected_accelerator_no_suffix,
                        summaries)

        # TODO: This has to be fixed
        if not is_element:
            beta4plot = start_bpm_horizontal_data[0]
            startpos = input_model.S[input_model.indx[start_bpm_name]]
            endpos = input_model.S[input_model.indx[end_bpm_name]]
            run4plot(save_path, startpos, endpos, beta4plot, options.bb, measurement_path, element_name, input_data.QXX, input_data.QYY, options.accel, input_data.couple_method)
        print "Everything done for", element_name, "\n"

    summaries.write_summaries_to_files()

    print "+++  Ended Segment by Segment   +++"
    return 0

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
        print >> sys.stderr, "Something strange ....Did you properly define the input ?"
        print >> sys.stderr, "Like: BPM1,BPM2,ARC12,IP1"
        print >> sys.stderr, "Structure must be for Segment => BPML,BPMR,NAME"
        print >> sys.stderr, "For Instrument just name"
        sys.exit()
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


def get_dispersion_info(input_data, startbpm, endbpm):
    if startbpm in input_data.dispersion_x.indx:
        dxx = input_data.dispersion_x.DX[input_data.dispersion_x.indx[startbpm]]
        dxp = input_data.dispersion_x.DPX[input_data.dispersion_x.indx[startbpm]]
    else:
        dxx = 0
        dxp = 0
        print "Start BPM for horizontal dispersion not found"
    if startbpm in input_data.dispersion_y.indx:
        dyy = input_data.dispersion_y.DY[input_data.dispersion_y.indx[startbpm]]
        dyp = input_data.dispersion_y.DPY[input_data.dispersion_y.indx[startbpm]]
    else:
        dyy = 0
        dyp = 0
        print "Start BPM for vertical dispersion not found"
    start_bpm_dispersion = [dxx,
          dxp,
          dyy,
          dyp]

    if endbpm in input_data.dispersion_x.indx:
        dxx = input_data.dispersion_x.DX[input_data.dispersion_x.indx[endbpm]]
        dxp = input_data.dispersion_x.DPX[input_data.dispersion_x.indx[endbpm]]
    else:
        dxx = 0
        dxp = 0
        print "End BPM for horizontal dispersion not found"
    if endbpm in input_data.dispersion_y.indx:
        dyy = input_data.dispersion_y.DY[input_data.dispersion_y.indx[endbpm]]
        dyp = input_data.dispersion_y.DPY[input_data.dispersion_y.indx[endbpm]]
    else:
        dyy = 0
        dyp = 0
        print "End BPM for vertical dispersion not found"
    end_bpm_dispersion = [dxx,
           dxp,
           dyy,
           dyp]

    return start_bpm_dispersion, end_bpm_dispersion


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
        if not startbpm in input_data.couple.indx:
            print "Start BPM ", startbpm, " not found in coupling measurement, will not compute coupling."
        elif not endbpm in input_data.couple.indx:
            print "End BPM ", endbpm, " not found in coupling measurement, will not compute coupling."
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
            print "Start and end BPMs found in coupling measurement."
    else:
        print "No coupling measurement"
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
            print "Start BPM ", startbpm, " not found in chromatic measurement, will not compute chromatic."
        elif not endbpm in input_data.wy.indx:
            print "End BPM ", endbpm, " not found in chromatic measurement, will not compute chromatic."
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
            print "Start and end BPMs found in chromatic measurement."
    return element_has_chrom, chrom_ini, chrom_end

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
            print >> sys.stderr, element_name, " Not found in model=> System exit"
            sys.exit()
        locations_list.append(element_s)
        translate[element_s] = [True, element_name]
        print "You selected an element"
    else:
        is_segment = True
        left_bpm_name = segment_bpms_names[0]
        right_bpm_name = segment_bpms_names[1]
        if left_bpm_name in model.indx and right_bpm_name in model.indx:
            left_bpm_s = model.S[model.indx[left_bpm_name]]
            right_bpm_s = model.S[model.indx[right_bpm_name]]
        else:
            print >> sys.stderr, left_bpm_name, right_bpm_name, " Not found in model=> System exit"
            sys.exit()
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

        print "You selected a segment"

    elements_names_in_model = model.NAME

    number_of_good_bpms = 0

    #filtering
    for current_element_name in elements_names_in_model:

        if "BPM" in current_element_name:
            current_element_s = model.S[model.indx[current_element_name]]

            if current_element_name in beta_y_twiss.indx and current_element_name in beta_x_twiss.indx:
                beta_y = beta_y_twiss.BETY[beta_y_twiss.indx[current_element_name]]
                beta_x = beta_x_twiss.BETX[beta_x_twiss.indx[current_element_name]]
                err_beta_x = beta_x_twiss.ERRBETX[beta_x_twiss.indx[current_element_name]]
                stdbetax = beta_x_twiss.STDBETX[beta_x_twiss.indx[current_element_name]]
                err_beta_y = beta_y_twiss.ERRBETY[beta_y_twiss.indx[current_element_name]]
                stdbetay = beta_y_twiss.STDBETY[beta_y_twiss.indx[current_element_name]]

                total_err_x = sqrt(err_beta_x ** 2 + stdbetax ** 2)
                total_err_y = sqrt(err_beta_y ** 2 + stdbetay ** 2)

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
        print >> sys.stderr, "Not enough BPMs! System exit"
        sys.exit()

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
        print left_bpm_is_good, right_bpm_is_good

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

    print selected_left_bpm, selected_right_bpm, "Will be used for the propagation"

    return [selected_left_bpm, selected_right_bpm]


def getAndWriteData(element_name, input_data, input_model, propagated_models, save_path, is_element, element_has_coupling, element_has_chrom, selected_accelerator, summaries):
    '''
    Function that returns the optics function at the given element

    :Parameters:
        # TODO: rewrite this
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''

    print "Start writing files for", element_name

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
     beta_y, err_beta_y, alfa_y, err_alfa_y) = sbs_writers.sbs_beta_writer.write_beta(element_name, is_element,
                                                                                      input_data.beta_x, input_data.beta_y,
                                                                                      input_model,
                                                                                      propagated_models,
                                                                                      save_path, beta_summary)
    if not is_element:
        sbs_writers.sbs_phase_writer.write_phase(element_name,
                                                 input_data.total_phase_x, input_data.total_phase_y, input_data.beta_x, input_data.beta_y,
                                                 propagated_models, save_path)
    if input_data.has_dispersion:
        sbs_writers.sbs_dispersion_writer.write_dispersion(element_name, is_element,
                                                           input_data.dispersion_x, input_data.dispersion_y, input_data.normalized_dispersion_x,
                                                           input_model,
                                                           propagated_models, save_path, disp_summary)
    if element_has_coupling:
        sbs_writers.sbs_coupling_writer.write_coupling(element_name, is_element, input_data.couple, input_model, propagated_models, save_path, coupling_summary)

    if element_has_chrom:
        sbs_writers.sbs_chromatic_writer.write_chromatic(element_name, is_element, input_data.wx, input_data.wy, input_model, propagated_models, save_path, chrom_summary)

    if "IP" in element_name and is_element:
        sbs_writers.sbs_special_element_writer.write_ip(input_data.beta_x, input_data.beta_y,
                                                        beta_x, err_beta_x, alfa_x, err_alfa_x,
                                                        beta_y, err_beta_y, alfa_y, err_alfa_y,
                                                        input_model, input_data.phase_x, input_data.phase_y, element_name,
                                                        selected_accelerator, save_path)
    # TODO: This need to be fixed before using (find BPMs)
    elif False and "ADT" in element_name and is_element:  # if False to avoid going inside
        sbs_writers.sbs_special_element_writer.write_transverse_damper(propagated_models, element_name, input_model,
                                                                       save_path, input_data.phase_x, input_data.phase_y,
                                                                       input_data.beta_x, input_data.beta_y,
                                                                       selected_accelerator)


def _run4mad(save_path,
             start_bpm_horizontal_data,
             start_bpm_vertical_data,
             end_bpm_horizontal_data,
             end_bpm_vertical_data,
             start_bpm_dispersion,
             end_bpm_dispersion,
             start_bpm_name,
             end_bpm_name,
             element_name,
             f_ini,
             f_end,
             chrom_ini,
             chrom_end,
             exppath,
             twiss_directory,
             coupling_method,
             accelerator,
             bb_path,
             madx_exe_path):

    copy_path = bb_path

    is_lhc_run_ii = False
    if accelerator.endswith("_II"):
        is_lhc_run_ii = True
    accelerator = accelerator.replace("_II", "")
    _copy_modifiers_and_corrections_locally(save_path, twiss_directory, element_name, accelerator)

    if accelerator == "LHCB2":
        start = "MKI.A5R8.B2"

    elif accelerator == "LHCB1":
        start = "MSIA.EXIT.B1"

    betx_ini = start_bpm_horizontal_data[0]
    bety_ini = start_bpm_vertical_data[0]
    betx_end = end_bpm_horizontal_data[0]
    bety_end = end_bpm_vertical_data[0]

    if betx_ini < 0. or bety_ini < 0.:
        print >> sys.stderr, "Negative betas in initial BPM of segment!"
        print >> sys.stderr, "Aborting..."
        sys.exit(1)
    if betx_end < 0. or bety_end < 0.:
        print >> sys.stderr, "Negative betas in last BPM of segment!"
        print >> sys.stderr, "Aborting..."
        sys.exit(1)

    alfx_ini = start_bpm_horizontal_data[1]
    alfy_ini = start_bpm_vertical_data[1]
    alfx_end = -end_bpm_horizontal_data[1]
    alfy_end = -end_bpm_vertical_data[1]

    ini_r11, ini_r12, ini_r21, ini_r22 = _get_R_terms(betx_ini, bety_ini, alfx_ini, alfy_ini, f_ini["f1001r"], f_ini["f1001i"], f_ini["f1010r"], f_ini["f1010i"])
    end_r11, end_r12, end_r21, end_r22 = _get_R_terms(betx_end, bety_end, alfx_end, alfy_end, f_end["f1001r"], f_end["f1001i"], f_end["f1010r"], f_end["f1010i"])

    mad_file_path, log_file_path = _get_files_for_mad(save_path, element_name)

    dict_for_replacing = dict(
            BBPATH=bb_path,
            SBSPATH=os.path.join(bb_path, "SegmentBySegment"),
            LHC_RUN=2 if is_lhc_run_ii else 1,
            STARTFROM=start_bpm_name.replace("-", "_"),
            ENDAT=end_bpm_name.replace("-", "_"),
            LABEL=element_name,
            ACCEL=accelerator,
            START=start,
            PATH=save_path,
            METHOD=coupling_method,
            EXP=exppath
            )

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

    maskfile = os.path.join(copy_path, 'SegmentBySegment', 'job.InterpolateBetas.mask')
    Utilities.iotools.replace_keywords_in_textfile(maskfile, dict_for_replacing, mad_file_path)

    with open(os.path.join(save_path, "measurement_" + element_name + ".madx"), "w") as measurement_file:
        for name, value in measurement_dict.iteritems():
            print >> measurement_file, name, "=", value, ";"

    _runmad(mad_file_path, log_file_path, madx_exe_path)

    _prepare_watchdog_file_command(save_path, element_name)


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

    gamma2 = 1. / (1. + 4. * (f1001r ** 2. + f1001i ** 2. - f1010r ** 2. - f1010i ** 2.))
    c11 = (f1001i + f1010i)
    c22 = (f1001i - f1010i)
    c12 = -(f1010r - f1001r)
    c21 = -(f1010r + f1001r)
    C = numpy.reshape(2 * gamma2 * numpy.array([c11, c12, c21, c22]), (2, 2))

    R = numpy.dot(J, numpy.dot(inv(Ga), numpy.dot(C, numpy.dot(Gb, -J))))
    R = numpy.sqrt(1 + det(R)) * R
    return numpy.ravel(numpy.transpose(R))


def _get_files_for_mad(save_path, element_name):
    mad_file_name = 't_' + str(element_name) + '.madx'
    log_file_name = element_name + "_mad.log"
    mad_file_path = os.path.join(save_path, mad_file_name)
    log_file_path = os.path.join(save_path, log_file_name)
    return mad_file_path, log_file_path


def _copy_modifiers_and_corrections_locally(save_path, twiss_directory, element_name, accelerator):
    modifiers_file_path = os.path.join(twiss_directory, 'modifiers.madx')
    if os.path.isfile(modifiers_file_path):
        Utilities.iotools.copy_item(modifiers_file_path, save_path)
    else:
        print "Cannot find modifiers.madx file, will create an empty file."
        open(os.path.join(save_path, 'modifiers.madx'), "a").close()

    correction_file_comments = ""
    if element_name.lower().startswith("ip") and accelerator.upper().startswith("LHCB"):
        correction_file_comments = _get_corrections_file_comments_for_ip(element_name, accelerator)

    output_corrections_file_path = os.path.join(save_path, "corrections_" + element_name + ".madx")
    if not os.path.isfile(output_corrections_file_path):
        corrections_file_path = os.path.join(twiss_directory, "corrections_" + element_name + ".madx")
        if os.path.isfile(corrections_file_path):
            Utilities.iotools.copy_item(corrections_file_path, save_path)
        else:
            print "Cannot find corrections file, will create an empty file."
            with open(os.path.join(save_path, "corrections_" + element_name + ".madx"), "a") as corrections_file:
                corrections_file.write(correction_file_comments)
    else:
        print "corrections file found in output path."


def _get_corrections_file_comments_for_ip(element_name, accelerator):
    this_file_path = os.path.abspath(os.path.dirname(__file__))
    try:
        ip = element_name.lower().replace("ip", "")
        beam = int(accelerator.lower().replace("lhcb", ""))
    except:
        return ""
    if beam == 1:
        all_list = json.load(open(os.path.join(this_file_path, "..", "MODEL", "LHCB", "fullresponse", "LHCB1", "AllLists.json")))
        all_list_couple = json.load(open(os.path.join(this_file_path, "..", "MODEL", "LHCB", "fullresponse", "LHCB1", "AllLists_couple.json")))
    elif beam == 2 or beam == 4:
        all_list = json.load(open(os.path.join(this_file_path, "..", "MODEL", "LHCB", "fullresponse", "LHCB2", "AllLists.json")))
        all_list_couple = json.load(open(os.path.join(this_file_path, "..", "MODEL", "LHCB", "fullresponse", "LHCB1", "AllLists_couple.json")))
    else:
        return ""
    comments_string = ""
    if ip in all_list["getListsByIR"][0]:
        for variable in all_list["getListsByIR"][0][ip]:
            comments_string += "! " + variable + "\n"
    if ip in all_list["getListsByIR"][1]:
        for variable in all_list["getListsByIR"][1][ip]:
            comments_string += "! " + variable + "\n"
    for variable in all_list_couple["Qs"]:
        if "r" + str(ip) in variable or "l" + str(ip) in variable:
            comments_string += "! " + variable + "\n"
    return comments_string


def _prepare_watchdog_file_command(save_path, element_name):
    corrections_file_name = "corrections_" + element_name + ".madx"
    sbs_command = "\"" + " ".join(sys.argv) + "\""  # Gets the full command to run in the watchfile
    watch_file_name = os.path.join(save_path, "watch_" + str(element_name))
    watch_file = open(watch_file_name, "w")
    print >> watch_file, "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/SegmentBySegment/watch.py " +\
                          corrections_file_name + " " +\
                          save_path + "/gplot_" + str(element_name) + " " +\
                          sbs_command
    watch_file.close()
    os.chmod(watch_file_name, 0777)


def _runmad(file_path, log_file_path, mad_exe_path):
    return_code = madxrunner.runForInputFile(file_path, mad_exe_path, open(log_file_path, "w"))
    if return_code != 0:
        print >> sys.stderr, "MAD execution failed, see log:", log_file_path
        print >> sys.stderr, "Aborting..."
        sys.exit(return_code)
    print "MAD done, log file:", log_file_path


def run4plot(save_path, start_point, end_point, beta4plot, beta_beat_path, measurements_path, element_name, qx, qy, accelerator, method):
    if method == "driven":
        method = ""   # patch to make it work at inj. rogelio

    dict_for_replacing = dict(PATH=save_path,
                              EndPoint=end_point,
                              StartPoint=start_point,
                              LABEL=element_name,
                              ACCEL=accelerator,
                              BETA=beta4plot,
                              QX=qx,
                              QY=qy,
                              METHOD=method,
                              MEA=measurements_path
                              )

    if "RHIC" in accelerator:
        maskfile = 'gplot_RHIC.mask'
    else:
        maskfile = 'gplot.mask'
    maskfile = os.path.join(beta_beat_path, 'SegmentBySegment', maskfile)

    plotscript = os.path.join(save_path, 'gplot_' + element_name)

    # read mask file, replace all keys and write to plot script:
    Utilities.iotools.replace_keywords_in_textfile(maskfile, dict_for_replacing, plotscript)

    os.system("gnuplot " + plotscript)


#delete  TODO delete?? can this be removed??
def reversetable(save_path, element_name):
    reverse_file = open(os.path.join(save_path, "twiss_" + element_name + "_back_rev.dat"), 'w')
    original_file = _try_to_load_twiss(os.path.join(save_path, "twiss_" + element_name + "_back.dat"))

    reverse_file.write("* NAME                                S               BETX               ALFX               BETY               ALFY    DX       DPX     DY     DPY   MUX   MUY\n")
    reverse_file.write("$ %s                                %le                %le                %le                %le                %le      %le                %le      %le   %le            %le                %le \n")

    bpms = original_file.NAME
    s = original_file.S
    endpos = original_file.S[original_file.indx[bpms[len(bpms) - 1]]]
    betx = original_file.BETX
    bety = original_file.BETY
    alfx = original_file.ALFX
    alfy = original_file.ALFY
    dx = original_file.DX
    dpx = original_file.DPX
    dy = original_file.DY
    dpy = original_file.DPY
    mux = original_file.MUX
    muy = original_file.MUY

    for i in range(len(bpms)):
        index = len(bpms) - 1 - i
        bpm = bpms[index]
        bex = betx[index]
        bey = bety[index]
        alx = alfx[index]
        aly = alfy[index]
        dex = dx[index]
        depx = dpx[index]
        dey = dy[index]
        depy = dpy[index]
        ss = s[index]
        muxx = mux[index]
        muyy = muy[index]

        reverse_file.write(str(bpm) + ' ' + str(endpos - ss) + ' ' + str(bex) + ' ' + str(alx) + ' ' + str(bey) + ' ' + str(aly) + '  ' + str(dex) + ' ' + str(depx) + ' ' + str(dey) + ' ' + str(depy) + ' ' + str(muxx) + ' ' + str(muyy) + '\n')

    reverse_file.close()


def _try_to_load_twiss(file_path):
    if not os.path.isfile(file_path):
        return None
    try:
        syserr = sys.stderr
        sys.stderr = open(os.devnull, "w")
        twiss_data = twiss(file_path)
        sys.stderr = syserr
    except ValueError:
        twiss_data = None
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
            print >> sys.stderr, "Cannot load", file_name, ". See mad log. Aborting"
            sys.exit(-1)
        return twiss_file


class _Summaries(object):

    def __init__(self, save_path):
        self.beta = sbs_writers.sbs_beta_writer.get_beta_summary_file(save_path)
        self.coupling = sbs_writers.sbs_coupling_writer.get_coupling_summary_file(save_path)
        self.dispersion = sbs_writers.sbs_dispersion_writer.get_dispersion_summary_file(save_path)
        self.chrom = sbs_writers.sbs_chromatic_writer.get_chrom_summary_file(save_path)

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

        self.amplitude_beta_x = twiss(_join_output_with("getampbetax.out"))
        self.amplitude_beta_y = twiss(_join_output_with("getampbetay.out"))

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
            print "Dispersion files OK", self.dispersion_x.DX[0], self.dispersion_x.NAME[0]
        else:
            self.has_dispersion = False
            print "No dispersion files... will continue without taking into account dispersion"

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
            print "Free coupling found"
        if not self.has_coupling:
            print "No coupling file... will continue without taking into account coupling"

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
            print "No chromatic files... will continue without taking into account chromatic"
            return False
        else:
            self.wx = wx_twiss
            self.wy = wy_twiss
            return True

    def __try_to_load_twiss_from_output(self, file_name):
        try:
            twiss_data = twiss(_join_output_with(file_name))
        except ValueError:
            print >> sys.stderr, "Imposible to read twiss file ", file_name
            return None
        return twiss_data


#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":
    (_options, _args) = parse_args()

    return_value = main(_options)

    sys.exit(return_value)
