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
from math import sqrt, cos, sin, pi, tan
import math

import numpy as np

import __init__  # @UnusedImport used for appending paths
import Utilities.iotools
from Python_Classes4MAD.metaclass import twiss
from Utilities import tfs_file_writer
import sbs_beta_writer
import sbs_phase_writer
import sbs_dispersion_writer
from SegmentBySegment import sbs_coupling_writer


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

    measurement_path = options.path
    input_data = _InputData(measurement_path)

    save_path = options.save + os.path.sep
    Utilities.iotools.create_dirs(save_path)

    elements_data = options.segf.split(',')
    error_cut = float(options.cuts)
    selected_accelerator = options.accel

    twiss_file = options.twiss
    if twiss_file == "./":
        twiss_file = os.path.join(options.bb, "MODEL", options.accel, "nominal.opt", "twiss_elements.dat")
    print "Twiss file ", twiss_file
    input_model = twiss(twiss_file)

    twiss_directory = os.path.dirname(twiss_file) + os.path.sep

    elements_names, start_bpms, end_bpms = structure_elements_info(elements_data)

    summaries = _Summaries(save_path)

    for element_name in elements_names:

        print element_name

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

        if input_data.has_coupling:
            f_coupling_parameters = get_coupling_parameters(input_data, start_bpm_name)
        else:
            f_coupling_parameters = [0, 0, 0, 0, 0, 0]
            print "No coupling"

        print "Madpass", options.madpass
        if str(options.madpass) == "0":
            print "Going to run4mad"
            run4mad(save_path,
                    start_bpm_horizontal_data,
                    start_bpm_vertical_data,
                    end_bpm_horizontal_data,
                    end_bpm_vertical_data,
                    start_bpm_dispersion,
                    end_bpm_dispersion,
                    start_bpm_name,
                    end_bpm_name,
                    element_name,
                    f_coupling_parameters,
                    options.path + "/",
                    twiss_directory,
                    input_data.couple_method)

        else:
            runmad(save_path, element_name)
            print "Just rerunning mad"

        reversetable(save_path, element_name)

        propagated_models = _PropagatedModels(save_path, element_name)

        print "Writing data for function ", element_name
        getAndWriteData(element_name,
                        input_data,
                        input_model,
                        propagated_models,
                        save_path,
                        is_element,
                        selected_accelerator,
                        summaries)

        # gnuplot  ### TODO: what is this? (jcoellod)
        if is_element == 0:
            beta4plot = start_bpm_horizontal_data[2]
            startpos = input_model.S[input_model.indx[start_bpm_name]]
            endpos = input_model.S[input_model.indx[end_bpm_name]]
            print startpos, endpos
            run4plot(save_path, startpos, endpos, beta4plot, options.bb, measurement_path, element_name, input_data.QXX, input_data.QYY, options.accel, input_data.couple_method)

    summaries.write_summaries_to_files()
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

        print elements_data[position], position
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
        print "Segment has been choosen"
        is_element = False
        segment = [start_bpms[element_name], end_bpms[element_name]]
        start_bpm_name, end_bpm_name = filterandfind(input_data.beta_x,
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
        print "Element has been choosen"
        is_element = True
        start_bpm_name, end_bpm_name = filterandfind(input_data.beta_x,
                                                     input_data.beta_y,
                                                     element_name,
                                                     [],
                                                     twiss_data,
                                                     errorcut)
    return start_bpm_name, end_bpm_name, is_element


def gather_data(input_data, startbpm, endbpm):
    start_bpm_horizontal_data = [input_data.beta_x.BETX[input_data.beta_x.indx[startbpm]],
        sqrt(input_data.beta_x.STDBETX[input_data.beta_x.indx[startbpm]] ** 2 + input_data.beta_x.ERRBETX[input_data.beta_x.indx[startbpm]] ** 2),
        input_data.beta_x.ALFX[input_data.beta_x.indx[startbpm]],
        sqrt(input_data.beta_x.ERRALFX[input_data.beta_x.indx[startbpm]] ** 2 + input_data.beta_x.STDALFX[input_data.beta_x.indx[startbpm]] ** 2)]
    start_bpm_vertical_data = [input_data.beta_y.BETY[input_data.beta_y.indx[startbpm]],
        sqrt(input_data.beta_y.STDBETY[input_data.beta_y.indx[startbpm]] ** 2 + input_data.beta_y.ERRBETY[input_data.beta_y.indx[startbpm]] ** 2),
        input_data.beta_y.ALFY[input_data.beta_y.indx[startbpm]],
        sqrt(input_data.beta_y.ERRALFY[input_data.beta_y.indx[startbpm]] ** 2 + input_data.beta_y.STDALFY[input_data.beta_y.indx[startbpm]] ** 2)]
    end_bpm_horizontal_data = [input_data.beta_x.BETX[input_data.beta_x.indx[endbpm]],
        sqrt(input_data.beta_x.STDBETX[input_data.beta_x.indx[endbpm]] ** 2 + input_data.beta_x.ERRBETX[input_data.beta_x.indx[endbpm]] ** 2),
        input_data.beta_x.ALFX[input_data.beta_x.indx[endbpm]],
        sqrt(input_data.beta_x.ERRALFX[input_data.beta_x.indx[endbpm]] ** 2 + input_data.beta_x.STDALFX[input_data.beta_x.indx[endbpm]] ** 2)]
    end_bpm_vertical_data = [input_data.beta_y.BETY[input_data.beta_y.indx[endbpm]],
        sqrt(input_data.beta_y.STDBETY[input_data.beta_y.indx[endbpm]] ** 2 + input_data.beta_y.ERRBETY[input_data.beta_y.indx[endbpm]] ** 2),
        input_data.beta_y.ALFY[input_data.beta_y.indx[endbpm]],
        sqrt(input_data.beta_y.ERRALFY[input_data.beta_y.indx[endbpm]] ** 2 + input_data.beta_y.STDALFY[input_data.beta_y.indx[endbpm]] ** 2)]
    return start_bpm_horizontal_data, start_bpm_vertical_data, end_bpm_horizontal_data, end_bpm_vertical_data


def get_dispersion_info(input_data, startbpm, endbpm):
    if startbpm in input_data.dispersion_x.indx:
        dxx = input_data.dispersion_x.DX[input_data.dispersion_x.indx[startbpm]]
        dxp = input_data.dispersion_x.DPX[input_data.dispersion_x.indx[startbpm]]
        dxe = input_data.dispersion_x.STDDX[input_data.dispersion_x.indx[startbpm]]
    else:
        dxx = 0
        dxp = 0
        dxe = 0
        print "Start BPM for horizontal dispersion not found"
    if startbpm in input_data.dispersion_y.indx:
        dyy = input_data.dispersion_y.DY[input_data.dispersion_y.indx[startbpm]]
        dyp = input_data.dispersion_y.DPY[input_data.dispersion_y.indx[startbpm]]
        dye = input_data.dispersion_y.STDDY[input_data.dispersion_y.indx[startbpm]]
    else:
        dyy = 0
        dyp = 0
        dye = 0
        print "Start BPM for vertical dispersion not found"
    start_bpm_dispersion = [dxx,
          dxp,
          dyy,
          dyp,
          dxe,
          dye]

    if endbpm in input_data.dispersion_x.indx:
        dxx = input_data.dispersion_x.DX[input_data.dispersion_x.indx[endbpm]]
        dxp = input_data.dispersion_x.DPX[input_data.dispersion_x.indx[endbpm]]
        dxe = input_data.dispersion_x.STDDX[input_data.dispersion_x.indx[endbpm]]
    else:
        dxx = 0
        dxp = 0
        dxe = 0
        print "End BPM for horizontal dispersion not found"
    if endbpm in input_data.dispersion_y.indx:
        dyy = input_data.dispersion_y.DY[input_data.dispersion_y.indx[endbpm]]
        dyp = input_data.dispersion_y.DPY[input_data.dispersion_y.indx[endbpm]]
        dye = input_data.dispersion_y.STDDY[input_data.dispersion_y.indx[endbpm]]
    else:
        dyy = 0
        dyp = 0
        dye = 0
        print "End BPM for vertical dispersion not found"
    end_bpm_dispersion = [dxx,
           dxp,
           dyy,
           dyp,
           dxe,
           dye]

    return start_bpm_dispersion, end_bpm_dispersion


def get_coupling_parameters(input_data, startbpm):
    print "Coupling"
    if startbpm in input_data.couple.indx:
        f1001r = input_data.couple.F1001R[input_data.couple.indx[startbpm]]
        f1001i = input_data.couple.F1001I[input_data.couple.indx[startbpm]]
        f1010r = input_data.couple.F1010R[input_data.couple.indx[startbpm]]
        f1010i = input_data.couple.F1010I[input_data.couple.indx[startbpm]]
        f1001std = input_data.couple.FWSTD1[input_data.couple.indx[startbpm]]
        f1010std = input_data.couple.FWSTD2[input_data.couple.indx[startbpm]]
    else:
        f1001r = 0
        f1001i = 0
        f1010r = 0
        f1010i = 0
        f1001std = 0
        f1010std = 0
        print "Start BPM ", startbpm, " not found in coupling measurement => values=0"
    f_coupling_parameters = [f1001r, f1001i, f1010r, f1010i, f1001std, f1010std]
    return f_coupling_parameters


#===================================================================================================
# helper-functions
#===================================================================================================


def modelIntersectgetf(exp, model):
    bpmsin = []
    for bpm in exp.NAME:
            try:
                    check = model.indx[bpm.upper()]
                    bpmsin.append([model.S[check], bpm])
            except:
                    print bpm, "Not in Model"
    if len(bpmsin) == 0:
            print "Zero intersection of Exp and Model"
            print "Please, provide a good Dictionary"
            print "Now we better leave!"
            sys.exit()
    bpmsin.sort()
    return bpmsin


def intersect(list_of_files):
    '''Pure intersection of all bpm names in all files '''
    if len(list_of_files) == 0:
        print "Nothing to intersect!!!!"
        sys.exit()
    z = list_of_files[0].NAME
    for b in list_of_files:
        z = filter(lambda x: x in z, b.NAME)
    #SORT by S
    result = []
    x0 = list_of_files[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))
    result.sort()
    return result


def filterandfind(beta_x_twiss, beta_y_twiss, element_name, segment_bpms_names, model, errorcut):
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

                #print (beta_x>0),(beta_y>0),(total_err_x<beta_x),(total_err_y<beta_y),(total_err_x>0),(total_err_y>0),(((total_err_x/beta_x)*100)<errorcut),(((total_err_y/beta_y)*100)<errorcut)

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

    print selected_left_bpm, selected_right_bpm, " Will be used for the propogation"

    return [selected_left_bpm, selected_right_bpm]


def TransverseDampers(twissp,twissb,element,model,savepath,phasex,phasey,errors):
    '''
    Function for getting the phase advance between the dampers and pick-ups

    :Parameters:
        'twissp': twiss?
            containing propagation from BPM to damper
        'twissb': twiss?
            containing back propagation from BPM to damper
        'element': string?
            element name
        'model': twiss
            twiss model
        'savepath': string
            where to save file
        'phasex':
            measured phase advances
        'phasey':
            measured phase advances
        'errors':
            twissmin,twissminb,twissmax,twissmaxb
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    # initial
    dampermap={}
    bpmmap={}
    path=savepath+"phases.out"
    accel=model.SEQUENCE
    if os.path.exists(path):
        writer=open(path,"a")
        firsttime=0
    else:
        writer=open(path,"w")
        firsttime=1
        print >> writer,"@ BEAM %s",accel
        print >> writer,"@ NIM %s Not In Measurement"
        print >> writer,"* NAME1 NAME2 S1 S2 PHX PHXe PHXM PHY PHYe PHYM"
        print >> writer,"$ %s %s %le %le %le %le %le %le %le %le "


    #beam1
    dampermap['ADTKH.D5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.C5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.B5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.A5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']

    dampermap['ADTKV.A5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.B5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.C5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.D5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']

    bpmmap['LHCB1_1']=['BPMC.9L4.B1','BPMC.7L4.B1','BPMC.8L4.B1'] #H
    bpmmap['LHCB1_2']=['BPMCA.7R4.B1','BPMC.9R4.B1','BPMC.8R4.B1'] #V

    #beam2
    dampermap['ADTKV.D5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.C5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.B5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.A5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']

    dampermap['ADTKH.A5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.B5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.C5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.D5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']


    bpmmap['LHCB2_1']=['BPMCA.7R4.B2','BPMC.9R4.B2','BPMC.8R4.B2'] #H
    bpmmap['LHCB2_2']=['BPMC.9L4.B2','BPMC.7L4.B2','BPMC.8L4.B2'] #V

    ## main
    if firsttime==1:

        #
        counter=2
        #

        for count in range(0,counter):

            count=count+1

            print "counter for loop ", count

            bpms=bpmmap[accel+"_"+str(count)]
            refBpm=bpms[0]
            eightbpm=bpms[2]
            endbpm=bpms[1]

            # hor bpms's
            try:
                in1x=phasex.indx[refBpm]

            except:
                in1x=-1

            if (in1x<len(phasex.indx)-1) and (in1x!=-1): # when BPM is not last
                if phasex.NAME[phasex.indx[refBpm]+1]==endbpm:
                    in2x=1
                    eightinphasex=0
                elif phasex.NAME[phasex.indx[refBpm]+2]==endbpm:
                    in2x=1
                    eightinphasex=1
                else: # not in measurement
                    in2x=-1
                    eightinphasex=0
            elif (in1x!=-1): # when BPM is last
                if phasex.NAME[0]==endbpm:
                    in2x=1
                    eightinphasex=0
                elif phasex.NAME[1]==endbpm:
                    in2x=1
                    eightinphasex=1
                else: # not in measurement
                    in2x=-1
                    eightinphasex=0
            else:
                in2x=-1
                eightinphasex=0


            # ver bpm's
            try:
                in1y=phasey.indx[refBpm]
            except:
                in1y=-1

            if (in1y<len(phasey.indx)-1) and (in1y!=-1): # when BPM is not last
                if phasey.NAME[phasey.indx[refBpm]+1]==endbpm:
                    in2y=1
                    eightinphasey=0
                elif phasey.NAME[phasey.indx[refBpm]+2]==endbpm:
                    in2y=1
                    eightinphasey=1
                else: # not in measurement
                    in2y=-1
                    eightinphasey=0
            elif in1y!=-1: # when BPM is last
                if phasey.NAME[0]==endbpm:
                    in2y=1
                    eightinphasey=0
                elif phasey.NAME[1]==endbpm:
                    in2y=1
                    eightinphasey=1
                else: # not in measurement
                    in2y=-1
                    eightinphasey=0
            else:
                in2y=-1
                eightinphasey=0





            ###### H plane
            if (in1x!=-1) and (in2x!=-1):
                if eightinphasex==1:
                    phaseh="%.5f" %float(phasex.PHASEX[phasex.indx[refBpm]]+phasex.PHASEX[phasex.indx[eightbpm]])
                    errphaseh="%.5f" %float(phasex.STDPHX[phasex.indx[refBpm]]+phasex.STDPHX[phasex.indx[eightbpm]])
                    phasemodelh="%.5f" %float(phasex.PHXMDL[phasex.indx[refBpm]]+phasex.PHXMDL[phasex.indx[eightbpm]])
                else:
                    phaseh="%.5f" %float(phasex.PHASEX[phasex.indx[refBpm]])
                    errphaseh="%.5f" %float(phasex.STDPHX[phasex.indx[refBpm]])
                    phasemodelh="%.5f" %float(phasex.PHXMDL[phasex.indx[refBpm]])
            else:
                print "Horizontal plane not found for transverse dampers pick-ups"
                phaseh='NIM'
                errphaseh='NIM'
                phasemodelh='NIM'

            ###### V plane
            print in1y,in2y,eightinphasey
            if (in1y!=-1) and (in2y!=-1):
                if eightinphasey==1:
                    phasev="%.5f" %float(phasey.PHASEY[phasey.indx[refBpm]]+phasey.PHASEY[phasey.indx[eightbpm]])
                    errphasev="%.5f" %float(phasey.STDPHY[phasey.indx[refBpm]]+phasey.STDPHY[phasey.indx[eightbpm]])
                    phasemodelv="%.5f" %float(phasey.PHYMDL[phasey.indx[refBpm]]+phasey.PHYMDL[phasey.indx[eightbpm]])
                else:

                    phasev="%.5f" %float(phasey.PHASEY[phasey.indx[refBpm]])
                    errphasev="%.5f" %float(phasey.STDPHY[phasey.indx[refBpm]])
                    phasemodelv="%.5f" %float(phasey.PHYMDL[phasey.indx[refBpm]])

            else:
                print "Vertical plane not found for transverse dampers pick-ups"
                phasev='NIM'
                errphasev='NIM'
                phasemodelv='NIM'


            print >> writer,refBpm,endbpm,model.S[model.indx[refBpm]],model.S[model.indx[endbpm]],phaseh,errphaseh,phasemodelh,phasev,errphasev,phasemodelv

    if element in dampermap:
        bpms=dampermap[element]
        passs=1
    else:
        print "WARN: Element ",element," Not found in dampermap"
        passs=0

    if passs==1:
        muxbpm1=twissp.MUX[twissp.indx[bpms[0]]]
        muxbpm2=twissp.MUX[twissp.indx[bpms[1]]]
        mux=twissp.MUX[twissp.indx[element]]
        muybpm1=twissp.MUY[twissp.indx[bpms[0]]]
        muybpm2=twissp.MUY[twissp.indx[bpms[1]]]
        muy=twissp.MUY[twissp.indx[element]]

        muxbpm1m=model.MUX[model.indx[bpms[0]]]
        muxbpm2m=model.MUX[model.indx[bpms[1]]]
        muxm=model.MUX[model.indx[element]]

        muybpm1m=model.MUY[model.indx[bpms[0]]]
        muybpm2m=model.MUY[model.indx[bpms[1]]]
        muym=model.MUY[model.indx[element]]

        pha1xp=abs(mux-muxbpm1)
        pha2xp=abs(muxbpm2-mux)
        pha1yp=abs(muy-muybpm1)
        pha2yp=abs(muybpm2-muy)

        muxbpm1=twissb.MUX[twissb.indx[bpms[0]]]
        muxbpm2=twissb.MUX[twissb.indx[bpms[1]]]
        mux=twissb.MUX[twissb.indx[element]]
        muybpm1=twissb.MUY[twissb.indx[bpms[0]]]
        muybpm2=twissb.MUY[twissb.indx[bpms[1]]]
        muy=twissb.MUY[twissb.indx[element]]


        pha1xb=abs(muxbpm1-mux)
        pha2xb=abs(mux-muxbpm2)
        pha1yb=abs(muybpm1-muy)
        pha2yb=abs(muy-muybpm2)


        pha1xm=abs(muxm-muxbpm1m)
        pha2xm=abs(muxbpm2m-muxm)
        pha1ym=abs(muym-muybpm1m)
        pha2ym=abs(muybpm2m-muym)




        #mean
        pha1x=(pha1xp+pha1xb)/2
        pha1y=(pha1yp+pha1yb)/2
        pha2x=(pha2xp+pha2xb)/2
        pha2y=(pha2yp+pha2yb)/2

        #error
        twissmin=errors[0]
        twissminb=errors[1]
        twissmax=errors[2]
        twissmaxb=errors[3]

        minn=abs(abs(twissmin.MUX[twissmin.indx[element]]-twissmin.MUX[twissmin.indx[bpms[0]]])-abs(twissminb.MUX[twissminb.indx[element]]-twissminb.MUX[twissminb.indx[bpms[0]]]))
        maxx=abs(abs(twissmax.MUX[twissmax.indx[element]]-twissmax.MUX[twissmax.indx[bpms[0]]])-abs(twissmaxb.MUX[twissmaxb.indx[element]]-twissmaxb.MUX[twissmaxb.indx[bpms[0]]]))
        pha1xe=minn+maxx


        minn=abs(abs(twissmin.MUX[twissmin.indx[bpms[1]]]-twissmin.MUX[twissmin.indx[element]])-abs(twissminb.MUX[twissminb.indx[bpms[1]]]-twissminb.MUX[twissminb.indx[element]]))
        maxx=abs(abs(twissmax.MUX[twissmax.indx[bpms[1]]]-twissmax.MUX[twissmax.indx[element]])-abs(twissmaxb.MUX[twissmaxb.indx[bpms[1]]]-twissmaxb.MUX[twissmaxb.indx[element]]))
        pha2xe=minn+maxx

        minn=abs(abs(twissmin.MUY[twissmin.indx[element]]-twissmin.MUY[twissmin.indx[bpms[0]]])-abs(twissminb.MUY[twissminb.indx[element]]-twissminb.MUY[twissminb.indx[bpms[0]]]))
        maxx=abs(abs(twissmax.MUY[twissmax.indx[element]]-twissmax.MUY[twissmax.indx[bpms[0]]])-abs(twissmaxb.MUY[twissmaxb.indx[element]]-twissmaxb.MUY[twissmaxb.indx[bpms[0]]]))
        pha1ye=minn+maxx

        minn=abs(abs(twissmin.MUY[twissmin.indx[bpms[1]]]-twissmin.MUY[twissmin.indx[element]])-abs(twissminb.MUY[twissminb.indx[bpms[1]]]-twissminb.MUY[twissminb.indx[element]]))
        maxx=abs(abs(twissmax.MUY[twissmax.indx[bpms[1]]]-twissmax.MUY[twissmax.indx[element]])-abs(twissmaxb.MUY[twissmaxb.indx[bpms[1]]]-twissmaxb.MUY[twissmaxb.indx[element]]))
        pha2ye=minn+maxx

        print >> writer,bpms[0],element,model.S[model.indx[bpms[0]]],model.S[model.indx[element]],"%.5f" %float(pha1x),"%.5f" %float(pha1xe),"%.5f" %float(pha1xm),"%.5f" %float(pha1y),"%.5f" %float(pha1ye),"%.5f" %float(pha1ym)
        print >> writer,element,bpms[1],model.S[model.indx[element]],model.S[model.indx[bpms[1]]],"%.5f" %float(pha2x),"%.5f" %float(pha2xe),"%.5f" %float(pha2xm),"%.5f" %float(pha2y),"%.5f" %float(pha2ye),"%.5f" %float(pha2ym)

    writer.close()


def getIP(betameA,basetwiss,betatwiss,alfatwiss,model,phasex,phasey,name,accel,path):
    '''
    Function calculating the optics parameters at the IP

    :Parameters:
        'betameA': list
            contaning horizontal and vertical measured amplitude functions
        'basetwiss': list
            list with propogated and back propogated twiss
        'betatwiss': list
            list of twisses for beta errors
        'alfatwiss': list
            list of twisses for alfa errors
        'model': twiss
            twiss model
        'phasex':
            measured phase advances
        'phasey':
            measured phase advances
        'name': string
            name of the IP
        'accel': string
            name of the accelerator
        'save_path': string
            where to save file
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    if os.path.isfile(path+"/IP_para.out"):
        ipfilepara=open(path+"/IP_para.out","a")
    else:
        ipfilepara=open(path+"/IP_para.out","w")
        print >> ipfilepara,"* NAME S BETX EBETX X EX BETY EBETY Y EY BETX_ph EBETX_ph BETY_ph EBETY_ph"
        print >> ipfilepara,"$ %s %le %le %le %le %le %le %le %le %le %le %le %le %le"


    if os.path.isfile(path+"/IP_pro_x.out"):
        ipfileprox=open(path+"/IP_pro_x.out","a")
    else:
        ipfileprox=open(path+"/IP_pro_x.out","w")
        print >> ipfileprox,"* NAME S BETSTARX EBETSTARX X[cm] EX[cm] BETIPX EBETIPX ALFIPX EALFIPX"
        print >> ipfileprox,"$ %s %le %le %le %le %le %le %le %le %le"

    if os.path.isfile(path+"/IP_pro_y.out"):
        ipfileproy=open(path+"/IP_pro_y.out","a")
    else:
        ipfileproy=open(path+"/IP_pro_y.out","w")
        print >> ipfileproy,"* NAME S BETY EBETY Y[cm] EY[cm] BETIPY EBETIPY ALFIPY EALFIPY"
        print >> ipfileproy,"$ %s %le %le %le %le %le %le %le %le %le"

    ## constructor
    if "B1" in accel:
        accelb="B1"
    else:
        accelb="B2"
    twissp=basetwiss[0]
    twissb=basetwiss[1]

    twissbmip=betatwiss[0]
    twissbmib=betatwiss[1]
    twissbmap=betatwiss[2]
    twissbmab=betatwiss[3]

    twissamip=alfatwiss[0]
    twissamib=alfatwiss[1]
    twissamap=alfatwiss[2]
    twissamab=alfatwiss[3]

    ampbetx=betameA[0]
    ampbety=betameA[1]

    ## beta from parabola calculation
    # Map bpm's
    bpmmap={}
    bpmmap['IP1']=["BPMSW.1L1.","BPMSW.1R1."]
    bpmmap['IP2']=["BPMSW.1L2.","BPMSW.1R2."]
    bpmmap['IP5']=["BPMSW.1L5.","BPMSW.1R5."]
    bpmmap['IP8']=["BPMSW.1L8.","BPMSW.1R8."]

    if(name in bpmmap):
        # Removed unusedphix, ephix, phiy, ephiy.
        # If needed --> git commit eca47b2d06ab7eab911c830f30dc222a8bcd91a0
        # --vimaier
        betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_p,ebetastar_h_p,betastar_v_p,ebetastar_v_p=getIPfrompara(bpmmap[name][0]+accelb,bpmmap[name][1]+accelb,ampbetx,ampbety,phasex,phasey)
    else:
        betastar_h_p=ebetastar_h_p=betastar_v_p=ebetastar_v_p=betastar_h=errbx=location_h=errsx=betastar_v=errby=location_v=errsy=0

    print >> ipfilepara,name,model.S[model.indx[name]],betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_p,ebetastar_h_p,betastar_v_p,ebetastar_v_p
    ipfilepara.close()

    ## beta from propogations
    #gathering info
    #error alfa

    dalfaxp=1/abs(twissamap.ALFX[twissamap.indx[name]]-twissamip.ALFX[twissamip.indx[name]])
    dalfaxb=1/abs(twissamab.ALFX[twissamab.indx[name]]-twissamib.ALFX[twissamib.indx[name]])

    dalfayp=1/abs(twissamap.ALFY[twissamap.indx[name]]-twissamip.ALFY[twissamip.indx[name]])
    dalfayb=1/abs(twissamab.ALFY[twissamab.indx[name]]-twissamib.ALFY[twissamib.indx[name]])

    #error beta's
    dbxp=1/abs(twissbmap.BETX[twissbmap.indx[name]]-twissbmip.BETX[twissbmip.indx[name]])
    dbxb=1/abs(twissbmab.BETX[twissbmab.indx[name]]-twissbmib.BETX[twissbmib.indx[name]])

    dbyp=1/abs(twissbmap.BETY[twissbmap.indx[name]]-twissbmip.BETY[twissbmip.indx[name]])
    dbyb=1/abs(twissbmab.BETY[twissbmab.indx[name]]-twissbmib.BETY[twissbmib.indx[name]])

    bb=twissb.ALFX[twissb.indx[name]]
    if bb<0.0:
        Balfx=abs(bb)
    else:
        Balfx=-bb

    bb=twissb.ALFY[twissb.indx[name]]
    if bb<0.0:
        Balfy=abs(bb)
    else:
        Balfy=-bb

    alfax=(1/(dalfaxp+dalfaxb))*(dalfaxp*twissp.ALFX[twissp.indx[name]]+Balfx*dalfaxb)

    ealfax=(sqrt((1/dalfaxp)**2+(1/dalfaxb)**2))/2

    alfay=(1/(dalfayp+dalfayb))*(dalfayp*twissp.ALFY[twissp.indx[name]]+Balfy*dalfayb)
    ealfay=(sqrt((1/dalfayp)**2+(1/dalfayb)**2))/2

    betxip=(1/(dbxp+dbxb))*(dbxp*twissp.BETX[twissp.indx[name]]+dbxb*twissb.BETX[twissb.indx[name]])

    betyip=(1/(dbyp+dbyb))*(dbyp*twissp.BETY[twissp.indx[name]]+dbyb*twissb.BETY[twissb.indx[name]])
    errbetxip=sqrt((1/dbxp)**2+(1/dbxb)**2)/2
    errbetyip=sqrt((1/dbyp)**2+(1/dbyb)**2)/2

    #betstar and waist
    #x
    betastar,ebetastar,waist,ewaist=getIPfromProp(betxip,errbetxip,alfax,ealfax)
    print >> ipfileprox,name,model.S[model.indx[name]],betastar,ebetastar,waist,ewaist,round(betxip,3),round(errbetxip,3),round(alfax,4),round(ealfax,4)
    ipfileprox.close()

    #y
    betastar,ebetastar,waist,ewaist=getIPfromProp(betyip,errbetyip,alfay,ealfay)
    print >> ipfileproy,name,model.S[model.indx[name]],betastar,ebetastar,waist,ewaist,round(betyip,3),round(errbetyip,3),round(alfay,4),round(ealfay,4)
    ipfileproy.close()

    ##dispersion


def getIPfromProp(betaip, errbetaip, alfaip, ealfaip):

    #values
    betastar = betaip / (1 + alfaip ** 2)
    waist = alfaip * betaip  # (sign flip)

    #errors
    ewaist=((ealfaip/abs(alfaip))+(errbetaip    /abs(betaip)))*abs(waist)
    ebetastar=sqrt(errbetaip**2+(((-2*alfaip)/(1+alfaip**2)**2)*alfaip)**2 )

    waist=waist*100 # transferring to CM!!!!
    ewaist=ewaist*100 # transferring to CM!!!!

    return round(betastar, 3), round(ebetastar, 3), round(waist, 3), round(ewaist, 3)


def getIPfrompara(bpmleft, bpmright, betax, betay, phasex, phasey):
    '''
    Function calculating beta at IP (based on Rioichy thesis)
    b(s)=b*+(s^2/b*)

    :Parameters:
        'bpmleft,bpmright': <dtype?>
            bpm's used to caculate the beta at the IP
        'ampbetx,ampbety': twiss
            twiss containing the horizontal and vertical beta from amp
        'phasex,phasey': <dtype?>
            measured phases
    :Return: <dtype?>
        beta at waist and offset
    '''
    #left horizontal
    try:
        blx=betax.BETX[betax.indx[bpmleft]]
        betax.BETXSTD[betax.indx[bpmleft]]
        hlp=1
    except:
        hlp=0

    #right horizontal
    try:
        brx=betax.BETX[betax.indx[bpmright]]
        betax.BETXSTD[betax.indx[bpmright]]
        hrp=1
    except:
        hrp=0


    #left vertical
    try:
        bly=betay.BETY[betay.indx[bpmleft]]
        betay.BETYSTD[betay.indx[bpmleft]]
        vlp=1
    except:
        vlp=0

    #right vertical
    try:
        bry=betay.BETY[betay.indx[bpmright]]
        betay.BETYSTD[betay.indx[bpmright]]
        vrp=1
    except:
        vrp=0

    #getting phaseadvance
    #inx1=phasex.indx[bpmleft]
    #inx2=phasex.indx[bpmright]
    try:
        inx1=phasex.indx[bpmleft]
        inx2=phasex.indx[bpmright]
        if (inx1==(inx2-1)) or (inx2==0 and inx1==len(phasex.NAME)-1):
            pxp=1
        else:
            pxp=0
    except:
        pxp=0

    try:
        iny1=phasey.indx[bpmleft]
        iny2=phasey.indx[bpmright]
        if iny1==(iny2-1) or (inx2==0 and inx1==len(phasex.NAME)-1):
            pyp=1
        else:
            pyp=0
    except:
        pyp=0


    #horizontal
    if(hrp==1) and (hlp==1) and (pxp==1):

        sxl=phasex.S[inx1]
        sxr=phasex.S[inx2]
        if (sxl==phasex.S[len(phasex.NAME)-1]): #if start point is in center of IP !
            L=sxr
        else:
            L=(sxr-sxl)/2

        phix=phasex.PHASEX[inx1]

        betastar_h=(2*sqrt(blx)*sqrt(brx)*sin(phix*2*pi))/(blx+brx-2*sqrt(blx)*sqrt(brx)*cos(2*pi*phix))*L
        location_h=((blx-brx)/(blx+brx-2*sqrt(blx)*sqrt(brx)*cos(2*pi*phix)))*L



        errbx=0
        errsx=0


        betastar_h_phase=(L)/tan(phix*pi)
        betastar_h_phase_e=(1+cos(2*phix*pi))*L/2


    else:
        print "WARN: Will not calculate horizontal IP"
        phix="NIM"
        betastar_h="NIM"
        location_h="NIM"
        errbx="NIM"
        errsx="NIM"
        betastar_h_phase="NIM"
        betastar_h_phase_e="NIM"
    #vertical
    if (vrp==1) and (vlp==1)and (pyp==1):

        syl=phasey.S[iny1]
        syr=phasey.S[iny2]
        if (syl==phasey.S[len(phasey.NAME)-1]): #if start point is in center of IP !
            L=syr
        else:
            L=(syr-syl)/2


        phiy=phasey.PHASEY[iny1]

        betastar_v=(2*sqrt(bly)*sqrt(bry)*sin(phiy*2*pi))/(bly+bry-2*sqrt(bly)*sqrt(bry)*cos(2*pi*phiy))*L
        location_v=((bly-bry)/(bly+bry-2*sqrt(bly)*sqrt(bry)*cos(2*pi*phiy)))*L

        betastar_v_phase=(L)/tan(phiy*pi)
        betastar_v_phase_e=(1+cos(2*phiy*pi))*L/2

        errby=0
        errsy=0
    else:
        print "WARN: Will not calculate vertical IP"
        phiy="NIM"
        betastar_v="NIM"
        location_v="NIM"
        errby="NIM"
        errsy="NIM"
        betastar_v_phase="NIM"
        betastar_v_phase_e="NIM"

    return betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_phase,betastar_h_phase_e,betastar_v_phase,betastar_v_phase_e


def getAndWriteData(element_name, input_data, input_model, propagated_models, save_path, is_element, selected_accelerator, summaries):
    '''
    Function that returns the optics function at the given element

    :Parameters:
        # TODO: rewrite this
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''

    print "INFO: Start writing files", is_element

    chromatic = []

    sbs_beta_writer.write_beta(element_name, is_element,
                               input_data.beta_x, input_data.beta_y,
                               input_model, propagated_models,
                               save_path, summaries.beta)

    sbs_phase_writer.write_phase(element_name,
                                 input_data.phase_x, input_data.phase_y, input_data.beta_x, input_data.beta_y,
                                 input_model, propagated_models, save_path)

    if input_data.has_dispersion:
        sbs_dispersion_writer.write_dispersion(element_name, is_element,
                                               input_data.dispersion_x, input_data.dispersion_y, input_data.normalized_dispersion_x,
                                               input_model, propagated_models, save_path)

    if input_data.has_coupling:
        sbs_coupling_writer.write_coupling(element_name, is_element, input_data.couple, input_model, propagated_models, save_path)


    ## to find IP
    if "IP" in element_name and is_element:
        getIP([betah[9],betav[9]],[modelp,modelb],[betah[1],betah[3],betah[2],betah[4]],[betah[5],betah[7],betah[6],betah[8]],model,phasex,phasey,namename,accel,path)
    elif ("ADT" in namename) and (switch==1):
        errors=[betah[1],betah[3],betah[2],betah[4]]
        TransverseDampers(modelp,modelb,namename,model,path,phases[0],phases[1],errors)


def weighted_average_for_SbS_elements(value1, sigma1, value2, sigma2):
    weighted_average =  (1/sigma1**2 * value1 + 1/sigma2**2 * value2) / (1/sigma1**2 + 1/sigma2**2)  # @IgnorePep8
    uncertainty_of_average = np.sqrt(1 / (1/sigma1**2 + 1/sigma2**2))  # @IgnorePep8
    weighted_rms = np.sqrt(2 * (1/sigma1**2 * (value1 - weighted_average)**2 + 1/sigma2**2 * (value2 - weighted_average)**2) / (1/sigma1**2 + 1/sigma2**2))  # @IgnorePep8
    final_error = np.sqrt(uncertainty_of_average**2 + weighted_rms**2)  # @IgnorePep8
    return weighted_average, final_error


def run4mad(save_path,
            start_bpm_horizontal_data,
            start_bpm_vertical_data,
            end_bpm_horizontal_data,
            end_bpm_vertical_data,
            start_bpm_dispersion,
            end_bpm_dispersion,
            start_bpm_name,
            end_bpm_name,
            element_name,
            f_coupling_parameters,
            exppath,
            twiss_directory,
            coupling_method):

    copy_path = options.bb
    _copy_getfterms_locally(copy_path, save_path)
    _copy_modifiers_locally(save_path, twiss_directory)

    if options.accel == "LHCB2":
        direction = -1
        start = "MKI.A5R8.B2"
        beam = "B2"

    elif options.accel == "LHCB1":
        direction = 1
        start = "MSIA.EXIT.B1"
        beam = "B1"

    wx_value, phi_x_value, wy_value, phi_y_value = _check_chromatic_functions_in_wpath(start_bpm_name)

    ### check on error propogation
    errbetx = start_bpm_horizontal_data[1]
    betx = start_bpm_horizontal_data[0]
    errbety = start_bpm_vertical_data[1]
    bety = start_bpm_vertical_data[0]
    errbetxb = end_bpm_horizontal_data[1]
    betxb = end_bpm_horizontal_data[0]
    errbetyb = end_bpm_vertical_data[1]
    betyb = end_bpm_vertical_data[0]
    f1001r = f_coupling_parameters[0]
    f1001i = f_coupling_parameters[1]
    f1010r = f_coupling_parameters[2]
    f1010i = f_coupling_parameters[3]
    f1001std = f_coupling_parameters[4]

    mad_file_name = os.path.join(save_path, 't_' + str(element_name) + '.madx')

#WARNING we are using f1001std here as before instead of f1010std which is not defined
    dict_for_replacing = dict(
            BETX=betx,
            BETY=bety,
            ERRBETX=errbetx,
            ERRBETY=errbety,
            ALFX=start_bpm_horizontal_data[2],
            ALFY=start_bpm_vertical_data[2],
            ERRALFX=start_bpm_horizontal_data[3],
            ERRALFY=start_bpm_vertical_data[3],
            DX=start_bpm_dispersion[0],
            ERRDX=start_bpm_dispersion[4],
            DY=start_bpm_dispersion[2],
            ERRDY=start_bpm_dispersion[5],
            DPX=start_bpm_dispersion[1],
            DPY=start_bpm_dispersion[3],
            ENDBX=betxb,
            ENDBY=betyb,
            ERRENDBX=errbetxb,
            ERRENDBY=errbetyb,
            ALFENDX=-end_bpm_horizontal_data[2],
            ALFENDY=-end_bpm_vertical_data[2],
            ERRALFENDX=end_bpm_horizontal_data[3],
            ERRALFENDY=end_bpm_vertical_data[3],
            DENDX=end_bpm_dispersion[0],
            ERRDENDX=end_bpm_dispersion[4],
            DENDY=end_bpm_dispersion[2],
            ERRDENDY=end_bpm_dispersion[5],
            DPENDX=-end_bpm_dispersion[1],
            DPENDY=-end_bpm_dispersion[3],
            STARTFROM=start_bpm_name.replace("-", "_"),
            ENDAT=end_bpm_name.replace("-", "_"),
            LABEL=element_name,
            ACCEL=options.accel,
            DIRE=direction,
            START=start,
            BEAM=beam,
            PATH=save_path,
            F1001R=f1001r,
            F1001I=f1001i,
            F1010R=f1010r,
            F1010I=f1010i,
            F1001maR=f1001r + f1001std,
            F1001maI=f1001i + f1001std,
            F1001miR=f1001r - f1001std,
            F1001miI=f1001i - f1001std,
            F1010maR=f1010r + f1001std,
            F1010maI=f1010i + f1001std,
            F1010miR=f1010r - f1001std,
            F1010miI=f1010i - f1001std,
            METHOD=coupling_method,
            EXP=exppath,
            WX=wx_value,
            PHIX=phi_x_value,
            WY=wy_value,
            PHIY=phi_y_value,
            WPATH=options.wpath
            )

    maskfile = os.path.join(copy_path, 'SegmentBySegment', 'job.InterpolateBetas.mask')

    Utilities.iotools.replace_keywords_in_textfile(maskfile, dict_for_replacing, mad_file_name)

    runmad(save_path, element_name)

    _prepare_watchdog_file_command(save_path, element_name, mad_file_name)


def _copy_modifiers_locally(save_path, twiss_directory):  # TODO: this has to be made platform independent
    if os.path.isfile(twiss_directory + '/modifiers.madx'):
        os.system('cp ' + twiss_directory + '/modifiers.madx' + ' ' + save_path + '/')
    else:
        os.system('touch ' + save_path + '/modifiers.madx')  # If the modifiers file does not exist create empty file


def _check_chromatic_functions_in_wpath(start_bpm_name):
    wx_value = 0
    phi_x_value = 0
    wy_value = 0
    phi_y_value = 0
    if options.wpath != "0":
        print "Chromatic save save_path,", options.wpath
        wx_twiss = _try_to_load_twiss(os.path.join(options.wpath, "wx_twiss.out"))
        wy_twiss = _try_to_load_twiss(os.path.join(options.wpath, "wy_twiss.out"))
        if not wx_twiss is None and not wy_twiss is None:
            if start_bpm_name in wx_twiss.indx:
                wx_value = wx_twiss.WX[wx_twiss.indx[start_bpm_name]]
                phi_x_value = wx_twiss.PHIX[wx_twiss.indx[start_bpm_name]]
            else:
                print "Start BPM, ", start_bpm_name, " not in WX file"
            if start_bpm_name in wy_twiss.indx:
                wy_value = wy_twiss.WY[wy_twiss.indx[start_bpm_name]]
                phi_y_value = wy_twiss.PHIY[wy_twiss.indx[start_bpm_name]]
            else:
                print "Start BPM, ", start_bpm_name, " not in WY file"
        else:
            print "No Chromatic functions (wx_twiss,wy_twiss) available at ", options.wpath
    return wx_value, phi_x_value, wy_value, phi_y_value


def _prepare_watchdog_file_command(save_path, element_name, mad_file_name):   # TODO: this has to be made platform independent
    watch_file_name = os.path.join(save_path, "watch_" + str(element_name))
    watch_file = open(watch_file_name, "w")
    print >> watch_file, "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/SegmentBySegment/watch.py " + mad_file_name + " " + save_path + "/gplot_" + str(element_name)
    watch_file.close()
    os.system("chmod +x " + watch_file_name)


def runmad(path, name):
    os.system(options.mad + ' < ' + path + 't_' + str(name) + '.madx')


def run4plot(save_path, start_point, end_point, beta4plot, beta_beat_path, measurements_path, element_name, qx, qy, accelerator, method):
    if method == "driven":
        method = ""   # patch to make it work at inj. rogelio

    dict_for_replacing = dict(PATH=save_path,
                              EndPoint=end_point,
                              StartPoint=start_point,
                              LABEL=element_name,
                              ACCEL=options.accel,
                              BETA=beta4plot,
                              QX=qx,
                              QY=qy,
                              METHOD=method,
                              MEA=measurements_path
                              )

    if (element_name == "IP8" and accelerator == "LHCB2") or (element_name == "IP2" and accelerator == "LHCB1"):
        maskfile = 'gplot.IP2IP8.mask'
    elif "RHIC" in options.accel:
        maskfile = 'gplot_RHIC.mask'
    else:
        maskfile = 'gplot.mask'
    maskfile = os.path.join(beta_beat_path, 'SegmentBySegment', maskfile)

    plotscript = os.path.join(save_path, 'gplot_' + element_name)

    # read mask file, replace all keys and write to plot script:
    Utilities.iotools.replace_keywords_in_textfile(maskfile, dict_for_replacing, plotscript)

    os.system("gnuplot " + plotscript)


#delete  TODO delete?? can this be removed??
def reversetable(path,name):
    newFile=open(path+"/twiss_"+name+"_back_rev.dat",'w')
    base=twiss(path+"/twiss_"+name+"_back.dat")

    newFile.write("* NAME                                S               BETX               ALFX               BETY               ALFY    DX       DPX     DY     DPY   MUX   MUY\n")
    newFile.write("$ %s                                %le                %le                %le                %le                %le      %le                %le      %le   %le            %le                %le \n")

    bpms=base.NAME
    s=base.S
    endpos=base.S[base.indx[bpms[len(bpms)-1]]]
    betx=base.BETX
    bety=base.BETY
    alfx=base.ALFX
    alfy=base.ALFY
    dx=base.DX
    dpx=base.DPX
    dy=base.DY
    dpy=base.DPY
    mux=base.MUX
    muy=base.MUY

    for i in range(len(bpms)):
        index=len(bpms)-1-i
        bpm=bpms[index]
        bex=betx[index]
        bey=bety[index]
        alx=alfx[index]
        aly=alfy[index]
        dex=dx[index]
        depx=dpx[index]
        dey=dy[index]
        depy=dpy[index]
        ss=s[index]
        muxx=mux[index]
        muyy=muy[index]

        #print bpm, dex

        newFile.write(str(bpm)+' '+str(endpos-ss)+' '+str(bex)+' '+str(alx)+' '+str(bey)+' '+str(aly)+'  '+str(dex)+' '+str(depx)+' '+str(dey)+' '+str(depy)+' '+str(muxx)+' '+str(muyy)+'\n')

    newFile.close()


def _try_to_load_twiss(file_path):
    try:
        twiss_data = twiss(file_path)
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
        self.back_propagation = self.__get_twiss_for_file('twiss_' + element_name + '_back_rev.dat')
        self.corrected_back_propagation = self.__get_twiss_for_file('twiss_' + element_name + '_back_rev_cor.dat')

    def __get_twiss_for_file(self, file_name):
        return twiss(os.path.join(self.__save_path, file_name))


class _Summaries(object):

    def __init__(self, save_path):
        self.beta = sbs_beta_writer.get_beta_summary_file(save_path)

        self.coupling = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_cou.out"))
        self.coupling.add_column_names("NAME", "S", "f1001", "f1001re", "f1001im", "f1010", "f1010re", "f1010im", "f1001_PLAY", "ef1001_play", "f1001re_PLAY", "f1001im_PLAY", "f1010_PLAY", "ef1010_play", "f1010re_PLAY", "f1010im_PLAY", "C11Mo", "C12Mo", "C21Mo", "C22Mo", "ANDMo", "C11_cor", "eC11_cor", "C12_cor", "eC12_cor", "C21_cor", "eC21_cor", "C22_cor", "eC22_cor", "ANG_cor", "eANG_cor", "S_MODEL")
        self.coupling.add_column_datatypes("%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le")

        self.dispersion = sbs_dispersion_writer.get_dispersion_summary_file(save_path)

        self.chrom = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_chrom.out"))
        self.chrom.add_column_names("NAME", "S", "WX_MDL", "WX_PLAY", "eWX_play", "PHIX_MDL", "PHIX_PLAY", "ePHIX_PLAY", "WY_MDL", "WY_PLAY", "eWY_play", "PHIY_MDL", "PHIY_PLAY", "ePHIY_PLAY", "S_MODEL")
        self.chrom.add_column_datatypes("%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le")

    def write_summaries_to_files(self):
        self.beta.write_to_file()
        self.coupling.write_to_file()
        self.dispersion.write_to_file()
        self.chrom.write_to_file()


class _InputData(object):

    def __init__(self, output_path):
        _set_output_path(output_path)
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
            print "No coupling file... will continue without taking into account coupling"

        if _exists(_join_output_with("getcouple_free.out")):
            self.couple = twiss(_join_output_with("getcouple_free.out"))
            self.couple_method = "_free"
            print "Free coupling found"

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
    (options, args) = parse_args()

    return_value = main(options)

    sys.exit(return_value)
