r"""
.. module: Correction.correct_coupleDy

Created on ??

Single beam correction of coupling resonances and vertical dispersion.
TODO: get better description (vimaier)

Usage example::

    python correct_coupleDy.py --accel=LHCB1
                                --path=C:\MyTemp\correct_coupleDy_out\
                                --cut=0.01
                                --errorcut=0.02,0.02
                                --modelcut=0.0,0.01
                                --rpath=H:\work/v/vimaier/public/Beta-Beat.src/
                                --MinStr=0.000001
                                --Dy=1,1,0,0,1
                                --opt=C:\MyTemp\correct_coupleDy_out\
                                --Variables=coupling_knobs

Hint: MinStr is not used in the script except of a print. The reason is possibly the author wanted to have the same set
of arguments for all correction scripts due to GUI compatibility(vimaier).

Options::

    -a ACCEL, --accelerator=ACCEL
                          What accelerator: LHCB1 LHCB2 SPS RHIC
    -p PATH, --path=PATH  Path to experimental files
    -c CUT, --cut=CUT     Singular value cut for the generalised inverse
    -e ERRORCUT, --errorcut=ERRORCUT
                          Maximum error allowed for the coupling measurement and
                          Dy measurement
    -m MODELCUT, --modelcut=MODELCUT
                          Maximum difference allowed between model and measured
                          phase
    -r RPATH, --rpath=RPATH
                          Path to BetaBeat repository (default is the afs
                          repository)
    -s MinStr, --MinStr=MinStr
                          Minimum strength of correctors in SVD correction
    -d Dy, --Dy=Dy        weight on corrections (f1001.re, f1001.im, f1010.re,
                          f1010.im, Dy)
    -o OPT, --opt=OPT     To specify the optics
    -v var, --Variables=var
                          variables split with , - if coupling_knobs is given,
                          the modelcut will be calculated automatically

.. moduleauthor:: Unknown
"""

import sys
import pickle
import os
import optparse
import re
import json

import numpy

import __init__ # @UnusedImport init will include paths
import Utilities.iotools
import Utilities.tfs_file_writer
import Utilities.math
import Python_Classes4MAD.GenMatrix_coupleDy as GenMatrix_coupleDy
import Python_Classes4MAD.metaclass as metaclass


# internal options
PRINT_DEBUG = False or sys.flags.debug  # If True, internal debug information will be printed (tbach)


#===================================================================================================
# _parse_args()-function
#===================================================================================================
def _parse_args():
    ''' Parses the arguments, checks for valid input and returns tupel '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accelerator",
                      help="What accelerator: LHCB1 LHCB2 SPS RHIC",
                      metavar="ACCEL", default="LHCB1", dest="ACCEL")
    parser.add_option("-p", "--path",
                      help="Path to experimental files",
                      metavar="PATH", default="./", dest="path")
    parser.add_option("-c", "--cut",
                      help="Singular value cut for the generalized inverse",
                      metavar="CUT", default=0.1 , dest="cut")
    parser.add_option("-e", "--errorcut",
                      help="Maximum error allowed for the coupling maesurement and Dy measurement",
                      metavar="ERRORCUT", default="0.1,0.1" , dest="errorcut")
    parser.add_option("-m", "--modelcut",
                      help="Maximum difference allowed between model and measured phase",
                      metavar="MODELCUT", default="0.1,0.1" , dest="modelcut")
    parser.add_option("-r", "--rpath",
                      help="Path to BetaBeat repository (default is the afs repository)",
                      metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
    parser.add_option("-s", "--MinStr",
                      help="Minimum strength of correctors in SVD correction",
                      metavar="MinStr", default=0.000001 , dest="MinStr")
    parser.add_option("-d", "--Dy",
                      help="weight on corrections (f1001.re, f1001.im, f1010.re, f1010.im, Dy)",
                      metavar="Dy", default="1,1,0,0,0", dest="Dy")
    parser.add_option("-o", "--opt",
                      help="To specify the optics",
                      metavar="OPT", default="nominal.opt", dest="opt")
    parser.add_option("-v", "--Variables",
                      help="variables split with , - if coupling_knobs is given, the modelcut will be calculated automatically",
                      metavar="var", default="MQSb1" , dest="var")

    (options, args) = parser.parse_args() # @UnusedVariable no args
    return options



#=======================================================================================================================
# main()-function
#=======================================================================================================================
def main(
         output_path,
         accel="LHCB1",
         singular_value_cut=0.1,
         errorcut="0.1,0.1",
         modelcut="0.1,0.1",
         beta_beat_root=Utilities.iotools.get_absolute_path_to_betabeat_root(),
         min_strength=0.000001,
         weights_on_corrections="1,1,0,0,0",
         path_to_optics_files_dir="nominal.opt",
         variables="MQSb1"
         ):

    _InputData.static_init(output_path, accel, singular_value_cut, errorcut, modelcut, beta_beat_root, min_strength, weights_on_corrections, path_to_optics_files_dir, variables)

    print "Start Correcting couple Dy"

    _generate_changeparameters_couple_file()

    print "handling data"
    if accel == "SPS":
        _handle_data_for_sps()

    if "LHC" in accel:
        _handle_data_for_lhc()

    print "Correcting couple Dy finished with weight " + str(_InputData.weights_list)


#=======================================================================================================================
# helper functions
#=======================================================================================================================

def _generate_changeparameters_couple_file():
    print "Starting loading Full Response optics"
    full_response = pickle.load(open(os.path.join(_InputData.path_to_optics_files_dir,"FullResponse_couple"), 'rb'))
    print "Loading ended"

    disp_y = _get_dispersion_y_if_set_and_available()

    path_all_lists_json_file = os.path.join(_InputData.accel_path, "AllLists_couple.json")
    knobsdict = json.load(file(path_all_lists_json_file, 'r'))
    print "Loaded json file: " + path_all_lists_json_file
    varslist = []
    for var in _InputData.variables_list:
        varslist = varslist + knobsdict[var]

    couple_twiss = _get_twiss_instance_of_getcouple()
    _InputData.model_cut_c = _calculate_automatic_model_cut_if_desired(couple_twiss)
    mad_twiss = full_response['0']
    mad_twiss.Cmatrix()
    mode = 'C'
    couple_list = GenMatrix_coupleDy.make_list(couple_twiss, mad_twiss, _InputData.model_cut_c, _InputData.error_cut_c, mode)
    mode = 'D'
    disp_y_list = GenMatrix_coupleDy.make_list(disp_y, mad_twiss, _InputData.model_cut_d, _InputData.error_cut_d, mode)

    if 0 == len(couple_list) and 0 == len(disp_y_list):
        print >> sys.stderr, "Couple and dispersion lists are empty. Are model and/or error cuts too strict?"
        sys.exit(1)

    print "entering couple input", len(couple_list)
    couple_inp = GenMatrix_coupleDy.CoupleInput(varslist, couple_list, disp_y_list, _InputData.weights_list)
    print "computing the sensitivity matrix"
    sensitivity_matrix = couple_inp.computeSensitivityMatrix(full_response)  # @UnusedVariable sensivity_matrix will be saved in couple_inp(vimaier)

    print "computing correct coupling "
    [deltas, varslist ] = GenMatrix_coupleDy.correctcouple(couple_twiss, disp_y, couple_inp, cut=_InputData.singular_value_cut, app=0, path=_InputData.output_path)
    print "deltas:", deltas


def _get_twiss_instance_of_getcouple():
    """ Loads either getcouple_free.out or getcouple.out or raises an Exception """
    file_path_free = os.path.join(_InputData.output_path, "getcouple_free.out")
    file_path = os.path.join(_InputData.output_path, "getcouple.out")
    if os.path.isfile(file_path_free):
        couple = metaclass.twiss(file_path_free)
        print "Will use free coupling"
    elif os.path.isfile(file_path):
        couple = metaclass.twiss(file_path)
        print "WARN: Free coupling not found!"
    else:
        raise ValueError("path", _InputData.output_path, " does not contain a getcouple_free.out or getcouple.out file")
    return couple


def _calculate_automatic_model_cut_if_desired(couple_twiss):
    if not"coupling_knobs" in _InputData.variables_list: # then it is not desired(vimaier)
        return _InputData.model_cut_c
    if 0 == couple_twiss.F1001W.size:
        return _InputData.model_cut_c
    print "coupling_knobs mode. Trying to do automatic correcting"
    # we want to have a value which indicates the worst 5 percent, so we sort and get the value from the 95% index (tbach)
    sortedF1001W = numpy.sort(couple_twiss.F1001W)
    five_percent_index = int(numpy.floor(sortedF1001W.size * 0.95)) # no requirement to be accurate (tbach)
    five_percent_value = sortedF1001W[five_percent_index]

    if PRINT_DEBUG:
        print "five_percent_index:", five_percent_index

    if numpy.allclose(five_percent_value, 0):
        print >> sys.stderr, "calculated value for modelcut is 0, exit"
        sys.exit(1)

    print "modelcutC will be changed from:", _InputData.model_cut_c, "to:", five_percent_value
    return five_percent_value


def _get_dispersion_y_if_set_and_available():
    dispersion_y_filenpath = os.path.join(_InputData.output_path, "getDy.out")
    if _InputData.is_dy_switch_set():
        if os.path.exists(dispersion_y_filenpath):
            return metaclass.twiss(dispersion_y_filenpath)
    return []


def _handle_data_for_sps():
    v = metaclass.twiss(os.path.join(_InputData.output_path, "changeparameters_couple.tfs"))  # @UnusedVariable
        # v will be used in VBumps.py. Don't do such things again!(vimaier)
    print '\nFor SPS, vertical bump file is imported:'
    print os.path.join(_InputData.accel_path, "Coupling/VBumps.py")
    vcorrs = None # Will be assigned in exec VBumps.py (vimaier)
    vcorrsYASP = None # Will be assigned in exec VBumpsYASP.py (vimaier)
    execfile(os.path.join(_InputData.accel_path, "Coupling", "VBumps.py")) # LOADS vcorrs
    execfile(os.path.join(_InputData.accel_path, "Coupling", "VBumpsYASP.py")) # LOADS vcorrsYASP
    #Output for YASP...
    f = open(os.path.join(_InputData.output_path, "Coupling", "changeparameters_couple.yasp"), "w")
    f.write("#PLANE V\n")
    f.write("#UNIT RAD\n")
    for vcorr in vcorrsYASP:
        print >> f, "#SETTING", vcorr,  vcorrsYASP[vcorr]
    f.close()

    #Output for Knob...
    g = open(os.path.join(_InputData.output_path, "Coupling", "changeparameters_couple.knob"), "w")
    h = open(os.path.join(_InputData.output_path, "Coupling", "changeparameters_couple.madx"), "w")
    g.write("* NAME  DELTA \n")
    g.write("$ %s    %le   \n")
    for vcorr in vcorrs:
        print >> g, "K"+vcorr, vcorrs[vcorr]
        print >> h, vcorr, "->KICK:=", vcorrs[vcorr], ";"
    h.write('return;')
    g.close()
    h.close()



def _handle_data_for_lhc():
    #.knob should always exist to be sent to LSA!
    src = os.path.join(os.path.join(_InputData.output_path, "changeparameters_couple.tfs"))
    dst = os.path.join(os.path.join(_InputData.output_path, "changeparameters_couple.knob"))
    Utilities.iotools.copy_item(src, dst)

    ##### for bumps
    if "bumps" in _InputData.variables_list:
        print "passing trough bumps loop"
        v = metaclass.twiss(os.path.join(_InputData.output_path, "changeparameters_couple.tfs"))  # @UnusedVariable
        # v will be used in Bumps.py. Don't do such things again!(vimaier)
        Utilities.iotools.delete_item(os.path.join(_InputData.output_path, "changeparameters_couple.tfs"))
        corrs = None # Will be assigned in exec Bumps.py (vimaier)
        execfile(os.path.join(_InputData.path_to_optics_files_dir, "Bumps.py"))
        execfile(os.path.join(_InputData.path_to_optics_files_dir, "mydictionary.py"))
        tfs_file_writer = Utilities.tfs_file_writer.TfsFileWriter("changeparameters_couple.tfs", _InputData.output_path)
        tfs_file_writer.add_column_names("NAME DELTA".split())
        tfs_file_writer.add_column_datatypes("%s %le".split())
        for vcorr in corrs:
            tfs_file_writer.add_table_row([vcorr, corrs[vcorr]])
        tfs_file_writer.write_to_file()
    #####

    v = metaclass.twiss(os.path.join(_InputData.output_path, "changeparameters_couple.tfs"))
    mad_script = open(os.path.join(_InputData.output_path, "changeparameters_couple.madx"), "w")
    names = getattr(v, "NAME", [])
    delta = getattr(v, "DELTA", [])

    for i in range(len(names)):
        if "bumps" in _InputData.variables_list:
            if cmp(delta[i], 0) == 1:
                mad_script.write(names[i] + "->KICK:=" + str(delta[i]) + ";\n")
            else:
                mad_script.write(names[i] + "->KICK:=" + str(delta[i]) + ";\n")
        else:
            if cmp(delta[i], 0) == 1:
                mad_script.write(names[i] + " = " + names[i] + " + " + str(delta[i]) + ";\n")
            else:
                mad_script.write(names[i] + " = " + names[i] + " " + str(delta[i]) + ";\n")

    mad_script.write("return;")
    mad_script.close()


#=======================================================================================================================
# helper class for script arguments
#=======================================================================================================================
class _InputData(object):
    """ Static class to access user input parameter. Necessary parameters will be checked. """
    output_path = ""
    accel_path = ""
    singular_value_cut = 0.0
    error_cut_c = 0.0
    error_cut_d = 0.0
    model_cut_c = 0.0
    model_cut_d = 0.0
    min_strength = 0.0 # Not used in code except of a print(vimaier)
    weights_list = []
    path_to_optics_files_dir = ""
    variables_list = []


    @staticmethod
    def static_init(output_path, accel, singular_value_cut, errorcut, modelcut, beta_beat_root, min_strength, weights_on_corrections, path_to_optics_files_dir, variables):
        if not Utilities.iotools.dirs_exist(output_path):
            raise ValueError("Output path does not exists. It has to contain getcouple[_free].out and getDy.out(when last flag in weights is 1.")
        _InputData.output_path = output_path

        _InputData._check_and_set_cuts(singular_value_cut, errorcut, modelcut)
        _InputData._check_and_set_accel_path(beta_beat_root, accel)

        if not Utilities.math.can_str_be_parsed_to_number(min_strength):
            raise ValueError("Given min strength is not a number: "+min_strength)
        _InputData.min_strength = float(min_strength)

        if re.match("^[01],[01],[01],[01],[01]$", weights_on_corrections) is None:
            raise ValueError("Wrong syntax of weigths: "+weights_on_corrections)
        weights = weights_on_corrections.split(',')
        _InputData.weights_list = [int(weights[0]), int(weights[1]), int(weights[2]), int(weights[3]), int(weights[4])]

        if not Utilities.iotools.dirs_exist(path_to_optics_files_dir):
            raise ValueError("Given path to optics files does not exist: "+path_to_optics_files_dir)
        _InputData.path_to_optics_files_dir = path_to_optics_files_dir
        _InputData.variables_list = variables.split(",")

        if PRINT_DEBUG:
            _InputData.print_input_data()


    @staticmethod
    def _check_and_set_cuts(singular_value_cut, errorcut, modelcut):
        if not Utilities.math.can_str_be_parsed_to_number(singular_value_cut):
            raise ValueError("Given cut is not a number: "+singular_value_cut)
        modelcuts = modelcut.split(",")
        if not Utilities.math.can_str_be_parsed_to_number(modelcuts[0]):
            raise ValueError("Given model cut is not a number: "+modelcuts[0]+" from "+modelcuts)
        if not Utilities.math.can_str_be_parsed_to_number(modelcuts[1]):
            raise ValueError("Given model cut is not a number: "+modelcuts[1]+" from "+modelcuts)
        errorcuts = errorcut.split(",")
        if not Utilities.math.can_str_be_parsed_to_number(errorcuts[0]):
            raise ValueError("Given error cut is not a number: "+errorcuts[0]+" from "+errorcuts)
        if not Utilities.math.can_str_be_parsed_to_number(errorcuts[1]):
            raise ValueError("Given error cut is not a number: "+errorcuts[1]+" from "+errorcuts)
        _InputData.singular_value_cut = float(singular_value_cut)
        _InputData.error_cut_c = float(errorcuts[0])
        _InputData.error_cut_d = float(errorcuts[1])
        _InputData.model_cut_c = float(modelcuts[0])
        _InputData.model_cut_d = float(modelcuts[1])

    @staticmethod
    def _check_and_set_accel_path(beta_beat_root, accel):
        if not Utilities.iotools.dirs_exist(beta_beat_root):
            print >> sys.stderr, "Given Beta-Beat.src path does not exist: "+beta_beat_root
            beta_beat_root = Utilities.iotools.get_absolute_path_to_betabeat_root()
            print >> sys.stderr, "Will take current Beta-Beat.src: "+beta_beat_root
        if accel not in ("LHCB1", "LHCB2", "SPS", "RHIC"):
            raise ValueError("Unknown/not supported accelerator: "+accel)
        if "LHC" in accel:
            accel_name = "LHCB"
        else:
            accel_name = "SPS"
        accel_path = os.path.join(beta_beat_root, "MODEL", accel_name, "fullresponse", accel)
        if not Utilities.iotools.dirs_exist(accel_path):
            raise ValueError("Acclelerator path does not exist: "+accel_path)
        _InputData.accel_path = accel_path

    @staticmethod
    def is_dy_switch_set():
        return 1 == _InputData.weights_list[4]

    @staticmethod
    def print_input_data():
        print "---------------------_InputData"
        print "Path to measurements:", _InputData.output_path
        print "Path to Accelerator model", _InputData.accel_path
        print "Path to optics files:", _InputData.path_to_optics_files_dir
        print "Minimum corrector strength", _InputData.min_strength
        print "Variables", _InputData.variables_list
        print "Singular value cut", _InputData.singular_value_cut
        print "Error cut C", _InputData.error_cut_c
        print "Error cut D", _InputData.error_cut_d
        print "Model cut C", _InputData.model_cut_c
        print "Model cut D", _InputData.model_cut_d
        print "Weights for correctors", _InputData.weights_list
        print "------------------------------"

    def __init__(self):
        raise NotImplementedError("static class _InputData cannot be instantiated")


#=======================================================================================================================
# main invocation
#=======================================================================================================================
def _start():
    options = _parse_args()
    main(
         output_path=options.path,
         accel=options.ACCEL,
         singular_value_cut=options.cut,
         errorcut=options.errorcut,
         modelcut=options.modelcut,
         beta_beat_root=options.rpath,
         min_strength=options.MinStr,
         weights_on_corrections=options.Dy,
         path_to_optics_files_dir=options.opt,
         variables=options.var
         )

if __name__ == "__main__":
    _start()
