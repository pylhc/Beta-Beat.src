r"""
.. module: Correction.correct

Created on ??

Single beam correction of phase, beta and horizontal dispersion.
TODO: get better description (vimaier)

Usage example::

    python correct_coupleDy.py  --accel=LHCB1
                                --Variables=Q
                                --ncorr=5
                                --rpath=C:\eclipse_4_2_2_python\workspace\Beta-Beat.src
                                --modelcut=0.3,0.3
                                --weight=1,1,0,0,0,10
                                --MinStr=0.00003
                                --path=C:\eclipse_4_2_2_python\workspace\Beta-Beat.src\Correction\test\correct\data\input\run1
                                --errorcut=0.5,0.5
                                --cut=0.01
                                --tech=SVD
                                --opt=C:\eclipse_4_2_2_python\workspace\Beta-Beat.src\Correction\test\correct\data\input

Hint: JustOneBeam is not used in the script except of a print. The reason is possibly the author wanted to have the same set
of arguments for all correction scripts due to GUI compatibility(vimaier).

Options::

      -h, --help            show this help message and exit
      -a ACCEL, --accel=ACCEL
                            What accelerator: LHCB1 LHCB2 SPS RHIC
      -t TECH, --tech=TECH  Which algorithm: SVD MICADO
      -n NCORR, --ncorr=NCORR
                            Number of Correctors for MICADO
      -p PATH, --path=PATH  Path to experimental files
      -c CUT, --cut=CUT     Singular value cut for the generalized inverse
      -e ERRORCUT, --errorcut=ERRORCUT
                            Maximum error allowed for the phase and dispersion
                            measurements, separated by commas; e.g. -e 0.013,0.2
      -m MODELCUT, --modelcut=MODELCUT
                            Maximum difference allowed between model and measured
                            phase and dispersion, separated by commas; e.g. -e
                            0.02,0.2
      -r RPATH, --rpath=RPATH
                            Path to BetaBeat repository (default is the afs
                            repository)
      -o opt, --optics=opt  optics
      -s MinStr, --MinStr=MinStr
                            Minimum strength of correctors in SVD correction
                            (default is 1e-6)
      -v var, --Variables=var
                            variables split with ,
      -j JUSTONEBEAM, --JustOneBeam=JUSTONEBEAM
                            0 Just quads from one beam are used, 1 all quads
                            (default is 0)
      -w WGT, --weight=WGT  Weighting factor (phasex, phasey, betax, betay,
                            dispersion, tunes)

.. moduleauthor:: Unknown
"""
import sys
import os
import optparse
import json
import pickle
import re


import __init__ # @UnusedImport init will include paths
import Python_Classes4MAD.GenMatrix
import Python_Classes4MAD.BCORR
import Python_Classes4MAD.metaclass
import Utilities.iotools
import Utilities.math

PRINT_DEBUG = False or sys.flags.debug # Change to 'True or...' or invoke python with -d option to activate further prints

#=======================================================================================================================
# _parse_args()-function
#=======================================================================================================================
def _parse_args():
    ''' Parses the arguments, checks for valid input and returns options '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                     help="What accelerator: LHCB1 LHCB2 SPS RHIC",
                     metavar="ACCEL", default="LHCB1",dest="ACCEL")
    parser.add_option("-t", "--tech",
                     help="Which algorithm: SVD MICADO",
                     metavar="TECH", default="SVD",dest="TECH")
    parser.add_option("-n", "--ncorr",
                     help="Number of Correctors for MICADO",
                     metavar="NCORR", default=5,dest="ncorr")
    parser.add_option("-p", "--path",
                     help="Path to experimental files",
                     metavar="PATH", default="./",dest="path")
    parser.add_option("-c", "--cut",
                      help="Singular value cut for the generalized inverse",
                      metavar="CUT", default=0.1 , dest="cut")
    parser.add_option("-e", "--errorcut",
                      help="Maximum error allowed for the phase and dispersion measurements, separated by commas; e.g. -e 0.013,0.2",
                      metavar="ERRORCUT", default="0.013,0.2" , dest="errorcut")
    parser.add_option("-m", "--modelcut",
                      help="Maximum difference allowed between model and measured phase and dispersion, separated by commas; e.g. -e 0.02,0.2",
                      metavar="MODELCUT", default="0.02,0.2" , dest="modelcut")
    parser.add_option("-r", "--rpath",
                      help="Path to BetaBeat repository (default is the afs repository)",
                      metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
    parser.add_option("-o", "--optics",
                      help="optics",
                      metavar="opt", default="nominal.opt" , dest="opt")
    parser.add_option("-s", "--MinStr",
                      help="Minimum strength of correctors in SVD correction (default is 1e-6)",
                      metavar="MinStr", default=0.000001 , dest="MinStr")
    parser.add_option("-v", "--Variables",
                      help="variables split with ,",
                      metavar="var", default="MQTb1" , dest="var")
    parser.add_option("-j", "--JustOneBeam",
                      help="0 Just quads from one beam are used, 1 all quads (default is 0)",
                      metavar="JUSTONEBEAM", default=0 , dest="JustOneBeam")
    parser.add_option("-w", "--weight",
                    help="Weighting factor (phasex, phasey, betax, betay, dispersion, tunes)",
                    metavar="WGT", default="1,1,0,0,1,10", dest="WGT")
    (options, args) = parser.parse_args() # @UnusedVariable no args
    return options

#=======================================================================================================================
# main()-function
#=======================================================================================================================
def main(
         output_path,
         accel="LHCB1",
         singular_value_cut=0.1,
         errorcut="0.013,0.2",
         modelcut="0.02,0.2",
         beta_beat_root=Utilities.iotools.get_absolute_path_to_betabeat_root(),
         min_strength=0.000001,
         weights_on_corrections="1,1,0,0,1,10",
         path_to_optics_files_dir="nominal.opt",
         variables="MQTb1",
         index_of_num_of_beams_in_gui=0,
         num_of_correctors=5,
         algorithm="SVD"
         ):

    _InputData.static_init(output_path, accel, singular_value_cut, errorcut, modelcut, beta_beat_root, min_strength,
                           weights_on_corrections, path_to_optics_files_dir, variables, index_of_num_of_beams_in_gui,
                           num_of_correctors, algorithm)


    _generate_changeparameters()

    _handle_data_for_accel(accel)


#=======================================================================================================================
# helper functions
#=======================================================================================================================
def _generate_changeparameters():
    full_response, phase_x, phase_y, beta_x, beta_y, dx, varslist = _load_input_files()

    phasexlist = Python_Classes4MAD.GenMatrix.MakePairs(phase_x, full_response['0'], modelcut=_InputData.model_cut, errorcut=_InputData.error_cut)
    phaseylist = Python_Classes4MAD.GenMatrix.MakePairs(phase_y, full_response['0'], modelcut=_InputData.model_cut, errorcut=_InputData.error_cut)
    betaxlist = make_beta_list(beta_x, full_response['0'], modelcut=_InputData.model_cut, errorcut=_InputData.error_cut)
    betaylist = make_beta_list(beta_y, full_response['0'], modelcut=_InputData.model_cut, errorcut=_InputData.error_cut)
    displist = Python_Classes4MAD.GenMatrix.MakeList(dx, full_response['0'])
    print "Input ready"

    beat_inp = Python_Classes4MAD.GenMatrix.beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, _InputData.weights_list)
    sensitivity_matrix = beat_inp.computeSensitivityMatrix(full_response) # @UnusedVariable sensitivity_matrix will be stored in beat_inp

    if _InputData.algorithm == "SVD":
        [deltas, varslist] = Python_Classes4MAD.GenMatrix.correctbeatEXP(phase_x, phase_y, dx, beat_inp, cut=_InputData.singular_value_cut, app=0, path=_InputData.output_path, xbet=beta_x, ybet=beta_y)
        iteration = 0 # Let's remove too low useless correctors
        while len([x for x in deltas if abs(x) < _InputData.min_strength]) > 0:
            iteration += 1
            il = len(varslist)
            varslist_t = []
            for i in range(il):
                if (abs(deltas[i]) > _InputData.min_strength):
                    varslist_t.append(varslist[i])

            varslist = varslist_t
            if len(varslist) == 0:
                print >> sys.stderr, "You want to correct with too high cut on the corrector strength"
                sys.exit(1)
            beat_inp = Python_Classes4MAD.GenMatrix.beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, _InputData.weights_list)
            sensitivity_matrix = beat_inp.computeSensitivityMatrix(full_response)  # @UnusedVariable sensitivity_matrix will be stored in beat_inp
            [deltas, varslist] = Python_Classes4MAD.GenMatrix.correctbeatEXP(phase_x, phase_y, dx, beat_inp, cut=_InputData.singular_value_cut, app=0, path=_InputData.output_path, xbet=beta_x, ybet=beta_y)
            print "Initial correctors:", il, ". Current: ", len(varslist), ". Removed for being lower than:", _InputData.min_strength, "Iteration:", iteration
        if PRINT_DEBUG:
            print deltas

    if _InputData.algorithm == "MICADO":
        Python_Classes4MAD.BCORR.bNCorrNumeric(phase_x, phase_y, dx, beat_inp, cut=_InputData.singular_value_cut, ncorr=_InputData.num_of_correctors, app=0, path=_InputData.output_path, beta_x=beta_x, beta_y=beta_y)


def _handle_data_for_accel(accel):
    if accel == "SPS":
        b = Python_Classes4MAD.metaclass.twiss(_InputData.output_path + "/changeparameters.tfs")
        corrs = []
        corrsYASP = []
        execfile(os.path.join(_InputData.accel_path, "Bumps.py")) # LOADS corrs
        execfile(os.path.join(_InputData.accel_path, "BumpsYASP.py")) # LOADS corrsYASP

        f = open(os.path.join(_InputData.output_path, "changeparameters.yasp"), "w") #Output for YASP...
        g = open(os.path.join(_InputData.output_path, "changeparameters.knob"), "w") #Output for Knob...
        f.write("#PLANE H\n")
        f.write("#UNIT RAD\n")
        g.write("* NAME  DELTA \n")
        g.write("$ %s    %le   \n")
        for corr in corrsYASP:
            print >> f, "#SETTING", corr, corrsYASP[corr]

        for corr in corrs:
            print >> g, "K" + corr, corrs[corr]

        f.close()
        g.close()
    if "LHC" in accel: #.knob should always exist to be sent to LSA!
        src = os.path.join(os.path.join(_InputData.output_path, "changeparameters.tfs"))
        dst = os.path.join(os.path.join(_InputData.output_path, "changeparameters.knob"))
        Utilities.iotools.copy_item(src, dst) # madx table
        b = Python_Classes4MAD.metaclass.twiss(os.path.join(_InputData.output_path, "changeparameters.tfs"))
        mad_script = open(os.path.join(_InputData.output_path, "changeparameters.madx"), "w")
        names = getattr(b, "NAME", [])
        delta = getattr(b, "DELTA", [])
        for i in range(len(names)):
            if cmp(delta[i], 0) == 1:
                mad_script.write(names[i] + " = " + names[i] + " + " + str(delta[i]) + ";\n")
            else:
                mad_script.write(names[i] + " = " + names[i] + " " + str(delta[i]) + ";\n")

        mad_script.write("return;")
        mad_script.close()


def _load_input_files():
    print "Starting loading Full Response optics"
    full_response = pickle.load(open(_InputData.path_to_optics_files_dir + '/FullResponse', 'rb'))
    print "Loading ended"

    path_to_phase_x = os.path.join(_InputData.output_path, "getphasex_free.out")
    path_to_phase_y = os.path.join(_InputData.output_path, "getphasey_free.out")
    if not os.path.exists(path_to_phase_x) or not os.path.exists(path_to_phase_y):
        path_to_phase_x = os.path.join(_InputData.output_path, "getphasex.out")
        path_to_phase_y = os.path.join(_InputData.output_path, "getphasey.out")

    path_to_beta_x = os.path.join(_InputData.output_path, "getbetax_free.out")
    path_to_beta_y = os.path.join(_InputData.output_path, "getbetay_free.out")
    if not os.path.exists(path_to_beta_x) or not os.path.exists(path_to_beta_y):
        path_to_beta_x = os.path.join(_InputData.output_path, "getbetax.out")
        path_to_beta_y = os.path.join(_InputData.output_path, "getbetay.out")

    phase_x = Python_Classes4MAD.metaclass.twiss(path_to_phase_x)
    phase_y = Python_Classes4MAD.metaclass.twiss(path_to_phase_y)
    beta_x = Python_Classes4MAD.metaclass.twiss(path_to_beta_x)
    beta_y = Python_Classes4MAD.metaclass.twiss(path_to_beta_y)

    path_to_ndx = os.path.join(_InputData.output_path, "getNDx.out")
    if os.path.exists(path_to_ndx):
        dx = Python_Classes4MAD.metaclass.twiss(path_to_ndx)
    else:
        print "WARNING: No good dispersion or inexistent file getDx"
        print "WARNING: Correction will not take into account NDx"
        dx = []

    _add_int_part_of_tunes_to_exp_data(full_response, phase_x, phase_y)

    # Load vars in AllLists
    var_source_file_path = os.path.join(_InputData.accel_path, "AllLists.json")
    knobsdict = json.load(file(var_source_file_path, 'r'))
    # extra depdency to be able to handle to different magnets group
    varslist = []
    for var in _InputData.variables_list:
        try:
            varslist += knobsdict[var]
        except KeyError:
            print >> sys.stderr, "Given variable({0}) not in dictionary({1})".format(var, var_source_file_path)

    return full_response, phase_x, phase_y, beta_x, beta_y, dx, varslist
def _add_int_part_of_tunes_to_exp_data(full_response, phase_x, phase_y):
    intqx = int(full_response['0'].Q1)
    intqy = int(full_response['0'].Q2)
    tune_x = getattr(phase_x, "Q1")
    tune_y = getattr(phase_y, "Q2")
    if tune_x < 0.0:
        tune_x += 1.0
    if tune_y < 0.0:
        tune_y += 1.0
    setattr(phase_x, "Q1", tune_x + intqx)
    setattr(phase_y, "Q2", tune_y + intqy)
    print "Integer part of tunes: ", intqx, intqy
    print "Experiment tunes: ", getattr(phase_x, "Q1"), getattr(phase_y, "Q2")


def  make_beta_list(x, m, modelcut=40, errorcut=20):   # Errors are in meters
    t = []
    cou = 0
    keys = x.__dict__.keys()
    if "BETY" in keys:
        bmdl = "BETYMDL"
        STD = x.STDBETY
        BET = x.BETY
    else:
        bmdl = "BETXMDL"
        STD = x.STDBETX
        BET = x.BETX
    print "Number of x BPMs", len(x.NAME)
    for i in range(len(x.NAME)):
        bm = x.__dict__[bmdl][i]
        if (STD[i] < errorcut and abs(BET[i]-bm) < modelcut):
            try:
                m.indx[x.NAME[i].upper()]
            except KeyError:
                print "Not in Response:", x.NAME[i].upper()
                cou += 1
            else:
                t.append(x.NAME[i])
        else:
            cou += 1
    if cou > 0:
        print "Warning in make_beta_list: ", cou, "BPM  removed from data for not beeing in the model or having too large error deviations: ", bmdl, modelcut, "STDPH", errorcut, "LEN", len(t)
    return t

#=======================================================================================================================
# helper class for script arguments
#=======================================================================================================================
class _InputData(object):
    """ Static class to access user input parameter. Necessary parameters will be checked. """
    output_path = ""
    accel_path = ""
    singular_value_cut = 0.0
    error_cut = 0.0
    error_cut_dx = 0.0
    model_cut = 0.0
    model_cut_dx = 0.0
    min_strength = 0.0
    weights_list = []
    path_to_optics_files_dir = ""
    variables_list = []
    use_two_beams = False # not used in code except of a print (vimaier)
    num_of_correctors = 0
    algorithm = ""


    @staticmethod
    def static_init(output_path, accel, singular_value_cut, errorcut, modelcut, beta_beat_root, min_strength,
                    weights_on_corrections, path_to_optics_files_dir, variables, index_of_num_of_beams_in_gui,
                    num_of_correctors, algorithm):
        if not Utilities.iotools.dirs_exist(output_path):
            raise ValueError("Output path does not exists. It has to contain getcouple[_free].out and getDy.out(when last flag in weights is 1.")
        _InputData.output_path = output_path

        _InputData._check_and_set_cuts(singular_value_cut, errorcut, modelcut)
        _InputData._check_and_set_accel_path(beta_beat_root, accel)

        if not Utilities.math.can_str_be_parsed_to_number(min_strength):
            raise ValueError("Given min strength is not a number: "+min_strength)
        if "SPS" != accel:
            _InputData.min_strength = float(min_strength)

        if re.match("^[01],[01],[01],[01],[01],\\d\\d$", weights_on_corrections) is None:
            raise ValueError("Wrong syntax of weigths: "+weights_on_corrections)
        weights = weights_on_corrections.split(',')
        _InputData.weights_list = [int(weights[0]), int(weights[1]), int(weights[2]), int(weights[3]), int(weights[4]), int(weights[5])]

        if not Utilities.iotools.dirs_exist(path_to_optics_files_dir):
            raise ValueError("Given path to optics files does not exist: "+path_to_optics_files_dir)
        _InputData.path_to_optics_files_dir = path_to_optics_files_dir
        _InputData.variables_list = variables.split(",")

        _InputData.use_two_beams = 1 == index_of_num_of_beams_in_gui # 0 index for 'One Beam'; 1 index for 'Two Beams'(vimaier)
        if not Utilities.math.can_str_be_parsed_to_number(num_of_correctors):
            raise ValueError("Given number of correctors is not a number: "+num_of_correctors)
        _InputData.num_of_correctors = int(num_of_correctors)
        if "SVD" != algorithm and "MICADO" != algorithm:
            raise ValueError("Wrong algorithm given: "+algorithm)
        _InputData.algorithm = algorithm



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
        _InputData.error_cut = float(errorcuts[0])
        _InputData.error_cut_dx = float(errorcuts[1])
        _InputData.model_cut = float(modelcuts[0])
        _InputData.model_cut_dx = float(modelcuts[1])

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
        print "Path to Accelerator model:", _InputData.accel_path
        print "Path to optics files:", _InputData.path_to_optics_files_dir
        print "Minimum corrector strength:", _InputData.min_strength
        print "Variables:", _InputData.variables_list
        print "Singular value cut:", _InputData.singular_value_cut
        print "Error cut:", _InputData.error_cut
        print "Error cut Dx:", _InputData.error_cut_dx
        print "Model cut:", _InputData.model_cut
        print "Model cut Dx:", _InputData.model_cut_dx
        print "Weights for correctors:", _InputData.weights_list
        print "Use two beams?:", _InputData.use_two_beams
        print "Number of correctors:", _InputData.num_of_correctors
        print "Chosen algorithm:", _InputData.algorithm
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
         weights_on_corrections=options.WGT,
         path_to_optics_files_dir=options.opt,
         variables=options.var,
         index_of_num_of_beams_in_gui=options.JustOneBeam,
         num_of_correctors=options.ncorr,
         algorithm=options.TECH
         )

if __name__ == "__main__":
    _start()
