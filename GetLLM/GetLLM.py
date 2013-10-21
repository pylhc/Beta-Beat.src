r'''
.. module: GetLLM.GetLLM

Created on 11/09/09

:author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)

:version: 3.00dev


GetLLM calculates a large collection of optics functions and parameters at the BPMs using the output from DRIVE.

The GetLLM output is distributed over different files according to the observable and the transverse plane:
 - getCO*; files containing the closed orbit at all BPMs, one file per plane.
 - getIP*; files containing beta and phase advance functions at the IPs.
 - getampbeta*; files containing the beta function as compputed from amplitude, one file per plane (in case of ACdipole free and driven data is also computed).
 - getbeta*; files containing the beta function as computed from phase advances between 3 BPMs. There is one file per plane (in case of AC-dipole free and driven data is computed).
 - getchi*; files containing the $\chi$ terms.
 - getcouple*; files containing the coupling resonances $f_{1001}$  and $f_{0101}$ (in case of AC-dipole free and driven data is computed).
 - getsex*; files containing the sextupolar resonance driving terms.
 - getphase*; files containing the phase advances between BPMs, one file per plane (in case of AC-dipole free and driven data is computed).
 - getphasetot*; files containing the total phase-advance using the first BPM as reference, one file per plane (in case of ACdipole free and driven data is computed).
 - getkick*; files containing the kick amplitudes and the linear invariants as derived from peak-to-peak values and spectral lines amplitudes.
 - getD*; files containing the dispersion, one per plane ( if off-momentum data is acquired).
 - getNDx*; files containing the normalized horizontal dispersion ( if off-momentum data is acquired).


Usage1::

    >pythonafs ../GetLLM.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./

Usage2::

    >pythonafs ../GetLLM.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./


Change history::

        git log GetLLM.py

'''
import sys
import traceback
import math


import __init__ # @UnusedImport init will include paths
import Python_Classes4MAD.metaclass
import utils.tfs_file
import algorithms.helper
import algorithms.phase
import algorithms.beta
import algorithms.compensate_ac_effect
import algorithms.dispersion
import algorithms.coupling
import algorithms.interaction_point
import algorithms.chi_terms
import Utilities.iotools



####
#######
#########
VERSION = 'V3.0.0 Dev'
#########
#######
####
DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# _parse_args()-function
#===================================================================================================
def _parse_args():
    ''' Parses command line arguments. '''
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-a", "--accel",
                    help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
                    metavar="ACCEL", default="LHCB1",dest="ACCEL")
    parser.add_option("-d", "--dictionary",
                    help="File with the BPM dictionary",
                    metavar="DICT", default="0", dest="dict")
    parser.add_option("-m", "--model",
                    help="Twiss free model file *.dat. For free oscillations, ONLY the file *.dat should be present in the path. For AC dipole, the path should also contain a file *_ac.dat (BOTH files a needed in this case).",
                    metavar="TwissFile", default="0", dest="Twiss")
    parser.add_option("-f", "--files",
                    help="Files from analysis, separated by comma",
                    metavar="TwissFile", default="0", dest="files")
    parser.add_option("-o", "--output",
                    help="Output Path",
                    metavar="OUT", default="./", dest="output")
    parser.add_option("-c", "--cocut",
                    help="Cut for closed orbit measurement [um]",
                    metavar="COCUT", default=4000, dest="COcut")
    parser.add_option("-n", "--nbcpl",
                    help="Analysis option for coupling, 1 bpm or 2 bpms",
                    metavar="NBCPL", default=2, dest="NBcpl")
    parser.add_option("-t", "--tbtana",
                    help="Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA",
                    metavar="TBTANA", default="SUSSIX", dest="TBTana")
    parser.add_option("-b", "--bpmu",
                    help="BPMunit: um, mm, cm, m (default um)",
                    metavar="BPMUNIT", default="um", dest="BPMUNIT")
    parser.add_option("-l", "--nonlinear",
                    help="Switch to output higher order resonance stuffs, on=1(default)/off=0",
                    metavar="HIGHER", default="1" , dest="higher")
    parser.add_option("-p", "--lhcphase",
                    help="Compensate phase shifts by tunes for the LHC experiment data, off=0(default)/on=1",
                    metavar="LHCPHASE", default="0" , dest="lhcphase")
    parser.add_option("-k", "--bbthreshold",
                    help="Set beta-beating threshold for action calculations, default = 0.15",
                    metavar="BBTHRESH", default="0.15" , dest="bbthreshold")
    parser.add_option("-e", "--errthreshold",
                    help="Set beta relative uncertainty threshold for action calculations, default = 0.1",
                    metavar="ERRTHRESH", default="0.15" , dest="errthreshold")
    parser.add_option("-g", "--threebpm",
                    help="Forces to use the 3 BPM method, yes=1/no=0, default = 0",
                    metavar="USE_ONLY_THREE_BPMS_FOR_BETA_FROM_PHASE", default="0" , dest="use_only_three_bpms_for_beta_from_phase")

    # Take index 0 since index 1(args) is not used (vimaier)
    options = parser.parse_args()[0]
    return options


#===================================================================================================
# main()-function
#===================================================================================================
def main(
         outputpath,
         files_to_analyse,
         model_filename,
         dict_file="0",
         accel="LHCB1",
         lhcphase="0",
         BPMU="um",
         COcut=4000,
         NBcpl=2,
         TBTana="SUSSIX",
         higher_order=1,
         bbthreshold="0.15",
         errthreshold="0.15",
         use_only_three_bpms_for_beta_from_phase=0
         ):
    '''
    GetLLM main function.

    :param string outputpath: The output path to store results
    :param string files_to_analyse: List of files, comma separated string.
    :param dict_file: Name of the script which will be executed. Should store dictionary with
                    mappings of BPM names.
    :param string accel: Type of accelerator. LHCB1, LHCB2, LHCB4, RHIC, SPS
    :param string lhcphase: "0" or "1" -- Compensate phase shifts by tunes for the LHC experiment data,
                             off=0(default)/on=1
    :param string BPMU: BPMunit: "um", "mm", "cm", "m" (default "um")
    :param int COcut: Cut for closed orbit measurement [um]
    :param int NBcpl: For selecting the coupling measurement method 1 bpm or 2 bpms
    :param string TBTana: Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA
    :param int higher_order': output higher order resonance stuff, on=1(default)/off=0

    :returns: int  -- 0 if the function run successfully otherwise !=0.
    '''
    return_code = 0
    print "Starting GetLLM ", VERSION

    # The following objects stores multiple variables for GetLLM to avoid having much local
    # variables. Identifiers supposed to be as short as possible.
    # --vimaier
    getllm_d = _GetllmData()
    twiss_d = _TwissData()
    tune_d = _TuneData()

    getllm_d, mad_twiss, mad_ac, bpm_dictionary, mad_elem = _intial_setup(getllm_d,
                                                                        outputpath,
                                                                        model_filename,
                                                                        dict_file,
                                                                        accel,
                                                                        BPMU,
                                                                        COcut,
                                                                        lhcphase,
                                                                        NBcpl)

    files_dict = _create_tfs_files(getllm_d, model_filename)

    twiss_d, files_dict = _analyse_src_files(getllm_d, twiss_d, files_to_analyse, TBTana, files_dict)


    tune_d.initialize_tunes(getllm_d.with_ac_calc, mad_twiss, mad_ac, twiss_d)

    # Construct pseudo-double plane BPMs
    if (getllm_d.accel=="SPS" or "RHIC" in getllm_d.accel) and twiss_d.has_zero_dpp_x() and twiss_d.has_zero_dpp_y():
        [pseudo_list_x, pseudo_list_y] = algorithms.helper.pseudo_double_plane_monitors(mad_twiss, twiss_d.zero_dpp_x, twiss_d.zero_dpp_y, bpm_dictionary)
    else:
        # Initialize variables otherwise calculate_coupling would raise an exception(vimaier)
        pseudo_list_x = None
        pseudo_list_y = None


    #-------- Check monitor compatibility between data and model
    _check_bpm_compatibility(twiss_d, mad_twiss)

    try:
        #-------- START Phase
        phase_d, tune_d = algorithms.phase.calculate_phase(getllm_d, twiss_d, tune_d, mad_twiss, mad_ac, mad_elem, files_dict)

        #-------- START Total Phase
        algorithms.phase.calculate_total_phase(getllm_d, twiss_d, tune_d, phase_d, mad_twiss, mad_ac, files_dict)

        #-------- START Beta
        beta_d = algorithms.beta.calculate_beta_from_phase(getllm_d, twiss_d, tune_d, phase_d, mad_twiss, mad_ac, files_dict, use_only_three_bpms_for_beta_from_phase)

        #------- START beta from amplitude
        beta_d = algorithms.beta.calculate_beta_from_amplitude(getllm_d, twiss_d, tune_d, phase_d, beta_d, mad_twiss, mad_ac, files_dict)

        #-------- START IP
        algorithms.interaction_point.calculate_ip(getllm_d, twiss_d, tune_d, phase_d, beta_d, mad_twiss, mad_ac, files_dict)

        #-------- START Orbit
        list_of_co_x, list_of_co_y, files_dict = _calculate_orbit(getllm_d, twiss_d, tune_d, mad_twiss, files_dict)

        #-------- START Dispersion
        algorithms.dispersion.calculate_dispersion(getllm_d, twiss_d, tune_d, mad_twiss, files_dict, beta_d.x_amp, list_of_co_x, list_of_co_y)

        #-------- START coupling.
        tune_d = algorithms.coupling.calculate_coupling(getllm_d, twiss_d, phase_d, tune_d, mad_twiss, mad_ac, files_dict, pseudo_list_x, pseudo_list_y)

        #-------- Phase, Beta and coupling for non-zero DPP
        _phase_and_beta_for_non_zero_dpp(getllm_d, twiss_d, tune_d, phase_d, bpm_dictionary, mad_twiss, files_dict, pseudo_list_x, pseudo_list_y, use_only_three_bpms_for_beta_from_phase)

        if higher_order:
            if TBTana == "SUSSIX":
                #------ Start getsextupoles @ Glenn Vanbavinckhove
                files_dict = _calculate_getsextupoles(twiss_d, phase_d, mad_twiss, files_dict, tune_d.q1f)

                #------ Start getchiterms @ Glenn Vanbavinckhove
                files_dict = algorithms.chi_terms.calculate_chiterms(getllm_d, twiss_d, mad_twiss, files_dict)

            #------ Start get Q,JX,delta
            files_dict = _calculate_kick(getllm_d, twiss_d, tune_d, phase_d, beta_d, mad_twiss, mad_ac, files_dict, bbthreshold, errthreshold)
        else:
            print "Not analysing higher order..."
    except:
        traceback.print_exc()
        return_code = 1
    finally:
        print "Writing files"
        for tfsfile in files_dict.itervalues():
            tfsfile.write_to_file(formatted=True)

    return return_code
# END main() ---------------------------------------------------------------------------------------


#===================================================================================================
# helper-functions
#===================================================================================================
def _intial_setup(getllm_d, outputpath, model_filename, dict_file, accel, bpm_unit, cut_co, lhcphase, num_beams_cpl):
    getllm_d.set_outputpath(outputpath)
    getllm_d.set_bpmu_and_cut_for_closed_orbit(cut_co, bpm_unit)
    getllm_d.lhc_phase = lhcphase
    getllm_d.num_beams_for_coupling = num_beams_cpl

    if dict_file == "0":
        bpm_dictionary = {}
    else:
        execfile(dict_file)
        bpm_dictionary = dictionary # temporaryly since presently name is not bpm_dictionary

    # Beam direction
    getllm_d.beam_direction = 1
    getllm_d.accel = accel
    if accel == "LHCB2":
        getllm_d.beam_direction = -1 # THIS IS CORRECT, be careful with tune sign in SUSSIX and eigenmode order in SVD
    elif accel == "LHCB4":
        getllm_d.beam_direction = 1 # IS THIS CORRECT? I (rogelio) try for Simon...
        getllm_d.accel = "LHCB2" #This is because elements later are named B2 anyway, not B4

    #-- finding base model
    try:
        mad_twiss = Python_Classes4MAD.metaclass.twiss(model_filename, bpm_dictionary) # model from MAD : Twiss instance
        print "Base model found!"
    except IOError:
        print >> sys.stderr, "Twiss file loading failed for:\n\t", model_filename
        print >> sys.stderr, "Provide a valid model file."
        sys.exit(1)

    #-- finding the ac dipole model
    getllm_d.with_ac_calc = False
    try:
        mad_ac = Python_Classes4MAD.metaclass.twiss(model_filename.replace(".dat", "_ac.dat")) # model with ac dipole : Twiss instance
        getllm_d.with_ac_calc = True
        print "Driven Twiss file found. AC dipole effects calculated with the effective model (get***_free2.out)"
    except IOError:
        mad_ac = mad_twiss
        print "WARN: AC dipole effects not calculated. Driven twiss file does not exsist !"

    #-- Test if the AC dipole (MKQA) is in the model of LHC
    mad_elem = None
    if getllm_d.with_ac_calc:
        if 'LHC' in accel:
            if 'MKQA.6L4.' + accel[3:] in getattr(mad_twiss, "NAME", []):
                print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
            else:
                try:
                    mad_elem = Python_Classes4MAD.metaclass.twiss(model_filename.replace(".dat", "_elements.dat"))
                    print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
                except IOError:
                    print 'WARN: AC dipoles not in the model. AC dipole effects not calculated with analytic equations !'
        else:
            print 'WARN: AC dipole effects calculated with analytic equations only for LHC for now'

    return getllm_d, mad_twiss, mad_ac, bpm_dictionary, mad_elem
# END _intial_setup ---------------------------------------------------------------------------------


def _create_tfs_files(getllm_d, model_filename):
    '''
    Creates the most tfs files and stores it in an dictionary whereby the key represents the file
    and the value is the corresponding GetllmTfsFile.

    :Return: dict: string --> GetllmTfsFile
            A dictionary of created GetllmTfsFile objects. Keys are the filenames and values are the
            GetllmTfsFile objects.
    '''
    # Static variable of GetllmTfsFile to save the outputfile, GetLLM version and model filename
    utils.tfs_file.GetllmTfsFile.s_output_path = getllm_d.outputpath
    utils.tfs_file.GetllmTfsFile.s_getllm_version = VERSION
    utils.tfs_file.GetllmTfsFile.s_mad_filename = model_filename

    files_dict = {}
    files_dict['getphasex.out'] = utils.tfs_file.GetllmTfsFile('getphasex.out')
    files_dict['getphasey.out'] = utils.tfs_file.GetllmTfsFile('getphasey.out')
    files_dict['getphasetotx.out'] = utils.tfs_file.GetllmTfsFile('getphasetotx.out')
    files_dict['getphasetoty.out'] = utils.tfs_file.GetllmTfsFile('getphasetoty.out')
    if getllm_d.with_ac_calc:
        files_dict['getphasex_free.out'] = utils.tfs_file.GetllmTfsFile('getphasex_free.out')
        files_dict['getphasey_free.out'] = utils.tfs_file.GetllmTfsFile('getphasey_free.out')
        files_dict['getphasex_free2.out'] = utils.tfs_file.GetllmTfsFile('getphasex_free2.out')
        files_dict['getphasey_free2.out'] = utils.tfs_file.GetllmTfsFile('getphasey_free2.out')
        files_dict['getphasetotx_free.out'] = utils.tfs_file.GetllmTfsFile('getphasetotx_free.out')
        files_dict['getphasetoty_free.out'] = utils.tfs_file.GetllmTfsFile('getphasetoty_free.out')
        files_dict['getphasetotx_free2.out'] = utils.tfs_file.GetllmTfsFile('getphasetotx_free2.out')
        files_dict['getphasetoty_free2.out'] = utils.tfs_file.GetllmTfsFile('getphasetoty_free2.out')
    files_dict['getbetax.out'] = utils.tfs_file.GetllmTfsFile('getbetax.out')
    files_dict['getbetay.out'] = utils.tfs_file.GetllmTfsFile('getbetay.out')
    if getllm_d.with_ac_calc:
        files_dict['getbetax_free.out'] = utils.tfs_file.GetllmTfsFile('getbetax_free.out')
        files_dict['getbetay_free.out'] = utils.tfs_file.GetllmTfsFile('getbetay_free.out')
        files_dict['getbetax_free2.out'] = utils.tfs_file.GetllmTfsFile('getbetax_free2.out')
        files_dict['getbetay_free2.out'] = utils.tfs_file.GetllmTfsFile('getbetay_free2.out')
    files_dict['getampbetax.out'] = utils.tfs_file.GetllmTfsFile('getampbetax.out')
    files_dict['getampbetay.out'] = utils.tfs_file.GetllmTfsFile('getampbetay.out')
    if getllm_d.with_ac_calc:
        files_dict['getampbetax_free.out'] = utils.tfs_file.GetllmTfsFile('getampbetax_free.out')
        files_dict['getampbetay_free.out'] = utils.tfs_file.GetllmTfsFile('getampbetay_free.out')
        files_dict['getampbetax_free2.out'] = utils.tfs_file.GetllmTfsFile('getampbetax_free2.out')
        files_dict['getampbetay_free2.out'] = utils.tfs_file.GetllmTfsFile('getampbetay_free2.out')
    files_dict['getCOx.out'] = utils.tfs_file.GetllmTfsFile('getCOx.out')
    files_dict['getCOy.out'] = utils.tfs_file.GetllmTfsFile('getCOy.out')
    files_dict['getNDx.out'] = utils.tfs_file.GetllmTfsFile('getNDx.out')
    files_dict['getDx.out'] = utils.tfs_file.GetllmTfsFile('getDx.out')
    files_dict['getDy.out'] = utils.tfs_file.GetllmTfsFile('getDy.out')
    files_dict['getcouple.out'] = utils.tfs_file.GetllmTfsFile('getcouple.out')
    if getllm_d.with_ac_calc:
        files_dict['getcouple_free.out'] = utils.tfs_file.GetllmTfsFile('getcouple_free.out')
        files_dict['getcouple_free2.out'] = utils.tfs_file.GetllmTfsFile('getcouple_free2.out')
    files_dict['getcoupleterms.out'] = utils.tfs_file.GetllmTfsFile('getcoupleterms.out')
    if "LHC" in getllm_d.accel:
        files_dict['getIP.out'] = utils.tfs_file.GetllmTfsFile('getIP.out')
        files_dict['getIPx.out'] = utils.tfs_file.GetllmTfsFile('getIPx.out')
        files_dict['getIPy.out'] = utils.tfs_file.GetllmTfsFile('getIPy.out')
        files_dict['getIPfromphase.out'] = utils.tfs_file.GetllmTfsFile('getIPfromphase.out')
        if getllm_d.with_ac_calc:
            files_dict['getIPx_free.out'] = utils.tfs_file.GetllmTfsFile('getIPx_free.out')
            files_dict['getIPy_free.out'] = utils.tfs_file.GetllmTfsFile('getIPy_free.out')
            files_dict['getIPx_free2.out'] = utils.tfs_file.GetllmTfsFile('getIPx_free2.out')
            files_dict['getIPy_free2.out'] = utils.tfs_file.GetllmTfsFile('getIPy_free2.out')
            files_dict['getIPfromphase_free.out'] = utils.tfs_file.GetllmTfsFile('getIPfromphase_free.out')
            files_dict['getIPfromphase_free2.out'] = utils.tfs_file.GetllmTfsFile('getIPfromphase_free2.out')

    files_dict["getsex3000.out"] = utils.tfs_file.GetllmTfsFile("getsex3000.out")
    files_dict['getchi3000.out'] = utils.tfs_file.GetllmTfsFile('getchi3000.out')
    files_dict['getchi1010.out'] = utils.tfs_file.GetllmTfsFile('getchi1010.out')
    files_dict['getkick.out'] = utils.tfs_file.GetllmTfsFile('getkick.out')
    files_dict['getkickphase.out'] = utils.tfs_file.GetllmTfsFile('getkickphase.out')
    files_dict['getkickac.out'] = utils.tfs_file.GetllmTfsFile('getkickac.out')

    return files_dict
# END _create_tfs_files -----------------------------------------------------------------------------


def _analyse_src_files(getllm_d, twiss_d, files_to_analyse, turn_by_turn_algo, files_dict):

    if turn_by_turn_algo == "SUSSIX":
        suffix_x = '_linx'
        suffix_y = '_liny'
    elif turn_by_turn_algo == 'SVD':
        suffix_x = '_svdx'
        suffix_y = '_svdy'
    elif turn_by_turn_algo == 'HA':
        suffix_x = '_hax'
        suffix_y = '_hay'

    for file_in in files_to_analyse.split(','):
        # x file
        if file_in.endswith(".gz"):
            file_x = file_in.replace(".gz", suffix_x + ".gz")
        else:
            file_x = file_in + suffix_x

        twiss_file_x = None
        try:
            twiss_file_x = Python_Classes4MAD.metaclass.twiss(file_x)
            if twiss_file_x.has_no_bpm_data():
                print >> sys.stderr, "Ignoring empty file:", twiss_file_x.filename
                twiss_file_x = None
        except IOError:
            print >> sys.stderr, "Cannot load file:", file_x
        except ValueError:
            pass # Information printed by metaclass already

        if None != twiss_file_x:
            try:
                dppi = getattr(twiss_file_x, "DPP", 0.0)
            except AttributeError:
                dppi = 0.0
            if type(dppi) != float:
                print >> sys.stderr, 'Warning: DPP may not be given as a number in ', file_x, '...trying to forcibly cast it as a number'
                try:
                    dppi = float(dppi)
                    print 'dppi= ', dppi
                except ValueError:
                    print >> sys.stderr, 'but failing. DPP in ', file_x, ' is something wrong. String? --- leaving GetLLM'
                    print >> sys.stderr, traceback.format_exc()
                    sys.exit(1)
            if dppi == 0.0:
                twiss_d.zero_dpp_x.append(twiss_file_x)
                files_dict['getphasex.out'].add_filename_to_getllm_header(file_x)
                files_dict['getphasetotx.out'].add_filename_to_getllm_header(file_x)
                files_dict['getbetax.out'].add_filename_to_getllm_header(file_x)
                files_dict['getampbetax.out'].add_filename_to_getllm_header(file_x)
                files_dict['getCOx.out'].add_filename_to_getllm_header(file_x)
                files_dict['getNDx.out'].add_filename_to_getllm_header(file_x)
                files_dict['getDx.out'].add_filename_to_getllm_header(file_x)
                files_dict['getcouple.out'].add_filename_to_getllm_header(file_in)
                if "LHC" in getllm_d.accel:
                    files_dict['getIPx.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPy.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPfromphase.out'].add_filename_to_getllm_header(file_in)
                    if getllm_d.with_ac_calc:
                        files_dict['getIPx_free.out'].add_filename_to_getllm_header(file_in)
                        files_dict['getIPy_free.out'].add_filename_to_getllm_header(file_in)
                        files_dict['getIPx_free2.out'].add_filename_to_getllm_header(file_in)
                        files_dict['getIPy_free2.out'].add_filename_to_getllm_header(file_in)
                        files_dict['getIPfromphase_free.out'].add_filename_to_getllm_header(file_in)
                        files_dict['getIPfromphase_free2.out'].add_filename_to_getllm_header(file_in)
                if getllm_d.with_ac_calc:
                    files_dict['getphasex_free.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getphasex_free2.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getphasetotx_free.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getphasetotx_free2.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getbetax_free.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getbetax_free2.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getampbetax_free.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getampbetax_free2.out'].add_filename_to_getllm_header(file_x)
                    files_dict['getcouple_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getcouple_free2.out'].add_filename_to_getllm_header(file_in)
            else:
                twiss_d.non_zero_dpp_x.append(twiss_file_x)
                files_dict['getNDx.out'].add_filename_to_getllm_header(file_x)
                files_dict['getDx.out'].add_filename_to_getllm_header(file_x)

        # y file
        if file_in.endswith(".gz"):
            file_y = file_in.replace(".gz", suffix_y + ".gz")
        else:
            file_y = file_in + suffix_y

        twiss_file_y = None
        try:
            twiss_file_y = Python_Classes4MAD.metaclass.twiss(file_y)
            if twiss_file_y.has_no_bpm_data():
                print >> sys.stderr, "Ignoring empty file:", twiss_file_y.filename
                twiss_file_y = None
        except IOError:
            print 'Warning: There seems no ' + str(file_y) + ' file in the specified directory.'
        except ValueError:
            pass # Information printed by metaclass already

        if None != twiss_file_y:
            try:
                dppi = getattr(twiss_file_y, "DPP", 0.0)
            except AttributeError:
                dppi = 0.0
            if type(dppi) != float:
                print >> sys.stderr, 'Warning: DPP may not be given as a number in ', file_y, '...trying to forcibly cast it as a number'
                try:
                    dppi = float(dppi)
                    print 'dppi= ', dppi
                except ValueError:
                    print >> sys.stderr, 'but failing. DPP in ', file_y, ' is something wrong. String? --- leaving GetLLM'
                    print >> sys.stderr, traceback.format_exc()
                    sys.exit(1)
            if dppi == 0.0:
                twiss_d.zero_dpp_y.append(twiss_file_y)
                files_dict['getphasey.out'].add_filename_to_getllm_header(file_y)
                files_dict['getphasetoty.out'].add_filename_to_getllm_header(file_y)
                files_dict['getbetay.out'].add_filename_to_getllm_header(file_y)
                files_dict['getampbetay.out'].add_filename_to_getllm_header(file_y)
                files_dict['getCOy.out'].add_filename_to_getllm_header(file_y)
                files_dict['getDy.out'].add_filename_to_getllm_header(file_y)
                if getllm_d.with_ac_calc:
                    files_dict['getphasey_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasey_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasetoty_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasetoty_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getbetay_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getbetay_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getampbetay_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getampbetay_free2.out'].add_filename_to_getllm_header(file_y)
            else:
                twiss_d.non_zero_dpp_y.append(twiss_file_y)
                files_dict['getDy.out'].add_filename_to_getllm_header(file_y)


    if not twiss_d.has_zero_dpp_x():
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'

        if twiss_d.has_non_zero_dpp_x():
            twiss_d.zero_dpp_x = twiss_d.non_zero_dpp_x
            twiss_d.zero_dpp_y = twiss_d.non_zero_dpp_y
            twiss_d.non_zero_dpp_x = []
            twiss_d.non_zero_dpp_y = []

            print "Previous warning suppressed, running in chromatic mode"
            files_dict['getphasex.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getbetax.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getampbetax.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getCOx.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getNDx.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getDx.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getcouple.out'].add_filename_to_getllm_header("chrommode")
            if getllm_d.with_ac_calc:
                files_dict['getcouple_free.out'].add_filename_to_getllm_header("chrommode")
                files_dict['getcouple_free2.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getphasey.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getbetay.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getampbetay.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getCOx.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getDy.out'].add_filename_to_getllm_header("chrommode")

    if twiss_d.has_no_input_files():
        print >> sys.stderr, "No parsed input files"
        sys.exit(1)

    return twiss_d, files_dict
# END _analyse_src_files ----------------------------------------------------------------------------

def _check_bpm_compatibility(twiss_d, mad_twiss):
    '''
    Checks the monitor compatibility between data and model. If a monitor will not be found in the
    model, a message will be print to sys.stderr.
    '''
    all_twiss_files = twiss_d.non_zero_dpp_x + twiss_d.zero_dpp_x + twiss_d.non_zero_dpp_y + twiss_d.zero_dpp_y
    for twiss_file in all_twiss_files:
        for bpm_name in twiss_file.NAME:
            try:
                mad_twiss.NAME[mad_twiss.indx[bpm_name]]
            except KeyError:
                try:
                    mad_twiss.NAME[mad_twiss.indx[str.upper(bpm_name)]]
                except KeyError:
                    print >> sys.stderr, 'Monitor ' + bpm_name + ' cannot be found in the model!'


def _calculate_orbit(getllm_d, twiss_d, tune_d, mad_twiss, files_dict):
    '''
    Calculates orbit and fills the following TfsFiles:
     - getCOx.out
     - getCOy.out
     - getCOx_dpp_' + str(k + 1) + '.out
     - getCOy_dpp_' + str(k + 1) + '.out

    :param _GetllmData getllm_d: accel is used(In-param, values will only be read)
    :param _TwissData twiss_d: Holds twiss instances of the src files. (In-param, values will only be read)
    :param _TuneData tune_d: Holds tunes and phase advances (In-param, values will only be read)

    :returns: (list, list, dict)
     - an list of dictionairies from horizontal computations
     - an list of dictionairies from vertical computations
     - the same dict as param files_dict to indicate that dict will be extended here.
    '''
    print 'Calculating orbit'
    list_of_co_x = []
    if twiss_d.has_zero_dpp_x():
        [cox, bpms] = algorithms.helper.calculate_orbit(mad_twiss, twiss_d.zero_dpp_x)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOx.out']
        tfs_file.add_string_descriptor("TABLE", 'ORBIT')
        tfs_file.add_string_descriptor("TYPE", 'ORBIT')
        tfs_file.add_string_descriptor("SEQUENCE", getllm_d.accel)
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "X", "STDX", "XMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), cox[bn1][0], cox[bn1][1], mad_twiss.X[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

        list_of_co_x.append(cox)
    list_of_co_y = []
    if twiss_d.has_zero_dpp_y():
        [coy, bpms] = algorithms.helper.calculate_orbit(mad_twiss, twiss_d.zero_dpp_y)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOy.out']
        tfs_file.add_string_descriptor("TABLE", 'ORBIT')
        tfs_file.add_string_descriptor("TYPE", 'ORBIT')
        tfs_file.add_string_descriptor("SEQUENCE", getllm_d.accel)
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "Y", "STDY", "YMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), coy[bn1][0], coy[bn1][1], mad_twiss.Y[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

        list_of_co_y.append(coy)
    #-------- Orbit for non-zero DPP
    if twiss_d.has_non_zero_dpp_x():
        k = 0
        for twiss_file in twiss_d.non_zero_dpp_x:
            list_with_single_twiss = []
            list_with_single_twiss.append(twiss_file)
            filename = 'getCOx_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.GetllmTfsFile(filename)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(twiss_file.filename)
            tfs_file.add_float_descriptor("DPP", float(twiss_file.DPP))
            tfs_file.add_float_descriptor("Q1", tune_d.q1)
            tfs_file.add_float_descriptor("Q2", tune_d.q2)
            [codpp, bpms] = algorithms.helper.calculate_orbit(mad_twiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "X", "STDX", "XMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), codpp[bn1][0], codpp[bn1][1], mad_twiss.X[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)

            list_of_co_x.append(codpp)
            k += 1

    if twiss_d.has_non_zero_dpp_y():
        k = 0
        for twiss_file in twiss_d.non_zero_dpp_y:
            list_with_single_twiss = []
            list_with_single_twiss.append(twiss_file)
            filename = 'getCOy_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.GetllmTfsFile(filename)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(twiss_file.filename)
            tfs_file.add_float_descriptor("DPP", float(twiss_file.DPP))
            tfs_file.add_float_descriptor("Q1", tune_d.q1)
            tfs_file.add_float_descriptor("Q2", tune_d.q2)
            [codpp, bpms] = algorithms.helper.calculate_orbit(mad_twiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "Y", "STDY", "YMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                #TODO: why twiss_d.zero_dpp_y.. above used twiss_d.non_zero_dpp_y(vimaier)
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), codpp[bn1][0], codpp[bn1][1], mad_twiss.Y[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)

            list_of_co_y.append(codpp)
            k += 1

    return list_of_co_x, list_of_co_y, files_dict
# END _calculate_orbit ------------------------------------------------------------------------------


def _phase_and_beta_for_non_zero_dpp(getllm_d, twiss_d, tune_d, phase_d, bpm_dictionary, mad_twiss, files_dict, pseudo_list_x, pseudo_list_y, use_only_three_bpms_for_beta_from_phase):
    '''
    Fills the following TfsFiles:
     - getphasex_dpp_' + str(k + 1) + '.out
     - getphasey_dpp_' + str(k + 1) + '.out
     - getbetax_dpp_' + str(k + 1) + '.out
     - getbetay_dpp_' + str(k + 1) + '.out

    :param _GetllmData getllm_d: lhc_phase, accel, beam_direction and num_beams_for_coupling are used (In-param, values will only be read)
    :param _TwissData twiss_d: Holds twiss instances of the src files. (In-param, values will only be read)
    :param _TuneData tune_d: Holds tunes and phase advances. q1 and q2 will be set if accel 'SPS' or 'RHIC' (In/Out-param, values will be read and set)

    :returns: _TuneData -- the same instance as param tune_d to indicate that tunes will be set.
    '''
    print "Calculating phase and Beta for non-zero DPP"
    if DEBUG:
        print "lenght of hor-files(linx) with non zero dpp: " + str(len(twiss_d.non_zero_dpp_x))
        print "lenght of ver-files(liny) with non zero dpp: " + str(len(twiss_d.non_zero_dpp_y))

    if twiss_d.has_non_zero_dpp_x():
        plane = 'H'
        k = 0
        phasexlist = []
        for twiss_file in twiss_d.non_zero_dpp_x:
            dpop = float(twiss_file.DPP)
            list_with_single_twiss = []
            list_with_single_twiss.append(twiss_file)
            dpp_twiss = algorithms.helper.construct_off_momentum_model(mad_twiss, dpop, bpm_dictionary)
            [phasex, q1_dpp, MUX, bpms] = algorithms.phase.get_phases(getllm_d, dpp_twiss, list_with_single_twiss, tune_d.q1, plane)
            phasex['DPP'] = dpop
            phasexlist.append(phasex)
            # 'getphasex_dpp_'+str(k+1)+'.out' was inteded for debugging (vimaier)
            if DEBUG:
                filename = 'getphasex_dpp_' + str(k + 1) + '.out'
                files_dict[filename] = utils.tfs_file.GetllmTfsFile(filename)
                tfs_file = files_dict[filename]
                tfs_file.add_filename_to_getllm_header(twiss_file.filename)
                tfs_file.add_float_descriptor("DPP", dpop)
                tfs_file.add_float_descriptor("Q1", tune_d.q1)
                tfs_file.add_float_descriptor("Q2", tune_d.q2)
                tfs_file.add_float_descriptor("Q1DPP", q1_dpp)
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASE", "STDPH", "PHXMDL", "MUXMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                length_bpms = len(bpms)
                for i in range(0, length_bpms):
                    bn1 = str.upper(bpms[i][1])
                    bns1 = bpms[i][0]
                    index = (i + 1) % length_bpms
                    bn2 = str.upper(bpms[index][1])
                    bns2 = bpms[index][0]
                    try:
                        phmdl = phase_d.ph_x[bn1][4]
                    except KeyError:
                        phmdl = 0.0
                    #phmdl=mad_twiss.MUX[mad_twiss.indx[bn2]]-mad_twiss.MUX[mad_twiss.indx[bn1]]
                    list_row_entries = ['"' + bn1 + '" ', '"' + bn2 + '"', bns1, bns2, 1, phasex[bn1][0], phasex[bn1][1], phmdl, mad_twiss.MUX[mad_twiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)

            betax = {}
            alfax = {}
            [betax, rmsbbx, alfax, bpms] = algorithms.beta.beta_from_phase(mad_twiss, list_with_single_twiss, phasex, plane, use_only_three_bpms_for_beta_from_phase)
            betax['DPP'] = dpop
            #betaxa = {}
            #[betaxa, rmsbbx, bpms, invJx] = algorithms.beta.beta_from_amplitude(mad_twiss, list_with_single_twiss, plane)
            #betaxa['DPP'] = dpop
            filename = 'getbetax_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.GetllmTfsFile(filename)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(twiss_file.filename)
            tfs_file.add_float_descriptor("DPP", dpop)
            tfs_file.add_float_descriptor("Q1", tune_d.q1)
            tfs_file.add_float_descriptor("Q2", tune_d.q2)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for bpm in bpms:
                bn1 = str.upper(bpm[1])
                bns1 = bpm[0]
                list_row_entries = ['"' + bn1 + '" ', bns1, len(twiss_d.zero_dpp_x), betax[bn1][0], betax[bn1][1], betax[bn1][2], alfax[bn1][0], alfax[bn1][1], alfax[bn1][2], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.ALFX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)

            k += 1

    if twiss_d.has_non_zero_dpp_y():
        plane = 'V'
        k = 0
        phaseylist = []
        for twiss_file in twiss_d.non_zero_dpp_y:
            dpop = float(twiss_file.DPP)
            list_with_single_twiss = []
            list_with_single_twiss.append(twiss_file)
            dpp_twiss = algorithms.helper.construct_off_momentum_model(mad_twiss, dpop, bpm_dictionary)
            [phasey, q2_dpp, MUY, bpms] = algorithms.phase.get_phases(getllm_d, dpp_twiss, list_with_single_twiss, tune_d.q2, plane)
            phasey['DPP'] = dpop
            phaseylist.append(phasey)
            # 'getphasex_dpp_'+str(k+1)+'.out' was inteded for debugging (vimaier)
            if DEBUG:
                filename = 'getphasey_dpp_' + str(k + 1) + '.out'
                files_dict[filename] = utils.tfs_file.GetllmTfsFile(filename)
                tfs_file = files_dict[filename]
                tfs_file.add_filename_to_getllm_header(twiss_file.filename)
                tfs_file.add_float_descriptor("Q1", tune_d.q1)
                tfs_file.add_float_descriptor("Q2", tune_d.q2)
                tfs_file.add_float_descriptor("Q2DPP", q2_dpp)
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASE", "STDPH", "PHYMDL", "MUYMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                length_bpms = len(bpms)
                for i in range(0, length_bpms):
                    bn1 = str.upper(bpms[i][1])
                    bns1 = bpms[i][0]
                    index = (i + 1) % length_bpms
                    bn2 = str.upper(bpms[index][1])
                    bns2 = bpms[index][0]
                    try:
                        phmdl = phase_d.ph_y[bn1][4]
                    except KeyError:
                        phmdl = 0.0
                    #phmdl=mad_twiss.MUY[mad_twiss.indx[bn2]]-mad_twiss.MUY[mad_twiss.indx[bn1]]
                    list_row_entries = ['"' + bn1 + '" ', '"' + bn2 + '" ', bns1, bns2, 1, phasey[bn1][0], phasey[bn1][1], phmdl, mad_twiss.MUY[mad_twiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)

            betay = {}
            alfay = {}
            [betay, rmsbby, alfay, bpms] = algorithms.beta.beta_from_phase(dpp_twiss, list_with_single_twiss, phasey, plane, use_only_three_bpms_for_beta_from_phase)
            #betay['DPP'] = dpop
            #betaya = {}
            #[betaya, rmsbby, bpms, invJy] = algorithms.beta.beta_from_amplitude(dpp_twiss, list_with_single_twiss, plane)
            #betaya['DPP'] = dpop
            filename = 'getbetay_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.GetllmTfsFile(filename)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(twiss_file.filename)
            tfs_file.add_float_descriptor("DPP", dpop)
            tfs_file.add_float_descriptor("Q1", tune_d.q1)
            tfs_file.add_float_descriptor("Q2", tune_d.q2)
            #TODO: check if it should be Y instead of X in the column names since it is getbetaY_dpp...out (vimaier)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '" ', bns1, len(twiss_d.zero_dpp_y), betay[bn1][0], betay[bn1][1], betay[bn1][2], alfay[bn1][0], alfay[bn1][1], alfay[bn1][2], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.ALFX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)

            k += 1

    # TODO: Branch useless except of Q1 and Q2 calculation? Results will neither be saved in a file nor returned(vimaier)
    if twiss_d.has_non_zero_dpp_y() and twiss_d.has_non_zero_dpp_x():
        if len(twiss_d.non_zero_dpp_x) != len(twiss_d.non_zero_dpp_y):
            raise ValueError("list of dppx is not equal list of dppy")
        for j in range(len(twiss_d.non_zero_dpp_x)):
            dpop = float(twiss_d.non_zero_dpp_x[j].DPP)
            list_with_single_twiss_x = []
            list_with_single_twiss_y = []
            list_with_single_twiss_x.append(twiss_d.non_zero_dpp_x[j])
            list_with_single_twiss_y.append(twiss_d.non_zero_dpp_y[j])
            ### coupling
            try:
                mad_twiss.Cmatrix()
            except:
                traceback.print_exc()

            if getllm_d.accel == "SPS" or "RHIC" in getllm_d.accel:
                plane = 'H'
                [phasexp, tune_d.q1, MUX, bpmsx] = algorithms.phase.get_phases(getllm_d, mad_twiss, pseudo_list_x, None, plane)
                plane = 'V'
                [phaseyp, tune_d.q2, MUY, bpmsy] = algorithms.phase.get_phases(getllm_d, mad_twiss, pseudo_list_y, None, plane)
                [fwqw, bpms] = algorithms.coupling.GetCoupling2(mad_twiss, pseudo_list_x, pseudo_list_y, tune_d.q1, tune_d.q2, phasexp, phaseyp, getllm_d.beam_direction, getllm_d.accel, getllm_d.outputpath)
            elif getllm_d.num_beams_for_coupling == 1:
                [fwqw, bpms] = algorithms.coupling.GetCoupling1(mad_twiss, list_with_single_twiss_x, list_with_single_twiss_y, tune_d.q1, tune_d.q2, getllm_d.outputpath)
            elif getllm_d.num_beams_for_coupling == 2:
                [fwqw, bpms] = algorithms.coupling.GetCoupling2(mad_twiss, list_with_single_twiss_x, list_with_single_twiss_y, tune_d.q1, tune_d.q2, phasexlist[j], phaseylist[j], getllm_d.beam_direction, getllm_d.accel, getllm_d.outputpath)
                if getllm_d.with_ac_calc:
                    [fwqw, bpms] = algorithms.coupling.getFreeCoupling(tune_d.q1f, tune_d.q2f, tune_d.q1, tune_d.q2, fwqw, mad_twiss, bpms)
            else:
                raise ValueError('Number of monitors for coupling analysis (option -n) should be 1 or 2.')
            fwqw['DPP'] = dpop
# END _phase_and_beta_for_non_zero_dpp --------------------------------------------------------------

def _calculate_getsextupoles(twiss_d, phase_d, mad_twiss, files_dict, q1f):
    '''
    Fills the following TfsFiles:
     - getsex3000.out

    :returns: dict string --> GetllmTfsFile -- The same instace of files_dict to indicate that the dict was extended.
    '''
    print "Calculating getsextupoles"
    # For getsex1200.out andgetsex2100.out take a look at older revisions. (vimaier)

    htot, afactor, pfactor = algorithms.helper.Getsextupole(mad_twiss, twiss_d.zero_dpp_x, phase_d.ph_x, q1f, 3, 0)

    tfs_file = files_dict["getsex3000.out"]
    tfs_file.add_float_descriptor("f2h_factor", afactor)
    tfs_file.add_float_descriptor("p_f2h_factor", pfactor)
    tfs_file.add_column_names(["NAME", "S", "AMP_20", "AMP_20std", "PHASE_20", "PHASE_20std", "f3000", "f3000std", "phase_f_3000", "phase_f_3000std", "h3000", "h3000_std", "phase_h_3000", "phase_h_3000_std"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    for bpm_key in htot:
        li = htot[bpm_key]
        list_row_entries = [li[0], li[1], li[2], li[3], li[4], li[5], li[6], li[7], li[8], li[9], li[10], li[11], li[12], li[13]]
        tfs_file.add_table_row(list_row_entries)

    return files_dict
# END _calculate_getsextupoles ----------------------------------------------------------------------


def _calculate_kick(getllm_d, twiss_d, tune_d, phase_d, beta_d, mad_twiss, mad_ac, files_dict, bbthreshold, errthreshold):
    '''
    Fills the following TfsFiles:
     - getkick.out
     - getkickac.out

    :returns: dict string --> GetllmTfsFile -- The same instace of files_dict to indicate that the dict was extended
    '''
    print "Calculating kick"
    files = [twiss_d.zero_dpp_x + twiss_d.non_zero_dpp_x, twiss_d.zero_dpp_y + twiss_d.non_zero_dpp_y]

    meansqrt_2jx = {}
    meansqrt_2jy = {}
    mean_2jx = {}
    mean_2jy = {}
    bpmrejx ={}
    bpmrejy ={}

    try:
        [meansqrt_2jx, meansqrt_2jy, mean_2jx, mean_2jy, tunes, dpp , bpmrejx, bpmrejy] = algorithms.helper.getkick(files, mad_twiss, beta_d, bbthreshold, errthreshold)
    except IndexError:# occurs if either no x or no y files exist
        return files_dict

    #mean_2j = mean{2J} and meansqrt_2j=mean{sqrt(2J)}

    tfs_file_model = files_dict['getkick.out']
    tfs_file_model.add_comment("Calculates the kick from the model beta function")
    column_names_list = ["DPP", "QX", "QXRMS", "QY", "QYRMS", "NATQX", "NATQXRMS", "NATQY", "NATQYRMS", "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX", "2JXSTD", "2JY", "2JYSTD"]
    column_types_list = ["%le", "%le", "%le", "%le", "%le",     "%le",      "%le",    "%le",      "%le", "%le",      "%le",        "%le",       "%le",    "%le",   "%le",  "%le",    "%le"]
    tfs_file_model.add_column_names(column_names_list)
    tfs_file_model.add_column_datatypes(column_types_list)
    for i in range(0, len(dpp)):
        list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i], tunes[6][i], tunes[7][i], meansqrt_2jx['model'][i][0], meansqrt_2jx['model'][i][1], meansqrt_2jy['model'][i][0], meansqrt_2jy['model'][i][1], (meansqrt_2jx['model'][i][0]**2), (2*meansqrt_2jx['model'][i][0]*meansqrt_2jx['model'][i][1]), (meansqrt_2jy['model'][i][0]**2), (2*meansqrt_2jy['model'][i][0]*meansqrt_2jy['model'][i][1])]
        tfs_file_model.add_table_row(list_row_entries)

    tfs_file_phase = files_dict['getkickphase.out']
    tfs_file_phase.add_float_descriptor("Threshold_for_abs(beta_d-beta_m)/beta_m", bbthreshold)
    tfs_file_phase.add_float_descriptor("Threshold_for_uncert(beta_d)/beta_d", errthreshold)
    tfs_file_phase.add_float_descriptor("X_BPMs_Rejected", bpmrejx['phase'][i])
    tfs_file_phase.add_float_descriptor("Y_BPMs_Rejected", bpmrejy['phase'][i])
    tfs_file_phase.add_column_names(column_names_list)
    tfs_file_phase.add_column_datatypes(column_types_list)
    for i in range(0, len(dpp)):
        list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i], tunes[6][i], tunes[7][i], meansqrt_2jx['phase'][i][0], meansqrt_2jx['phase'][i][1], meansqrt_2jy['phase'][i][0], meansqrt_2jy['phase'][i][1], (meansqrt_2jx['model'][i][0]**2), (2*meansqrt_2jx['model'][i][0]*meansqrt_2jx['model'][i][1]), (meansqrt_2jy['model'][i][0]**2), (2*meansqrt_2jy['model'][i][0]*meansqrt_2jy['model'][i][1])]
        tfs_file_phase.add_table_row(list_row_entries)

    if getllm_d.with_ac_calc:
        tfs_file = files_dict['getkickac.out']
        tfs_file.add_float_descriptor("RescalingFactor_for_X", beta_d.x_ratio_f)
        tfs_file.add_float_descriptor("RescalingFactor_for_Y", beta_d.y_ratio_f)
        tfs_file.add_column_names(column_names_list+["sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES", "2JXSTDRES", "2JYRES", "2JYSTDRES"])
        tfs_file.add_column_datatypes(column_types_list+["%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        [inv_jx, inv_jy, tunes, dpp] = algorithms.compensate_ac_effect.getkickac(mad_ac, files, tune_d.q1, tune_d.q2, tune_d.q1f, tune_d.q2f, phase_d.acphasex_ac2bpmac, phase_d.acphasey_ac2bpmac, getllm_d.beam_direction, getllm_d.lhc_phase)
        for i in range(0, len(dpp)):
            #TODO: in table will be the ratio without f(beta_d.x_ratio) used but rescaling factor is f version(beta_d.x_ratio_f). Check it (vimaier)
            list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i], tunes[6][i], tunes[7][i], inv_jx[i][0], inv_jx[i][1], inv_jy[i][0], inv_jy[i][1], (inv_jx[i][0] ** 2), (2 * inv_jx[i][0] * inv_jx[i][1]), (inv_jy[i][0] ** 2), (2 * inv_jy[i][0] * inv_jy[i][1]), (inv_jx[i][0] / math.sqrt(beta_d.x_ratio)), (inv_jx[i][1] / math.sqrt(beta_d.x_ratio)), (inv_jy[i][0] / math.sqrt(beta_d.y_ratio)), (inv_jy[i][1] / math.sqrt(beta_d.y_ratio)), (inv_jx[i][0] ** 2 / beta_d.x_ratio), (2 * inv_jx[i][0] * inv_jx[i][1] / beta_d.x_ratio), (inv_jy[i][0] ** 2 / beta_d.y_ratio), (2 * inv_jy[i][0] * inv_jy[i][1] / beta_d.y_ratio)]
            tfs_file.add_table_row(list_row_entries)

    return files_dict
# END _calculate_kick -------------------------------------------------------------------------------


#===================================================================================================
# helper classes for data structures
#===================================================================================================
class _GetllmData(object):
    ''' Holds some data from parameters of main function. '''

    def __init__(self):
        '''Constructor'''
        self.outputpath = ""
        self.list_of_input_files = []

        self.accel = ""
        self.beam_direction = 0
        self.lhc_phase = ""
        self.bpm_unit = ""
        self.cut_for_closed_orbit = 0
        self.num_beams_for_coupling = 0

        self.with_ac_calc = False

    def set_outputpath(self, outputpath):
        ''' Sets the outputpath and creates directories if they not exist.

        :param string outputpath: Path to output dir. If dir(s) to output do(es) not exist, it/they will be created.
        '''
        Utilities.iotools.create_dirs(outputpath)
        self.outputpath = outputpath

    def set_bpmu_and_cut_for_closed_orbit(self, cut_co, bpm_unit):
        ''' Calculates and sets the cut and bpm unit.
        :param int cut_co: Cut in um(micrometer).
        :param string bpm_unit: Indicates used unit. um, mm, cm or m
        '''
        self.bpm_unit = bpm_unit

        if bpm_unit == 'um':
            self.cut_for_closed_orbit = cut_co
        elif bpm_unit == 'mm':
            self.cut_for_closed_orbit = cut_co / 1.0e3
        elif bpm_unit == 'cm':
            self.cut_for_closed_orbit = cut_co / 1.0e4
        elif bpm_unit == 'm':
            self.cut_for_closed_orbit = cut_co / 1.0e6
        else:
            print >> sys.stderr, "Wrong BPM unit:", bpm_unit


class _TwissData(object):
    ''' Holds twiss instances of all src files. '''
    def __init__(self):
        '''Constructor'''
        self.zero_dpp_x = [] # List of src files which have dpp==0.0
        self.non_zero_dpp_x = [] # List of src files which have dpp!=0.0
        self.zero_dpp_y = [] # List of src files which have dpp==0.0
        self.non_zero_dpp_y = [] # List of src files which have dpp!=0.0

    def has_zero_dpp_x(self):
        ''' Returns True if _linx file(s) exist(s) with dpp==0 '''
        return 0 != len(self.zero_dpp_x)

    def has_non_zero_dpp_x(self):
        ''' Returns True if _linx file(s) exist(s) with dpp!=0 '''
        return 0 != len(self.non_zero_dpp_x)

    def has_zero_dpp_y(self):
        ''' Returns True if _liny file(s) exist(s) with dpp==0 '''
        return 0 != len(self.zero_dpp_y)

    def has_non_zero_dpp_y(self):
        ''' Returns True if _liny file(s) exist(s) with dpp!=0 '''
        return 0 != len(self.non_zero_dpp_y)

    def has_no_input_files(self):
        return not self.has_zero_dpp_x() and not self.has_zero_dpp_y() and not self.has_non_zero_dpp_x() and not self.has_non_zero_dpp_y()


class _TuneData(object):
    ''' Used as data structure to hold tunes and phase advances. '''
    def __init__(self):
        '''Constructor'''
        self.q1 = 0.0 # Driven horizontal tune
        self.q2 = 0.0 # Driven vertical tune
        self.mux = 0.0 # Driven horizontal phase advance
        self.muy = 0.0 # Driven vertical phase advance

        # Free is from analytic equation
        self.q1f = 0.0 # Free horizontal tune
        self.q2f = 0.0 # Free vertical tune
        self.muxf = 0.0 # Free horizontal phase advance
        self.muyf = 0.0 # Free vertical phase advance

        # Free2 is using the effective model
        self.muxf2 = 0.0 # Free2 horizontal phase advance
        self.muyf2 = 0.0 # Free2 vertical phase advance

        self.delta1 = None # Used later to calculate free Q1. Only if with ac calculation.
        self.delta2 = None # Used later to calculate free Q2. Only if with ac calculation.

    def initialize_tunes(self, with_ac_calc, mad_twiss, mad_ac, twiss_d):
        ''' Calculates and sets the initial tunes. '''
        if with_ac_calc:
            # Get fractional part: frac(62.23) = 0.23; 62.23 % 1 ==> 0.23 (vimaier)
            self.q1f = abs(mad_twiss.Q1) % 1 #-- Free Q1 (tempolarlly, overwritten later)
            self.q2f = abs(mad_twiss.Q2) % 1 #-- Free Q2 (tempolarlly, overwritten later)
            self.q1 = abs(mad_ac.Q1) % 1 #-- Drive Q1 (tempolarlly, overwritten later)
            self.q2 = abs(mad_ac.Q2) % 1 #-- Drive Q2 (tempolarlly, overwritten later)
            self.delta1 = self.q1-self.q1f #-- Used later to calculate free Q1
            self.delta2 = self.q2-self.q2f #-- Used later to calculate free Q2
        else:
            try:
                self.q1f = twiss_d.zero_dpp_x[0].Q1
                self.q2f = twiss_d.zero_dpp_y[0].Q2
            except IndexError:
                pass

#===================================================================================================
# main invocation
#===================================================================================================
def _start():
    '''
    Starter function to avoid polluting global namespace with variables options,args.
    Before the following code was after 'if __name__=="__main__":'
    '''
    options = _parse_args()
    main(outputpath=options.output,
         dict_file=options.dict,
         files_to_analyse=options.files,
         model_filename=options.Twiss,
         accel=options.ACCEL,
         lhcphase=options.lhcphase,
         BPMU=options.BPMUNIT,
         COcut=float(options.COcut),
         NBcpl=int(options.NBcpl),
         TBTana=options.TBTana,
         higher_order=options.higher,
         bbthreshold=options.bbthreshold,
         errthreshold=options.errthreshold)

if __name__ == "__main__":
    _start()


