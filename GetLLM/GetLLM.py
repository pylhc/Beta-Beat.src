'''
Created on 11/09/09

@author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)

@version: 2.38b

TODO: more detailed description(vimaier)

Python script to obtain Linear Lattice functions and More -> GetLLM



Usage1 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./
Usage2 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./




Some rules for variable name:
- Dictionary is used to contain the output of function
- Variable containing 'm' is a value directly obtained from measurment data
- Variable containing 'mdl' is a value related to model

Change history:
V1.0, 11/Feb/2008 by Masa. Aiba
V1.1, 18/Feb/2008:
- Debugging, add model phase and tunes to output
- add function to obtain DY
- add chromatic parameter (phase for non zero DPP)
V1.2, 22/Feb/2008:
- test version for beta with all BPM
V1.3, 29/Feb/2008:
- beta from phases is improved, averaging beta1, 2 and 3
V1.31, 12/Mar/2008:
- debugged on alpha3
V1.4, 12/Mar/2008:
- modify output to fit latest TSF format and to meet requests from Rogelio
- fix buggs in r.m.s. beta-beat and added to the output of getbetax/y.out
V1.5 Rogelio, 13 March 2008:
- Update to option parser, include BPMdictionary to filter BPMs not in Model
V1.51, 13/Mar/2008:
- Modify output to fit latest TSF format again. Add STD to beta.
V1.6, 15/Jul/2008:
- Add the integer part of tunes - assuming that the phase advance is always
  less than 1.0.
V1.71 27/Jul/2008:
- Add GetCO. Filter in dispersion calculation to exclude bad bpms.
V1.8, 13/Aug/2008 Ref. note by A. Franchi, R. T. Garcia, G. Vanbavinckhove:
- Add GetCoupling
- "Computation of the Coupling Resonance Driving term f1001 and the coupling
  coefficient C from turn-by-turn single-BPM data", 28/May/2008
- The GetCoupling.py is initiated by Glenn V. and imported into GetLLM /
  finalized by Masa. Aiba
- Some bugs are fixed - you can run even without liny file and find the
  results for horizontal plane only.
V1.81, 1/Sep/2008:
- For an accelerator in which the beam goes opposite direction to the model
  as in LHCB2, the beam direction parameter (bd) is added to GetPhases.
- Bug in the phi13 for the last monitor and the last but one is fixed.
V1.9, 21/Oct/2008:
- Add the beta from spectrum height.
V1.91, 08/Dec/2008:
- Add option - SUSSIX or SVD for input file
V1.99, 13/Feb/09:
- Add DPX! The MEMORIAL version for Version 1.** since all the linear lattice
  parameters become available!
- Add the option for the Harmonic analysis.
- Refine coding, add several lines of comment
V2.0, 17/Feb/2009:
- Major version up to output "More", that is, other than the linear lattice
  parameters.
- Add off-momentum lattice (dbeta/beta)/(dp/p)
- Add SPS coupling with "Pseudo-double plane BPM"-it need at this moment a
  list containing
- pre-paired H-V monitor. This should be replaced by a clever algorithm to
  find good pairs automatically.
V2.01, 10/Mar/2009:
- Fix bug on SPS double plane BPM monitor, in which missing BPM could cause
  an error.
- Modify BetaFromAmplitude to output invariant J (Rogelio / finalised by MA)
V2.02, 10/Mar/2009:
- Fix bug in getcoupling, bad input for phasex and phasey
V2.10, 13/Mar/2009, Glenn Vanbavinckhove:
- Added function for finding sextupole lines (amp and phases) + chiterms amp
V2.11. 26/Mar/2009
- Fix bug in Getphase (model phase advance from the last to first monitor).
V2.12, 28/03/2009
- Following bugs fixed : avoid negative square, change option -h to -l
- Add -r option as in correct.py
- Change the way to import BPM pair file for SPS. (import -> execfile)
V2.13, 06/Apl/2009
- Fix bug in weight function to accept negative dpp
V2.14, 08/Apl/2009
- Fix bug in Normalized dispersion to treat COcut correctly.
V2.15:
- Enable coupling in RHIC as in SPS
V2.16, 28/May/2009:
- Add STDBET for the beta from amplitude.
- Add option for off momentum beta-beating to choose algorithm, that is, beta
  from phase or amp
- Add a routine to detect wrong data having two lines in linx/y file with the
  same BPM name.
- Add a routine to avoid zero division due to exactly n*pi phase advance in
  beta from phase (see the last part of GetPhases).
V2.21, 23/June/2009:
- Add STDBET Model for off momentum beta beat phase.
V2.25:
- Fixed bd flag (must be -1 for beam2)
V2.25:
- Adding VERSION variable to be always output and modified in subsequent
  versions, do not forget!!!
V2.26:
- Adding F2000 (two different methods linear and non-linear)
V2.27:
- Adding new method for IP calculation
V2.28:
- Changing the rejection of bad BPM for the coupling phase - averaging the
  phase over the sets of data first , then cut if the q1 and q2 are very
  different.
- 24/Feb/2010 the change is not yet checked.
- Hyphens in the @ field of tfs files is not allowed:  The previous label
  "RMS-beta-beat" has been  moved to "RMSbetabeat"
V2.29 5/Mar/2010
- Change the default value for COcut from 1000 to 4000 as it was too small
V2.30 13/Sept/2010:
- Updating for AC-Dipole, implementing chromatic coupling (not done for RHIC
  or SPS)
V2.31 8/November/2010:
- Update for AC-Dipole: gives free beta,phase,coupling(global factor)
V2.32 15/January/2010:
- Taking models
V2.33 7/February/2011:
- implementing coupling correction for AC-dipole.
V2.34 7/04/2011:
- Updating to deal with chromatic twiss files.
V2.35 9/06/2011:
- Functions to cancel the AC dipole effect for beta, phase, and total phase,
  based on equations, are added.
- A function to calculate beta* from the phase advance between Q1s is added.
- Phase shift by tune is compensated for the LHC experiment data.
V2.36 30/09/2011:
- Rescaling algorithm for BetaFromAmp and action (by Andy) is implemented.
- 2nd function to calculate IP parameters from Q1s is added.
- The compensation of the phase shift by tune for LHC exp data is modified.
- A function to calculate action of the AC dipole excitation is added.
V2.38 08/03/2012:
- added main() function
- using raise ValueError() instead of sys.exit() some places
13/09/2012:
- merged in patch 2.35-2.37
V2.38b 03/dec/2012, tbach:
- reformatted comments
- changed all parts of code, where the program exits inside an exception to
  display the exception, because the messages are not always helpful
- removed ";" from all over the code (still some left)
- 207 errors, 983 warning, 574 infos left from static code analysis...
 - x.xxx, vimaier  16th Apr 2013:
    deleted functions function and GetCoupling2
    Changed GetCoupling2b to GetCoupling2
    Set some TODOs
    Changed wolin* and acswitch (* := x|y|x2|y2)
    Defined variables Q1,Q2,Q1f,Q2f,MUX,MUY,muxf,muyf with standard values
      Saves a lot of try/excepts while writing files
    Reformatted a lot of code
    Introduced tfs_file for all output files in main() --> Formatted output files
    Extracted functions modelIntersect and intersect to Utilities.bpm
    Major refactoring:
        Extracted helper functions into GetLLM.algorithms.helper.py
        extracted main into smaller functions

'''
import sys
if "/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/" not in sys.path: # add internal path for python scripts to current environment (tbach)
    sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')
if "/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" not in sys.path: # added for Utilities.bpm (vimaier)
    sys.path.append('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/')
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import string
import math

import numpy as np

import metaclass
import traceback
import utils.tfs_file
import algorithms.helper
# tentative solution for SPS pseudo double plane BPM
# from SPSBPMpair import *




####
#######
#########
VERSION='V2.38b PRO'
#########
#######
####
DEBUG = True 
''' For now only used to decide if getphase(x|y)_dpp_' + str(k + 1) + '.out will be created or not. '''


#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
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

    (options, args) = parser.parse_args()
    return options,args

#===================================================================================================
# main()-function
#===================================================================================================
def main(outputpath,files_to_analyse,model_filename,dict_file="0",accel="LHCB1",lhcphase="0",BPMU="um",COcut=4000,NBcpl=2,TBTana="SUSSIX",higher_order=1):
    '''
     GetLLM main function.
     
     :Parameters:
        'outputpath': string
            The output path to store results
        'files_to_analyse': string
            List of files, comma separated string.
        'dict_file': string
            Name of the script which will be executed. Should store dictionary with mappings of 
            BPM names.
        'accel': string
            Type of accelerator. LHCB1, LHCB2, LHCB4, RHIC, SPS
        'lhcphase': string "0" or "1"
            Compensate phase shifts by tunes for the LHC experiment data, off=0(default)/on=1
        'BPMU': string
            BPMunit: um, mm, cm, m (default um)
        'COcut': int
            Cut for closed orbit measurement [um]
        'NBcpl': int
            For selecting the coupling measurement method 1 bpm or 2 bpms
        'TBTana': string
            Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA
        'higher_order': int
            output higher order resonance stuff
                
        :Return: int
            0 if the function run successfully otherwise !=0. 
    '''

    print "Starting GetLLM ", VERSION

    outputpath, accel, with_ac_calc, listOfInputFiles, Suffix_x, Suffix_y, MADTwiss, MADTwiss_ac, BPMdictionary, beam_direction, MADTwissElem, COcut = intial_setup(outputpath, files_to_analyse, model_filename, dict_file, accel, BPMU, COcut, TBTana)

    #-------- create TfsFile
    files_dict = create_tfs_files(outputpath, model_filename, accel, with_ac_calc)

    ListOfZeroDPPX, ListOfZeroDPPY, with_linx_for_zero_dppx, with_liny_for_zero_dppy, ListOfNonZeroDPPX, ListOfNonZeroDPPY, FileOfNonZeroDPPX, FileOfNonZeroDPPY, with_liny_for_nonzero_dppy, with_linx_for_nonzero_dppx = analyse_src_files(accel, with_ac_calc, listOfInputFiles, Suffix_x, Suffix_y, files_dict)

    Q1 = 0.0
    Q2 = 0.0
    MUX = 0.0
    MUY = 0.0
    
    Q1f = 0.0
    Q2f = 0.0
    muxf = 0.0
    muyf = 0.0
    
    d1 = None
    d2 = None
    
    if with_ac_calc:
        # Get fractional part: frac(62.23) = 0.23; 62.23 % 1 ==> 0.23 (vimaier)
        Q1f = abs(MADTwiss.Q1) % 1 #-- Free Q1 (tempolarlly, overwritten later)
        Q2f = abs(MADTwiss.Q2) % 1 #-- Free Q2 (tempolarlly, overwritten later)
        Q1 = abs(MADTwiss_ac.Q1) % 1 #-- Drive Q1 (tempolarlly, overwritten later)
        Q2 = abs(MADTwiss_ac.Q2) % 1 #-- Drive Q2 (tempolarlly, overwritten later)
        d1 =Q1-Q1f #-- Used later to calculate free Q1
        d2 =Q2-Q2f #-- Used later to calculate free Q2
    else:
        Q1f=ListOfZeroDPPX[0].Q1
        Q2f=ListOfZeroDPPY[0].Q2
        
    #TODO: initialize variables otherwise calculate_coupling would raie an exception(vimaier)
    PseudoListX = None
    PseudoListY = None
    # Construct pseudo-double plane BPMs
    if (accel=="SPS" or "RHIC" in accel) and with_linx_for_zero_dppx and with_liny_for_zero_dppy:
        execfile(model_filename.replace("twiss.dat","BPMpair.py"))
        [PseudoListX,PseudoListY] = algorithms.helper.PseudoDoublePlaneMonitors(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, BPMdictionary)


    #-------- Check monitor compatibility between data and model
    all_twiss_files = ListOfNonZeroDPPX+ListOfZeroDPPX+ListOfNonZeroDPPY+ListOfZeroDPPY
    for twiss_file in all_twiss_files:
        for bpm_name in twiss_file.NAME:
            #TODO: maybe easier with the usage of intersect?(vimaier)
            # Check if all BPMs are in the model(vimaier)
            try:
                MADTwiss.NAME[MADTwiss.indx[bpm_name]]
            except:
                try:
                    MADTwiss.NAME[MADTwiss.indx[string.upper(bpm_name)]]
                except:
                    print 'Monitor '+bpm_name+' cannot be found in the model!'
                    #exit()

    #-------- START Phase
    Q1,  MUX, Q2, MUY, Q1f, acphasex_ac2bpmac, muxf, Q2f, muyf, muxf2, muyf2, acphasey_ac2bpmac, phasex, phasexf, phasey, phaseyf, bpmsx, bpmsy, phasexf2, phaseyf2, phasexlist, phaseylist = calculate_phase(outputpath, accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, MADTwissElem, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, MUX, MUY, Q1f, Q2f, muxf, muyf, d1, d2)

    #-------- START Total Phase
    calculate_total_phase(accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, MUX, MUY, Q1f, Q2f, muxf, muyf, acphasex_ac2bpmac, muxf2, muyf2, acphasey_ac2bpmac)

    #-------- START Beta
    bpms, betax, betaxf, betay, betayf = calculate_beta(MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, Q1f, Q2f, phasex, phasexf, phasey, phaseyf)

    #------- START beta from amplitude
    betax, betay, bpms, beta2_save, betaxalist, betayalist, betax_ratio, betay_ratio, betaxf_ratio, betayf_ratio = calculate_beta_from_amplitude(accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, Q1f, Q2f, acphasex_ac2bpmac, acphasey_ac2bpmac, bpms, betax, betaxf, betay, betayf)

    #-------- START IP
    calculate_ip(accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, Q1f, Q2f, bpm_name, acphasex_ac2bpmac, acphasey_ac2bpmac, phasex, phasexf, phasey, phaseyf, bpmsx, bpmsy, phasexf2, phaseyf2, betax, betay)

    #-------- START Orbit
    ListOfCOX, ListOfCOY = calculate_orbit(model_filename, accel, MADTwiss, files_dict, FileOfNonZeroDPPX, FileOfNonZeroDPPY, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, with_liny_for_zero_dppy, with_liny_for_nonzero_dppy, with_linx_for_zero_dppx, with_linx_for_nonzero_dppx, Q1, Q2, bpms)

    #-------- START Dispersion
    calculate_dispersion(BPMU, COcut, MADTwiss, files_dict, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, with_liny_for_zero_dppy, with_liny_for_nonzero_dppy, with_linx_for_zero_dppx, with_linx_for_nonzero_dppx, Q1, Q2, bpms, beta2_save, ListOfCOX, ListOfCOY)

    #-------- START coupling.
    Q1, Q2 = calculate_coupling(outputpath, accel, lhcphase, NBcpl, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, Q1f, Q2f, PseudoListX, PseudoListY, acphasex_ac2bpmac, acphasey_ac2bpmac, phasexlist, phaseylist, bpms)

    #-------- Phase, Beta and coupling for non-zero DPP
    Q1, Q2 = phase_and_beta_for_non_zero_dpp(outputpath, model_filename, accel, lhcphase, NBcpl, BPMdictionary, beam_direction, MADTwiss, with_ac_calc, files_dict, FileOfNonZeroDPPX, FileOfNonZeroDPPY, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, with_liny_for_nonzero_dppy, with_linx_for_nonzero_dppx, Q1, Q2, Q1f, Q2f, PseudoListX, PseudoListY, phasex, phasey, phasexlist, phaseylist, bpms, betax, betay, betaxalist, betayalist)

    if not higher_order:
        print "Not analysing higher order..."
        return 0
    
    if TBTana == "SUSSIX":
        #------ Start getsextupoles @ Glenn Vanbavinckhove
        calculate_getsextupoles(model_filename, TBTana, MADTwiss, files_dict, ListOfZeroDPPX, Q1f, phasexlist)
    
        #------ begin getchiterms @ Glenn Vanbavinckhove
        calculate_chiterms(model_filename, accel, MADTwiss, files_dict, ListOfZeroDPPX, ListOfZeroDPPY)

    #----------------------------- begin get Q,JX,delta
    calculate_kick(model_filename, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, Q1, Q2, Q1f, Q2f, acphasex_ac2bpmac, acphasey_ac2bpmac, betax_ratio, betay_ratio, betaxf_ratio, betayf_ratio)


    #-------- Write files
    for tfsfile in files_dict.itervalues():
        tfsfile.write_to_file(formatted=True)
        
# END main() ---------------------------------------------------------------------------------------

    
#===================================================================================================
# helper-functions
#===================================================================================================
def intial_setup(outputpath, files_to_analyse, model_filename, dict_file, accel, BPMU, COcut, TBTana):
    outputpath = algorithms.helper.fix_output(outputpath)
    if dict_file == "0":
        BPMdictionary = {}
    else:
        execfile(dict_file)
        BPMdictionary = dictionary # temporaryly since presently name is not BPMdictionary
    listOfInputFiles = files_to_analyse.split(",")
# Beam direction
    beam_direction = 1
    if accel == "LHCB2":
        beam_direction = -1 # THIS IS CORRECT, be careful with tune sign in SUSSIX and eigenmode order in SVD
    elif accel == "LHCB4":
        beam_direction = 1 # IS THIS CORRECT? I (rogelio) try for Simon...
        accel = "LHCB2" #This is because elements later are named B2 anyway, not B4
#-- finding base model
    try:
        MADTwiss = metaclass.twiss(model_filename, BPMdictionary) # MODEL from MAD
        print "Base model found!"
    except:
        print >> sys.stderr, "twiss file loading failed for:", model_filename
        print >> sys.stderr, traceback.format_exc()
        #-- finding the ac dipole model
        sys.exit(1)
    with_ac_calc = False
    try:
        MADTwiss_ac = metaclass.twiss(model_filename.replace(".dat", "_ac.dat"))
        with_ac_calc = True
        print "Driven Twiss file found. AC dipole effects calculated with the effective model (get***_free2.out)"
    except:
        MADTwiss_ac = MADTwiss
        print "WARN: AC dipole effects not calculated. Driven twiss file does not exsist !"
#-- Test if the AC dipole (MKQA) is in the model of LHC
    MADTwissElem = None
    if with_ac_calc:
        if 'LHC' in accel:
            if 'MKQA.6L4.' + accel[3:] in MADTwiss.NAME:
                print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
            else:
                try:
                    MADTwissElem = metaclass.twiss(model_filename.replace(".dat", "_elements.dat"))
                    print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
                except:
                    print 'WARN: AC dipoles not in the model. AC dipole effects not calculated with analytic equations !'
        else:
            print 'WARN: AC dipole effects calculated with analytic equations only for LHC for now'
    if BPMU == 'um':
        COcut = COcut
    elif BPMU == 'mm':
        COcut = COcut / 1.0e3
    elif BPMU == 'cm':
        COcut = COcut / 1.0e4
    elif BPMU == 'm':
        COcut = COcut / 1.0e6
    if TBTana == "SUSSIX":
        Suffix_x = '_linx'
        Suffix_y = '_liny'
    elif TBTana == 'SVD':
        Suffix_x = '_svdx'
        Suffix_y = '_svdy'
    elif TBTana == 'HA':
        Suffix_x = '_hax'
        Suffix_y = '_hay'
    return outputpath, accel, with_ac_calc, listOfInputFiles, Suffix_x, Suffix_y, MADTwiss, MADTwiss_ac, BPMdictionary, beam_direction, MADTwissElem, COcut
# END intial_setup ---------------------------------------------------------------------------------


def create_tfs_files(outputpath, model_filename, accel, with_ac_calc):
    '''
    Creates the most tfs files and stores it in an dictionary whereby the key represents the file
    and the value is the corresponding TfsFile.

    :Return: dict: string --> TfsFile
            A dictionary of created TfsFile objects. Keys are the filenames and values are the 
            TfsFile objects.
    '''
    # Static variable of TfsFile to save the outputfile
    utils.tfs_file.TfsFile.s_output_path = outputpath
    files_dict = {}
    files_dict['getphasex.out'] = utils.tfs_file.TfsFile('getphasex.out').add_getllm_header(VERSION, model_filename)
    files_dict['getphasey.out'] = utils.tfs_file.TfsFile('getphasey.out').add_getllm_header(VERSION, model_filename)
    files_dict['getphasetotx.out'] = utils.tfs_file.TfsFile('getphasetotx.out').add_getllm_header(VERSION, model_filename)
    files_dict['getphasetoty.out'] = utils.tfs_file.TfsFile('getphasetoty.out').add_getllm_header(VERSION, model_filename)
    if with_ac_calc:
        files_dict['getphasex_free.out'] = utils.tfs_file.TfsFile('getphasex_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasey_free.out'] = utils.tfs_file.TfsFile('getphasey_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasex_free2.out'] = utils.tfs_file.TfsFile('getphasex_free2.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasey_free2.out'] = utils.tfs_file.TfsFile('getphasey_free2.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasetotx_free.out'] = utils.tfs_file.TfsFile('getphasetotx_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasetoty_free.out'] = utils.tfs_file.TfsFile('getphasetoty_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasetotx_free2.out'] = utils.tfs_file.TfsFile('getphasetotx_free2.out').add_getllm_header(VERSION, model_filename)
        files_dict['getphasetoty_free2.out'] = utils.tfs_file.TfsFile('getphasetoty_free2.out').add_getllm_header(VERSION, model_filename)
    files_dict['getbetax.out'] = utils.tfs_file.TfsFile('getbetax.out').add_getllm_header(VERSION, model_filename)
    files_dict['getbetay.out'] = utils.tfs_file.TfsFile('getbetay.out').add_getllm_header(VERSION, model_filename)
    if with_ac_calc:
        files_dict['getbetax_free.out'] = utils.tfs_file.TfsFile('getbetax_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getbetay_free.out'] = utils.tfs_file.TfsFile('getbetay_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getbetax_free2.out'] = utils.tfs_file.TfsFile('getbetax_free2.out').add_getllm_header(VERSION, model_filename)
        files_dict['getbetay_free2.out'] = utils.tfs_file.TfsFile('getbetay_free2.out').add_getllm_header(VERSION, model_filename)
    files_dict['getampbetax.out'] = utils.tfs_file.TfsFile('getampbetax.out').add_getllm_header(VERSION, model_filename)
    files_dict['getampbetay.out'] = utils.tfs_file.TfsFile('getampbetay.out').add_getllm_header(VERSION, model_filename)
    if with_ac_calc:
        files_dict['getampbetax_free.out'] = utils.tfs_file.TfsFile('getampbetax_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getampbetay_free.out'] = utils.tfs_file.TfsFile('getampbetay_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getampbetax_free2.out'] = utils.tfs_file.TfsFile('getampbetax_free2.out').add_getllm_header(VERSION, model_filename)
        files_dict['getampbetay_free2.out'] = utils.tfs_file.TfsFile('getampbetay_free2.out').add_getllm_header(VERSION, model_filename)
    files_dict['getCOx.out'] = utils.tfs_file.TfsFile('getCOx.out').add_getllm_header(VERSION, model_filename)
    files_dict['getCOy.out'] = utils.tfs_file.TfsFile('getCOy.out').add_getllm_header(VERSION, model_filename)
    files_dict['getNDx.out'] = utils.tfs_file.TfsFile('getNDx.out').add_getllm_header(VERSION, model_filename)
    files_dict['getDx.out'] = utils.tfs_file.TfsFile('getDx.out').add_getllm_header(VERSION, model_filename)
    files_dict['getDy.out'] = utils.tfs_file.TfsFile('getDy.out').add_getllm_header(VERSION, model_filename)
    files_dict['getcouple.out'] = utils.tfs_file.TfsFile('getcouple.out').add_getllm_header(VERSION, model_filename)
    if with_ac_calc:
        files_dict['getcouple_free.out'] = utils.tfs_file.TfsFile('getcouple_free.out').add_getllm_header(VERSION, model_filename)
        files_dict['getcouple_free2.out'] = utils.tfs_file.TfsFile('getcouple_free2.out').add_getllm_header(VERSION, model_filename)
    files_dict['getcoupleterms.out'] = utils.tfs_file.TfsFile('getcoupleterms.out').add_getllm_header(VERSION, model_filename)
    if "LHC" in accel:
        files_dict['getIP.out'] = utils.tfs_file.TfsFile('getIP.out').add_getllm_header(VERSION, model_filename)
        files_dict['getIPx.out'] = utils.tfs_file.TfsFile('getIPx.out').add_getllm_header(VERSION, model_filename)
        files_dict['getIPy.out'] = utils.tfs_file.TfsFile('getIPy.out').add_getllm_header(VERSION, model_filename)
        files_dict['getIPfromphase.out'] = utils.tfs_file.TfsFile('getIPfromphase.out').add_getllm_header(VERSION, model_filename)
        if with_ac_calc:
            files_dict['getIPx_free.out'] = utils.tfs_file.TfsFile('getIPx_free.out').add_getllm_header(VERSION, model_filename)
            files_dict['getIPy_free.out'] = utils.tfs_file.TfsFile('getIPy_free.out').add_getllm_header(VERSION, model_filename)
            files_dict['getIPx_free2.out'] = utils.tfs_file.TfsFile('getIPx_free2.out').add_getllm_header(VERSION, model_filename)
            files_dict['getIPy_free2.out'] = utils.tfs_file.TfsFile('getIPy_free2.out').add_getllm_header(VERSION, model_filename)
            files_dict['getIPfromphase_free.out'] = utils.tfs_file.TfsFile('getIPfromphase_free.out').add_getllm_header(VERSION, model_filename)
            files_dict['getIPfromphase_free2.out'] = utils.tfs_file.TfsFile('getIPfromphase_free2.out').add_getllm_header(VERSION, model_filename)
    return files_dict
# END create_tfs_files -----------------------------------------------------------------------------


def analyse_src_files(accel, with_ac_calc, listOfInputFiles, Suffix_x, Suffix_y, files_dict):
    FileOfNonZeroDPPX = []
    FileOfNonZeroDPPY = []
    ListOfZeroDPPX = []
    ListOfNonZeroDPPX = []
    ListOfZeroDPPY = []
    ListOfNonZeroDPPY = []
    for file_in in listOfInputFiles:
        if file_in.endswith(".gz"):
            file_x = file_in.replace(".gz", Suffix_x + ".gz")
        else:
            file_x = file_in + Suffix_x
        twiss_file_x = metaclass.twiss(file_x)
        try:
            dppi = twiss_file_x.DPP
        except:
            dppi = 0.0
        if type(dppi) != float:
            print 'Warning: DPP may not be given as a number in ', file_x, '...trying to forcibly cast it as a number'
            try:
                dppi = float(dppi)
                print 'dppi= ', dppi
            except:
                print >> sys.stderr, 'but failing. DPP in ', file_x, ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)
        if dppi == 0.0:
            ListOfZeroDPPX.append(twiss_file_x)
            files_dict['getphasex.out'].add_filename_to_getllm_header(file_x)
            files_dict['getphasetotx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getbetax.out'].add_filename_to_getllm_header(file_x)
            files_dict['getampbetax.out'].add_filename_to_getllm_header(file_x)
            files_dict['getCOx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getNDx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getDx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getcouple.out'].add_filename_to_getllm_header(file_in)
            if "LHC" in accel:
                files_dict['getIPx.out'].add_filename_to_getllm_header(file_in)
                files_dict['getIPy.out'].add_filename_to_getllm_header(file_in)
                files_dict['getIPfromphase.out'].add_filename_to_getllm_header(file_in)
                if with_ac_calc:
                    files_dict['getIPx_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPy_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPx_free2.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPy_free2.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPfromphase_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPfromphase_free2.out'].add_filename_to_getllm_header(file_in)
            if with_ac_calc:
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
            ListOfNonZeroDPPX.append(twiss_file_x)
            FileOfNonZeroDPPX.append(file_x)
            files_dict['getNDx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getDx.out'].add_filename_to_getllm_header(file_x)
        try:
            if file_in.endswith(".gz"):
                file_y = file_in.replace(".gz", Suffix_y + ".gz")
            else:
                file_y = file_in + Suffix_y
            twiss_file_y = metaclass.twiss(file_y)
            try:
                dppi = twiss_file_y.DPP
            except:
                dppi = 0.0
            if type(dppi) != float:
                print 'Warning: DPP may not be given as a number in ', file_y, '...trying to forcibly cast it as a number'
                try:
                    dppi = float(dppi)
                    print 'dppi= ', dppi
                except:
                    print 'but failing. DPP in ', file_y, ' is something wrong. String?'
                    raise ValueError('leaving GetLLM')
            if dppi == 0.0:
                ListOfZeroDPPY.append(twiss_file_y)
                files_dict['getphasey.out'].add_filename_to_getllm_header(file_y)
                files_dict['getphasetoty.out'].add_filename_to_getllm_header(file_y)
                files_dict['getbetay.out'].add_filename_to_getllm_header(file_y)
                files_dict['getampbetay.out'].add_filename_to_getllm_header(file_y)
                files_dict['getCOy.out'].add_filename_to_getllm_header(file_y)
                files_dict['getDy.out'].add_filename_to_getllm_header(file_y)
                if with_ac_calc:
                    files_dict['getphasey_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasey_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasetoty_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasetoty_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getbetay_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getbetay_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getampbetay_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getampbetay_free2.out'].add_filename_to_getllm_header(file_y)
            else:
                ListOfNonZeroDPPY.append(twiss_file_y)
                FileOfNonZeroDPPY.append(file_y)
                files_dict['getDy.out'].add_filename_to_getllm_header(file_y)
        except:
            print 'Warning: There seems no ' + str(file_y) + ' file in the specified directory.'
    
    with_liny_for_zero_dppy = True #FLAG meaning there is _liny file for zero DPPY!
    with_liny_for_nonzero_dppy = True #FLAG meaning there is _liny file for non-zero DPPY!
    with_linx_for_zero_dppx = True
    with_linx_for_nonzero_dppx = True
    if len(ListOfZeroDPPY) == 0:
        with_liny_for_zero_dppy = False
    if len(ListOfNonZeroDPPY) == 0:
        with_liny_for_nonzero_dppy = False
    if len(ListOfZeroDPPX) == 0:
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'
        with_linx_for_zero_dppx = False
    if len(ListOfNonZeroDPPX) == 0:
        with_linx_for_nonzero_dppx = False
    if (len(ListOfNonZeroDPPX) != 0) and (len(ListOfZeroDPPX) == 0):
        ListOfZeroDPPX = ListOfNonZeroDPPX
        ListOfZeroDPPY = ListOfNonZeroDPPY
        with_linx_for_zero_dppx = True
        with_liny_for_zero_dppy = True
        with_liny_for_nonzero_dppy = False
        with_linx_for_nonzero_dppx = False
        print "Previous warning suppressed, running in chromatic mode"
        files_dict['getphasex.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getbetax.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getampbetax.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getCOx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getNDx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getDx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getcouple.out'].add_filename_to_getllm_header("chrommode")
        if with_ac_calc:
            files_dict['getcouple_free.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getcouple_free2.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getphasey.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getbetay.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getampbetay.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getCOx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getDy.out'].add_filename_to_getllm_header("chrommode")
    return ListOfZeroDPPX, ListOfZeroDPPY, with_linx_for_zero_dppx, with_liny_for_zero_dppy, ListOfNonZeroDPPX, ListOfNonZeroDPPY, FileOfNonZeroDPPX, FileOfNonZeroDPPY, with_liny_for_nonzero_dppy, with_linx_for_nonzero_dppx
# END analyse_src_files ----------------------------------------------------------------------------


def calculate_phase(outputpath, accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, MADTwissElem, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, MUX, MUY, Q1f, Q2f, muxf, muyf, d1, d2):
    #TODO: initialize variables otherwise return stmt of this function would raie an exception(vimaier)
    acphasex_ac2bpmac = None
    acphasey_ac2bpmac = None
    muxf2 = None
    muyf2 = None
    phasexf = None
    phasexf2 = None
    phaseyf = None
    phaseyf2 = None
    
    print 'Calculating phase' 
#---- Calling GetPhases first to save tunes
    #TODO:redundant code here? Better: if x: doX(); if y: doY();? Check it (vimaier)
    if with_linx_for_zero_dppx and with_liny_for_zero_dppy:
        #-- Calculate temp value of tune
        Q1temp = []
        Q2temp = []
        for twiss_f in ListOfZeroDPPX:
            Q1temp.append(np.mean(twiss_f.TUNEX))
        
        for twiss_f in ListOfZeroDPPY:
            Q2temp.append(np.mean(twiss_f.TUNEY))
        
        Q1temp = np.mean(Q1temp)
        Q2temp = np.mean(Q2temp)
        if len(ListOfZeroDPPX[0].NAME) == 0:
            print "No BPMs in linx file"
            sys.exit(1)
        if len(ListOfZeroDPPY[0].NAME) == 0:
            print "No BPMs in liny file"
            sys.exit(1)
        [phasex, Q1, MUX, bpmsx] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPX, Q1temp, 'H', outputpath, beam_direction, accel, lhcphase)
        [phasey, Q2, MUY, bpmsy] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPY, Q2temp, 'V', outputpath, beam_direction, accel, lhcphase) 
        #TODO: what is KK??(vimaier)
        print "KK end"
    elif with_linx_for_zero_dppx: #-- Calculate temp value of tune
        Q1temp = []
        for i in ListOfZeroDPPX:
            Q1temp.append(np.mean(i.TUNEX))
        
        Q1temp = np.mean(Q1temp)
        [phasex, Q1, MUX, bpmsx] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPX, Q1temp, 'H', outputpath, beam_direction, accel, lhcphase)
        print 'liny missing and output x only ...'
    elif with_liny_for_zero_dppy: #-- Calculate temp value of tune
        Q2temp = []
        for i in ListOfZeroDPPY:
            Q2temp.append(np.mean(i.TUNEY))
        
        Q2temp = np.mean(Q2temp)
        [phasey, Q2, MUY, bpmsy] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPY, Q2temp, 'V', outputpath, beam_direction, accel, lhcphase)
        print 'linx missing and output y only ...'

#---- Re-run GetPhase to fix the phase shift by Q for exp data of LHC
    if lhcphase == "1":
        if with_linx_for_zero_dppx and with_liny_for_zero_dppy:
            [phasex, Q1, MUX, bpmsx] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPX, Q1, 'H', outputpath, beam_direction, accel, lhcphase)
            [phasey, Q2, MUY, bpmsy] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPY, Q2, 'V', outputpath, beam_direction, accel, lhcphase)
        elif with_linx_for_zero_dppx:
            [phasex, Q1, MUX, bpmsx] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPX, Q1, 'H', outputpath, beam_direction, accel, lhcphase)
        elif with_liny_for_zero_dppy:
            [phasey, Q2, MUY, bpmsy] = algorithms.helper.GetPhases(MADTwiss_ac, ListOfZeroDPPY, Q2, 'V', outputpath, beam_direction, accel, lhcphase)

#---- ac to free phase from eq and the model
    if with_ac_calc:
        if with_linx_for_zero_dppx:
            Q1f = Q1 - d1 #-- Free H-tune
            try:
                acphasex_ac2bpmac = algorithms.helper.GetACPhase_AC2BPMAC(MADTwissElem, Q1, Q1f, 'H', accel)
            except:
                acphasex_ac2bpmac = algorithms.helper.GetACPhase_AC2BPMAC(MADTwiss, Q1, Q1f, 'H', accel)
            [phasexf, muxf, bpmsxf] = algorithms.helper.GetFreePhase_Eq(MADTwiss, ListOfZeroDPPX, Q1, Q1f, acphasex_ac2bpmac, 'H', beam_direction, lhcphase)
            [phasexf2, muxf2, bpmsxf2] = algorithms.helper.getfreephase(phasex, Q1, Q1f, bpmsx, MADTwiss_ac, MADTwiss, "H")
        if with_liny_for_zero_dppy:
            Q2f = Q2 - d2 #-- Free V-tune
            try:
                acphasey_ac2bpmac = algorithms.helper.GetACPhase_AC2BPMAC(MADTwissElem, Q2, Q2f, 'V', accel)
            except:
                acphasey_ac2bpmac = algorithms.helper.GetACPhase_AC2BPMAC(MADTwiss, Q2, Q2f, 'V', accel)
            [phaseyf, muyf, bpmsyf] = algorithms.helper.GetFreePhase_Eq(MADTwiss, ListOfZeroDPPY, Q2, Q2f, acphasey_ac2bpmac, 'V', beam_direction, lhcphase)
            [phaseyf2, muyf2, bpmsyf2] = algorithms.helper.getfreephase(phasey, Q2, Q2f, bpmsy, MADTwiss_ac, MADTwiss, "V")

#---- H plane result
    if with_linx_for_zero_dppx:
        phasexlist = []
        phasex['DPP'] = 0.0
        phasexlist.append(phasex)
        tfs_file = files_dict['getphasex.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEX", "STDPHX", "PHXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpmsx)):
            bn1 = string.upper(bpmsx[i][1])
            bns1 = bpmsx[i][0]
            phmdl = phasex[bn1][4]
            if i == len(bpmsx) - 1:
                bn2 = string.upper(bpmsx[0][1])
                bns2 = bpmsx[0][0]
            else:
                bn2 = string.upper(bpmsx[i + 1][1])
                bns2 = bpmsx[i + 1][0]
            list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPX), phasex[bn1][0], phasex[bn1][1], phmdl, MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free phase
        if with_ac_calc:
            #-- from eq
            try:
                tfs_file = files_dict['getphasex_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEX", "STDPHX", "PHXMDL", "MUXMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(len(bpmsxf)):
                    bn1 = string.upper(bpmsxf[i][1])
                    bns1 = bpmsxf[i][0]
                    phmdlf = phasexf[bn1][4]
                    if i == len(bpmsxf) - 1:
                        bn2 = string.upper(bpmsxf[0][1])
                        bns2 = bpmsxf[0][0]
                    else:
                        bn2 = string.upper(bpmsxf[i + 1][1])
                        bns2 = bpmsxf[i + 1][0]
                    list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPX), phasexf[bn1][0], phasexf[bn1][1], phmdlf, MADTwiss.MUX[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            
            except:
                pass
            #-- from the model
            tfs_file = files_dict['getphasex_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEX", "STDPHX", "PHXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsxf2)):
                bn1 = string.upper(bpmsxf2[i][1])
                bns1 = bpmsxf2[i][0]
                phmdlf2 = phasexf2[bn1][2]
                bn2 = phasexf2[bn1][3]
                bns2 = phasexf2[bn1][4]
                list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPX), phasexf2[bn1][0], phasexf2[bn1][1], phmdlf2, MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
#---- V plane result
    if with_liny_for_zero_dppy:
        phaseylist = []
        phasey['DPP'] = 0.0
        phaseylist.append(phasey)
        tfs_file = files_dict['getphasey.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEY", "STDPHY", "PHYMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpmsy)):
            bn1 = string.upper(bpmsy[i][1])
            bns1 = bpmsy[i][0]
            phmdl = phasey[bn1][4]
            # TODO easier with modulo(vimaier)
            if i == len(bpmsy) - 1:
                bn2 = string.upper(bpmsy[0][1])
                bns2 = bpmsy[0][0]
            else:
                bn2 = string.upper(bpmsy[i + 1][1])
                bns2 = bpmsy[i + 1][0]
            list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPY), phasey[bn1][0], phasey[bn1][1], phmdl, MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free phase
        if with_ac_calc:
            #-- from eq
            try:
                tfs_file = files_dict['getphasey_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEY", "STDPHY", "PHYMDL", "MUYMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(len(bpmsyf)):
                    bn1 = string.upper(bpmsyf[i][1])
                    bns1 = bpmsyf[i][0]
                    phmdlf = phaseyf[bn1][4]
                    if i == len(bpmsyf) - 1:
                        bn2 = string.upper(bpmsyf[0][1])
                        bns2 = bpmsyf[0][0]
                    else:
                        bn2 = string.upper(bpmsyf[i + 1][1])
                        bns2 = bpmsyf[i + 1][0]
                    list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPY), phaseyf[bn1][0], phaseyf[bn1][1], phmdlf, MADTwiss.MUY[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            
            except:
                pass
            #-- from the model
            tfs_file = files_dict['getphasey_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEY", "STDPHY", "PHYMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsyf2)):
                bn1 = string.upper(bpmsyf2[i][1])
                bns1 = bpmsyf2[i][0]
                phmdlf2 = phaseyf2[bn1][2]
                bn2 = phaseyf2[bn1][3]
                bns2 = phaseyf2[bn1][4]
                list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPY), phaseyf2[bn1][0], phaseyf2[bn1][1], phmdlf2, MADTwiss.MUY[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
    return Q1, MUX, Q2, MUY, Q1f, acphasex_ac2bpmac, muxf, Q2f, muyf, muxf2, muyf2, acphasey_ac2bpmac, phasex, phasexf, phasey, phaseyf, bpmsx, bpmsy, phasexf2, phaseyf2, phasexlist, phaseylist
# END calculate_phase ------------------------------------------------------------------------------


def calculate_total_phase(accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, MUX, MUY, Q1f, Q2f, muxf, muyf, acphasex_ac2bpmac, muxf2, muyf2, acphasey_ac2bpmac):
    print 'Calculating total phase' 
#---- H plane result
    if with_linx_for_zero_dppx:
        [phasexT, bpmsxT] = algorithms.helper.GetPhasesTotal(MADTwiss_ac, ListOfZeroDPPX, Q1, 'H', beam_direction, accel, lhcphase)
        tfs_file = files_dict['getphasetotx.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEX", "STDPHX", "PHXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpmsxT)):
            bn1 = string.upper(bpmsxT[i][1])
            bns1 = bpmsxT[i][0]
            phmdl = phasexT[bn1][2]
            bn2 = string.upper(bpmsxT[0][1])
            bns2 = bpmsxT[0][0]
            list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPX), phasexT[bn1][0], phasexT[bn1][1], phmdl, MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free total phase
        if with_ac_calc:
            #-- from eq
            try:
                [phasexTf, bpmsxTf] = algorithms.helper.GetFreePhaseTotal_Eq(MADTwiss, ListOfZeroDPPX, Q1, Q1f, acphasex_ac2bpmac, 'H', beam_direction, lhcphase)
                tfs_file = files_dict['getphasetotx_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEX", "STDPHX", "PHXMDL", "MUXMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsxTf)):
                    bn1 = string.upper(bpmsxTf[i][1])
                    bns1 = bpmsxTf[i][0]
                    phmdlf = phasexTf[bn1][2]
                    bn2 = string.upper(bpmsxTf[0][1])
                    bns2 = bpmsxTf[0][0]
                    list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPX), phasexTf[bn1][0], phasexTf[bn1][1], phmdlf, MADTwiss.MUX[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except:
                
                print #TODO: debug remove
#                 pass
            #-- from the model
            [phasexTf2, bpmsxTf2] = algorithms.helper.getfreephaseTotal(phasexT, bpmsxT, "H", MADTwiss, MADTwiss_ac)
            tfs_file = files_dict['getphasetotx_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEX", "STDPHX", "PHXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsxTf2)):
                bn1 = string.upper(bpmsxTf2[i][1])
                bns1 = bpmsxTf2[i][0]
                phmdlf2 = phasexTf2[bn1][2]
                bn2 = string.upper(bpmsxTf2[0][1])
                bns2 = bpmsxTf2[0][0]
                list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPX), phasexTf2[bn1][0], phasexTf2[bn1][1], phmdlf2, MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
#---- V plane result
    if with_liny_for_zero_dppy:
        [phaseyT, bpmsyT] = algorithms.helper.GetPhasesTotal(MADTwiss_ac, ListOfZeroDPPY, Q2, 'V', beam_direction, accel, lhcphase)
        tfs_file = files_dict['getphasetoty.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEY", "STDPHY", "PHYMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpmsyT)):
            bn1 = string.upper(bpmsyT[i][1])
            bns1 = bpmsyT[i][0]
            phmdl = phaseyT[bn1][2]
            bn2 = string.upper(bpmsyT[0][1])
            bns2 = bpmsyT[0][0]
            list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPY), phaseyT[bn1][0], phaseyT[bn1][1], phmdl, MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free total phase
        if with_ac_calc:
            #-- from eq
            try:
                [phaseyTf, bpmsyTf] = algorithms.helper.GetFreePhaseTotal_Eq(MADTwiss, ListOfZeroDPPY, Q2, Q2f, acphasey_ac2bpmac, 'V', beam_direction, lhcphase)
                tfs_file = files_dict['getphasetoty_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEY", "STDPHY", "PHYMDL", "MUYMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsyTf)):
                    bn1 = string.upper(bpmsyTf[i][1])
                    bns1 = bpmsyTf[i][0]
                    phmdlf = phaseyTf[bn1][2]
                    bn2 = string.upper(bpmsyTf[0][1])
                    bns2 = bpmsyTf[0][0]
                    list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPY), phaseyTf[bn1][0], phaseyTf[bn1][1], phmdlf, MADTwiss.MUY[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except:
                pass
            #-- from the model
            [phaseyTf2, bpmsyTf2] = algorithms.helper.getfreephaseTotal(phaseyT, bpmsyT, "V", MADTwiss, MADTwiss_ac)
            tfs_file = files_dict['getphasetoty_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASEY", "STDPHY", "PHYMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsyTf2)):
                bn1 = string.upper(bpmsyTf2[i][1])
                bns1 = bpmsyTf2[i][0]
                phmdlf2 = phaseyTf2[bn1][2]
                bn2 = string.upper(bpmsyTf2[0][1])
                bns2 = bpmsyTf2[0][0]
                list_row_entries = ['"' + bn1 + '"', '"' + bn2 + '"', bns1, bns2, len(ListOfZeroDPPY), phaseyTf2[bn1][0], phaseyTf2[bn1][1], phmdlf2, MADTwiss.MUY[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
# END calculate_total_phase ------------------------------------------------------------------------

def calculate_beta(MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, Q1f, Q2f, phasex, phasexf, phasey, phaseyf):
    #TODO: initialize variables otherwise return stmt of this function would raie an exception(vimaier)
    betaxf = None
    betayf = None
    
    print 'Calculating beta'
#---- H plane
    if with_linx_for_zero_dppx:
        [betax, rmsbbx, alfax, bpms] = algorithms.helper.BetaFromPhase(MADTwiss_ac, ListOfZeroDPPX, phasex, 'H')
        betax['DPP'] = 0
        tfs_file = files_dict['getbetax.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbx))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), betax[bn1][0], betax[bn1][1], betax[bn1][2], alfax[bn1][0], alfax[bn1][1], alfax[bn1][2], MADTwiss_ac.BETX[MADTwiss_ac.indx[bn1]], MADTwiss_ac.ALFX[MADTwiss_ac.indx[bn1]], MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free beta
        if with_ac_calc:
            #-- from eq
            try:
                [betaxf, rmsbbxf, alfaxf, bpmsf] = algorithms.helper.BetaFromPhase(MADTwiss, ListOfZeroDPPX, phasexf, 'H')
                tfs_file = files_dict['getbetax_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = string.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), betaxf[bn1][0], betaxf[bn1][1], betaxf[bn1][2], alfaxf[bn1][0], alfaxf[bn1][1], alfaxf[bn1][2], MADTwiss.BETX[MADTwiss.indx[bn1]], MADTwiss.ALFX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except:
                pass
            #-- from the model
            [betaxf2, rmsbbxf2, alfaxf2, bpmsf2] = algorithms.helper.getFreeBeta(MADTwiss_ac, MADTwiss, betax, rmsbbx, alfax, bpms, 'H')
            tfs_file = files_dict['getbetax_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf2))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = string.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), betaxf2[bn1][0], betaxf2[bn1][1], betaxf2[bn1][2], alfaxf2[bn1][0], alfaxf2[bn1][1], alfaxf2[bn1][2], MADTwiss.BETX[MADTwiss.indx[bn1]], MADTwiss.ALFX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
#---- V plane
    if with_liny_for_zero_dppy:
        [betay, rmsbby, alfay, bpms] = algorithms.helper.BetaFromPhase(MADTwiss_ac, ListOfZeroDPPY, phasey, 'V')
        betay['DPP'] = 0
        tfs_file = files_dict['getbetay.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbby))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "ERRBETY", "STDBETY", "ALFY", "ERRALFY", "STDALFY", "BETYMDL", "ALFYMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), betay[bn1][0], betay[bn1][1], betay[bn1][2], alfay[bn1][0], alfay[bn1][1], alfay[bn1][2], MADTwiss_ac.BETY[MADTwiss_ac.indx[bn1]], MADTwiss_ac.ALFY[MADTwiss_ac.indx[bn1]], MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free beta
        if with_ac_calc:
            #-- from eq
            try:
                [betayf, rmsbbyf, alfayf, bpmsf] = algorithms.helper.BetaFromPhase(MADTwiss, ListOfZeroDPPY, phaseyf, 'V')
                tfs_file = files_dict['getbetay_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "ERRBETY", "STDBETY", "ALFY", "ERRALFY", "STDALFY", "BETYMDL", "ALFYMDL", "MUYMDL"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = string.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), betayf[bn1][0], betayf[bn1][1], betayf[bn1][2], alfayf[bn1][0], alfayf[bn1][1], alfayf[bn1][2], MADTwiss.BETY[MADTwiss.indx[bn1]], MADTwiss.ALFY[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except:
                pass
            #-- from the model
            [betayf2, rmsbbyf2, alfayf2, bpmsf2] = algorithms.helper.getFreeBeta(MADTwiss_ac, MADTwiss, betay, rmsbby, alfay, bpms, 'V')
            tfs_file = files_dict['getbetay_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf2))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "ERRBETY", "STDBETY", "ALFY", "ERRALFY", "STDALFY", "BETYMDL", "ALFYMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = string.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), betayf2[bn1][0], betayf2[bn1][1], betayf2[bn1][2], alfayf2[bn1][0], alfayf2[bn1][1], alfayf2[bn1][2], MADTwiss.BETY[MADTwiss.indx[bn1]], MADTwiss.ALFY[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
    return bpms, betax, betaxf, betay, betayf 
# END calculate_beta -------------------------------------------------------------------------------


def calculate_beta_from_amplitude(accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, Q1f, Q2f, acphasex_ac2bpmac, acphasey_ac2bpmac, bpms, betax, betaxf, betay, betayf):
    #TODO: initialize variables otherwise return stmt of this function would raise an exception(vimaier)
    betaxf_ratio = None
    betayf_ratio = None
    
    print 'Calculating beta from amplitude'
    betaxalist = []
    betayalist = []
    betaxPhaseCopy = betax #-- For Andy's BetaFromAmp re-scaling
    betayPhaseCopy = betay #-- For Andy's BetaFromAmp re-scaling
    betaxPhaseCopyf = betaxf #-- For Andy's BetaFromAmp re-scaling
    betayPhaseCopyf = betayf #-- For Andy's BetaFromAmp re-scaling
    #---- H plane
    if with_linx_for_zero_dppx:
        [betax, rmsbbx, bpms, invJx] = algorithms.helper.BetaFromAmplitude(MADTwiss_ac, ListOfZeroDPPX, 'H')
        betax['DPP'] = 0
        beta2_save = betax
        betaxalist.append(betax) #-- Rescaling
        betax_ratio = 0
        skipped_bpmx = []
        arcbpms = algorithms.helper.filterbpm(bpms)
        for bpm in arcbpms:
            name = string.upper(bpm[1]) # second entry is the name
        #Skip BPM with strange data
            if abs(betaxPhaseCopy[name][0] / betax[name][0]) > 100:
                skipped_bpmx.append(name)
            elif (betax[name][0] < 0 or betaxPhaseCopy[name][0] < 0):
                skipped_bpmx.append(name)
            else:
                betax_ratio = betax_ratio + (betaxPhaseCopy[name][0] / betax[name][0])
        
        try:
            betax_ratio = betax_ratio / (len(arcbpms) - len(skipped_bpmx))
        except:
            betax_ratio = 1
        betax_rescale = {}
        for bpm in map(string.upper, zip(*bpms)[1]):
            betax_rescale[bpm] = [betax_ratio * betax[bpm][0], betax_ratio * betax[bpm][1], betax[bpm][2]]
        
        tfs_file = files_dict['getampbetax.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbx))
        tfs_file.add_descriptor("RescalingFactor", "%le", str(betax_ratio))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), betax[bn1][0], betax[bn1][1], MADTwiss_ac.BETX[MADTwiss_ac.indx[bn1]], MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]], betax_rescale[bn1][0], betax_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries) #-- ac to free amp beta
        
        if with_ac_calc: #-- from eq
            try:
                [betaxf, rmsbbxf, bpmsf] = algorithms.helper.GetFreeBetaFromAmp_Eq(MADTwiss_ac, ListOfZeroDPPX, Q1, Q1f, acphasex_ac2bpmac, 'H', beam_direction, lhcphase)[:3]
                #-- Rescaling
                betaxf_ratio = 0
                skipped_bpmxf = []
                arcbpms = algorithms.helper.filterbpm(bpmsf)
                for bpm in arcbpms:
                    name = string.upper(bpm[1]) # second entry is the name
                #Skip BPM with strange data
                    if abs(betaxPhaseCopyf[name][0] / betaxf[name][0]) > 10:
                        skipped_bpmxf.append(name)
                    elif abs(betaxPhaseCopyf[name][0] / betaxf[name][0]) < 0.1:
                        skipped_bpmxf.append(name)
                    elif (betaxf[name][0] < 0 or betaxPhaseCopyf[name][0] < 0):
                        skipped_bpmxf.append(name)
                    else:
                        betaxf_ratio = betaxf_ratio + (betaxPhaseCopyf[name][0] / betaxf[name][0])
                
                try:
                    betaxf_ratio = betaxf_ratio / (len(arcbpms) - len(skipped_bpmxf))
                except:
                    betaxf_ratio = 1
                tfs_file = files_dict['getampbetax_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf))
                tfs_file.add_descriptor("RescalingFactor", "%le", str(betaxf_ratio))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = string.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), betaxf[bn1][0], betaxf[bn1][1], MADTwiss.BETX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]], betaxf_ratio * betaxf[bn1][0], betaxf_ratio * betaxf[bn1][1]]
                    tfs_file.add_table_row(list_row_entries) # Since invJxf(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            
            except:
                #-- from the model
                traceback.print_exc()
            # Since invJxf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betaxf2, rmsbbxf2, bpmsf2] = algorithms.helper.getFreeAmpBeta(betax, rmsbbx, bpms, invJx, MADTwiss_ac, MADTwiss, 'H')[:3]
            betaxf2_rescale = algorithms.helper.getFreeAmpBeta(betax_rescale, rmsbbx, bpms, invJx, MADTwiss_ac, MADTwiss, 'H')[0]
            tfs_file = files_dict['getampbetax_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf2))
            tfs_file.add_descriptor("RescalingFactor", "%le", str(betax_ratio))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = string.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), betaxf2[bn1][0], betaxf2[bn1][1], MADTwiss.BETX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]], betaxf2_rescale[bn1][0], betaxf2_rescale[bn1][1]]
                tfs_file.add_table_row(list_row_entries) #---- V plane
    
    if with_liny_for_zero_dppy:
        [betay, rmsbby, bpms, invJy] = algorithms.helper.BetaFromAmplitude(MADTwiss_ac, ListOfZeroDPPY, 'V')
        betay['DPP'] = 0
        betayalist.append(betay)
        #-- Rescaling
        betay_ratio = 0
        skipped_bpmy = []
        arcbpms = algorithms.helper.filterbpm(bpms)
        for bpm in arcbpms:
            name = string.upper(bpm[1]) # second entry is the name
        #Skip BPM with strange data
            if abs(betayPhaseCopy[name][0] / betay[name][0]) > 100:
                skipped_bpmy.append(name)
            elif (betay[name][0] < 0 or betayPhaseCopy[name][0] < 0):
                skipped_bpmy.append(name)
            else:
                betay_ratio = betay_ratio + (betayPhaseCopy[name][0] / betay[name][0])
        
        try:
            betay_ratio = betay_ratio / (len(arcbpms) - len(skipped_bpmy))
        except:
            betay_ratio = 1
        betay_rescale = {}
        for bpm in map(string.upper, zip(*bpms)[1]):
            betay_rescale[bpm] = [betay_ratio * betay[bpm][0], betay_ratio * betay[bpm][1], betay[bpm][2]]
        
        tfs_file = files_dict['getampbetay.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbby))
        tfs_file.add_descriptor("RescalingFactor", "%le", str(betay_ratio))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), betay[bn1][0], betay[bn1][1], MADTwiss_ac.BETY[MADTwiss_ac.indx[bn1]], MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]], betay_rescale[bn1][0], betay_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries) #-- ac to free amp beta
        
        if with_ac_calc: #-- from eq
            try:
                [betayf, rmsbbyf, bpmsf] = algorithms.helper.GetFreeBetaFromAmp_Eq(MADTwiss_ac, ListOfZeroDPPY, Q2, Q2f, acphasey_ac2bpmac, 'V', beam_direction, accel)[:3] #-- Rescaling
                betayf_ratio = 0
                skipped_bpmyf = []
                arcbpms = algorithms.helper.filterbpm(bpmsf)
                for bpm in arcbpms:
                    name = string.upper(bpm[1]) # second entry is the name
                #Skip BPM with strange data
                    if abs(betayPhaseCopyf[name][0] / betayf[name][0]) > 10:
                        skipped_bpmyf.append(name)
                    elif (betayf[name][0] < 0 or betayPhaseCopyf[name][0] < 0):
                        skipped_bpmyf.append(name)
                    elif abs(betayPhaseCopyf[name][0] / betayf[name][0]) < 0.1:
                        skipped_bpmyf.append(name)
                    else:
                        betayf_ratio = betayf_ratio + (betayPhaseCopyf[name][0] / betayf[name][0])
                
                try:
                    betayf_ratio = betayf_ratio / (len(arcbpms) - len(skipped_bpmyf))
                except:
                    betayf_ratio = 1
                tfs_file = files_dict['getampbetay_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf))
                tfs_file.add_descriptor("RescalingFactor", "%le", str(betayf_ratio))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = string.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), betayf[bn1][0], betayf[bn1][1], MADTwiss.BETY[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]], (betayf_ratio * betayf[bn1][0]), (betayf_ratio * betayf[bn1][1])]
                    tfs_file.add_table_row(list_row_entries) # 'except ALL' catched a SystemExit from filterbpm().(vimaier)
            
                # Since invJyf(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            except SystemExit:
                sys.exit(1)
            except:
                #-- from the model
                traceback.print_exc()
            # Since invJyf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betayf2, rmsbbyf2, bpmsf2] = algorithms.helper.getFreeAmpBeta(betay, rmsbby, bpms, invJy, MADTwiss_ac, MADTwiss, 'V')[:3]
            betayf2_rescale = algorithms.helper.getFreeAmpBeta(betay_rescale, rmsbby, bpms, invJy, MADTwiss_ac, MADTwiss, 'V')[0]
            tfs_file = files_dict['getampbetay_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf2))
            tfs_file.add_descriptor("RescalingFactor", "%le", str(betay_ratio))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = string.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), betayf2[bn1][0], betayf2[bn1][1], MADTwiss.BETY[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]], betayf2_rescale[bn1][0], betayf2_rescale[bn1][1]]
                tfs_file.add_table_row(list_row_entries)
    
    return betax, betay, bpms, beta2_save, betaxalist, betayalist, betax_ratio, betay_ratio, betaxf_ratio, betayf_ratio
# END calculate_beta_from_amplitude ----------------------------------------------------------------


def calculate_ip(accel, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, Q1f, Q2f, bpm_name, acphasex_ac2bpmac, acphasey_ac2bpmac, phasex, phasexf, phasey, phaseyf, bpmsx, bpmsy, phasexf2, phaseyf2, betax, betay):
    print 'Calculating IP'
    if "LHC" in accel:
        tfs_file = files_dict['getIP.out']
        tfs_file.add_column_names(["NAME", "BETASTARH", "BETASTARHMDL", "H", "PHIH", "PHIXH", "PHIHMDL", "BETASTARV", "BETASTARVMDL", "V", "PHIV", "PHIYV", "PHIVMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        ips = ["1", "2", "3", "4", "5", "6", "7", "8"]
        try:
            measured = [betax, betay]
            phases = [phasex, phasey]
            bpmss = [bpmsx, bpmsy]
        except:
            pass
        for ip in ips:
            try:
                betahor, betaver = algorithms.helper.getIP(ip, measured, MADTwiss, phases, bpmss) #print "entering"
            except:
                betahor = [0, 0, 0, 0, 0, 0, 0]
                betaver = [0, 0, 0, 0, 0, 0, 0]
            #print str(betahor[6])
            list_row_entries = ['"IP' + ip + '"', betahor[1], betahor[4], betahor[2], betahor[3], betahor[6], betahor[5], betaver[1], betaver[4], betaver[2], betaver[3], betaver[6], betaver[5]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- Parameters at IP1, IP2, IP5, and IP8
        IPx = algorithms.helper.GetIP2(MADTwiss_ac, ListOfZeroDPPX, Q1, 'H', beam_direction, accel, lhcphase)
        IPy = algorithms.helper.GetIP2(MADTwiss_ac, ListOfZeroDPPY, Q2, 'V', beam_direction, accel, lhcphase)
        tfs_file_x = files_dict['getIPx.out']
        tfs_file_x.add_column_names(["NAME", "BETX", "BETXSTD", "BETXMDL", "ALFX", "ALFXSTD", "ALFXMDL", "BETX*", "BETX*STD", "BETX*MDL", "SX*", "SX*STD", "SX*MDL", "rt(2JX)", "rt(2JX)STD"])
        tfs_file_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        tfs_file_y = files_dict['getIPy.out']
        tfs_file_y.add_column_names(["NAME", "BETY", "BETYSTD", "BETYMDL", "ALFY", "ALFYSTD", "ALFYMDL", "BETY*", "BETY*STD", "BETY*MDL", "SY*", "SY*STD", "SY*MDL", "rt(2JY)", "rt(2JY)STD"])
        tfs_file_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for bpm_name in 'IP1', 'IP5', 'IP8', 'IP2':
            try:
                list_row_entries = ['"' + bpm_name + '"']
                for k in IPx[bpm_name]:
                    list_row_entries.append(k)
                
                tfs_file_x.add_table_row(list_row_entries)
            except:
                pass
            try:
                list_row_entries = ['"' + bpm_name + '"']
                for k in IPy[bpm_name]:
                    list_row_entries.append(k)
                
                tfs_file_y.add_table_row(list_row_entries)
            except:
                pass
        #-- ac to free parameters at IP1, IP2, IP5, and IP8
        if with_ac_calc:
            #-- From Eq
            IPxf = algorithms.helper.GetFreeIP2_Eq(MADTwiss, ListOfZeroDPPX, Q1, Q1f, acphasex_ac2bpmac, 'H', beam_direction, accel, lhcphase)
            IPyf = algorithms.helper.GetFreeIP2_Eq(MADTwiss, ListOfZeroDPPY, Q2, Q2f, acphasey_ac2bpmac, 'V', beam_direction, accel, lhcphase)
            tfs_file_x = files_dict['getIPx_free.out']
            tfs_file_x.add_column_names(["NAME", "BETX", "BETXSTD", "BETXMDL", "ALFX", "ALFXSTD", "ALFXMDL", "BETX*", "BETX*STD", "BETX*MDL", "SX*", "SX*STD", "SX*MDL", "rt(2JXD)", "rt(2JXD)STD"])
            tfs_file_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            tfs_file_y = files_dict['getIPy_free.out']
            tfs_file_y.add_column_names(["NAME", "BETY", "BETYSTD", "BETYMDL", "ALFY", "ALFYSTD", "ALFYMDL", "BETY*", "BETY*STD", "BETY*MDL", "SY*", "SY*STD", "SY*MDL", "rt(2JYD)", "rt(2JYD)STD"])
            tfs_file_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for bpm_name in 'IP1', 'IP5', 'IP8', 'IP2':
                try:
                    list_row_entries = ['"' + bpm_name + '"']
                    for k in IPxf[bpm_name]:
                        list_row_entries.append(k)
                    
                    tfs_file_x.add_table_row(list_row_entries)
                except:
                    pass
                try:
                    list_row_entries = ['"' + bpm_name + '"']
                    for k in IPyf[bpm_name]:
                        list_row_entries.append(k)
                    
                    tfs_file_y.add_table_row(list_row_entries)
                except:
                    pass
            #-- From model
            IPxf2 = algorithms.helper.GetFreeIP2(MADTwiss, MADTwiss_ac, IPx, 'H', accel)
            IPyf2 = algorithms.helper.GetFreeIP2(MADTwiss, MADTwiss_ac, IPy, 'V', accel)
            tfs_file_x = files_dict['getIPx_free2.out']
            tfs_file_x.add_column_names(["NAME", "BETX", "BETXSTD", "BETXMDL", "ALFX", "ALFXSTD", "ALFXMDL", "BETX*", "BETX*STD", "BETX*MDL", "SX*", "SX*STD", "SX*MDL", "rt(2JXD)", "rt(2JXD)STD"])
            tfs_file_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            tfs_file_y = files_dict['getIPy_free2.out']
            tfs_file_y.add_column_names(["NAME", "BETY", "BETYSTD", "BETYMDL", "ALFY", "ALFYSTD", "ALFYMDL", "BETY*", "BETY*STD", "BETY*MDL", "SY*", "SY*STD", "SY*MDL", "rt(2JYD)", "rt(2JYD)STD"])
            tfs_file_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for bpm_name in 'IP1', 'IP5', 'IP8', 'IP2':
                try:
                    list_row_entries = ['"' + bpm_name + '"']
                    for k in IPxf2[bpm_name]:
                        list_row_entries.append(k)
                    
                    tfs_file_x.add_table_row(list_row_entries)
                except:
                    pass
                try:
                    list_row_entries = ['"' + bpm_name + '"']
                    for k in IPyf2[bpm_name]:
                        list_row_entries.append(k)
                    
                    tfs_file_y.add_table_row(list_row_entries)
                except:
                    pass

        #-- IP beta* and phase from phase only
        try:
            IPfromphase = algorithms.helper.GetIPFromPhase(MADTwiss_ac, phasex, phasey, accel)
        except:
            print 'No output from IP from phase. H or V file missing?'
        tfs_file = files_dict['getIPfromphase.out']
        tfs_file.add_column_names(["NAME", "2L", "BETX*", "BETX*STD", "BETX*MDL", "BETY*", "BETY*STD", "BETY*MDL", "PHX", "PHXSTD", "PHXMDL", "PHY", "PHYSTD", "PHYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for bpm_name in 'IP1', 'IP5', 'IP8', 'IP2':
            list_row_entries = ['"' + bpm_name + '"']
            try:
                for k in IPfromphase[bpm_name]:
                    list_row_entries.append(k)
                
                tfs_file.add_table_row(list_row_entries)
            except:
                pass
        #-- ac to free beta*
        if with_ac_calc:
            #-- from eqs
            try:
                IPfromphasef = algorithms.helper.GetIPFromPhase(MADTwiss, phasexf, phaseyf, accel)
            except:
                pass
            tfs_file = files_dict['getIPfromphase_free.out']
            tfs_file.add_column_names(["NAME", "2L", "BETX*", "BETX*STD", "BETX*MDL", "BETY*", "BETY*STD", "BETY*MDL", "PHX", "PHXSTD", "PHXMDL", "PHY", "PHYSTD", "PHYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for bpm_name in 'IP1', 'IP5', 'IP8', 'IP2':
                list_row_entries = ['"' + bpm_name + '"']
                try:
                    for k in IPfromphasef[bpm_name]:
                        list_row_entries.append(k)
                    
                    tfs_file.add_table_row(list_row_entries)
                except:
                    pass
            #-- from the model
            try:
                IPfromphasef2 = algorithms.helper.GetIPFromPhase(MADTwiss, phasexf2, phaseyf2, accel)
            except:
                pass
            tfs_file = files_dict['getIPfromphase_free2.out']
            tfs_file.add_column_names(["NAME", "2L", "BETX*", "BETX*STD", "BETX*MDL", "BETY*", "BETY*STD", "BETY*MDL", "PHX", "PHXSTD", "PHXMDL", "PHY", "PHYSTD", "PHYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for bpm_name in 'IP1', 'IP5', 'IP8', 'IP2':
                list_row_entries = ['"' + bpm_name + '"']
                try:
                    for k in IPfromphasef2[bpm_name]:
                        list_row_entries.append(k)
                    
                    tfs_file.add_table_row(list_row_entries)
                except:
                    pass
# END calculate_ip ---------------------------------------------------------------------------------


def calculate_orbit(twiss_model_file, accel, MADTwiss, files_dict, FileOfNonZeroDPPX, FileOfNonZeroDPPY, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, with_liny_for_zero_dppy, with_liny_for_nonzero_dppy, with_linx_for_zero_dppx, with_linx_for_nonzero_dppx, Q1, Q2, bpms):
    ListOfCOX = []
    if with_linx_for_zero_dppx:
        [cox, bpms] = algorithms.helper.GetCO(MADTwiss, ListOfZeroDPPX)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOx.out']
        tfs_file.add_descriptor("TABLE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("TYPE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("SEQUENCE", "%05s", '"' + accel + '"')
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "X", "STDX", "XMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), cox[bn1][0], cox[bn1][1], MADTwiss.X[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        ListOfCOX.append(cox)
    ListOfCOY = []
    if with_liny_for_zero_dppy:
        [coy, bpms] = algorithms.helper.GetCO(MADTwiss, ListOfZeroDPPY)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOy.out']
        tfs_file.add_descriptor("TABLE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("TYPE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("SEQUENCE", "%05s", '"' + accel + '"')
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "Y", "STDY", "YMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), coy[bn1][0], coy[bn1][1], MADTwiss.Y[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        ListOfCOY.append(coy)
#-------- Orbit for non-zero DPP
    if with_linx_for_nonzero_dppx:
        k = 0
        for j in ListOfNonZeroDPPX:
            list_with_single_twiss = []
            list_with_single_twiss.append(j)
            filename = 'getCOx_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPX[k])
            tfs_file.add_descriptor("DPP", "%le", str(float(j.DPP)))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            [codpp, bpms] = algorithms.helper.GetCO(MADTwiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "X", "STDX", "XMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), codpp[bn1][0], codpp[bn1][1], MADTwiss.X[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
            
            ListOfCOX.append(codpp)
            k += 1
    
    if with_liny_for_nonzero_dppy:
        k = 0
        for j in ListOfNonZeroDPPY:
            list_with_single_twiss = []
            list_with_single_twiss.append(j)
            filename = 'getCOy_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPY[k])
            tfs_file.add_descriptor("DPP", "%le", str(float(j.DPP)))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            [codpp, bpms] = algorithms.helper.GetCO(MADTwiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "Y", "STDY", "YMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                #TODO: why ListOfZeroDPPY.. above used ListOfNonZeroDPPY(vimaier)
                list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPY), codpp[bn1][0], codpp[bn1][1], MADTwiss.Y[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
            
            ListOfCOY.append(codpp)
            k += 1
    
    return ListOfCOX, ListOfCOY
# END calculate_orbit ------------------------------------------------------------------------------

def calculate_dispersion(BPMU, COcut, MADTwiss, files_dict, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, with_liny_for_zero_dppy, with_liny_for_nonzero_dppy, with_linx_for_zero_dppx, with_linx_for_nonzero_dppx, Q1, Q2, bpms, beta2_save, ListOfCOX, ListOfCOY):
    if with_linx_for_zero_dppx and with_linx_for_nonzero_dppx:
        [nda, Dx, DPX, bpms] = algorithms.helper.NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, beta2_save, COcut)
        tfs_file = files_dict['getNDx.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "NDX", "STDNDX", "DX", "DPX", "NDXMDL", "DXMDL", "DPXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            ndmdl = MADTwiss.DX[MADTwiss.indx[bn1]] / math.sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfNonZeroDPPX), nda[bn1][0], nda[bn1][1], Dx[bn1][0], DPX[bn1], ndmdl, MADTwiss.DX[MADTwiss.indx[bn1]], MADTwiss.DPX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        [dxo, bpms] = algorithms.helper.DispersionfromOrbit(ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, COcut, BPMU)
        DPX = algorithms.helper.GetDPX(MADTwiss, dxo, bpms)
        tfs_file = files_dict['getDx.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "DX", "STDDX", "DPX", "DXMDL", "DPXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfNonZeroDPPX), dxo[bn1][0], dxo[bn1][1], DPX[bn1], MADTwiss.DX[MADTwiss.indx[bn1]], MADTwiss.DPX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
    
    if with_liny_for_zero_dppy and with_liny_for_nonzero_dppy:
        [dyo, bpms] = algorithms.helper.DispersionfromOrbit(ListOfZeroDPPY, ListOfNonZeroDPPY, ListOfCOY, COcut, BPMU)
        DPY = algorithms.helper.GetDPY(MADTwiss, dyo, bpms)
        tfs_file = files_dict['getDy.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "DY", "STDDY", "DPY", "DYMDL", "DPYMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfNonZeroDPPY), dyo[bn1][0], dyo[bn1][1], DPY[bn1], MADTwiss.DY[MADTwiss.indx[bn1]], MADTwiss.DPY[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
# END calculate_dispersion -------------------------------------------------------------------------

def calculate_coupling(outputpath, accel, lhcphase, NBcpl, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfZeroDPPY, with_liny_for_zero_dppy, with_linx_for_zero_dppx, Q1, Q2, Q1f, Q2f, PseudoListX, PseudoListY, acphasex_ac2bpmac, acphasey_ac2bpmac, phasexlist, phaseylist, bpms):
    print "Calculating coupling"

    if with_linx_for_zero_dppx and with_liny_for_zero_dppy:
        #-- Coupling in the model
        try:
            MADTwiss.Cmatrix()
        except:
            pass
        #-- Main part
        if NBcpl == 1:
            # Avoids crashing the programm(vimaier)
            fwqwf = None
            fwqwf2 = None
            [fwqw, bpms] = algorithms.helper.GetCoupling1(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2)
            tfs_file = files_dict['getcouple.out']
            tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
            tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD", "Q1001W", "QWSTD", "MDLF1001R", "MDLF1001I"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                try:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), (math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2)), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, MADTwiss.f1001[MADTwiss.indx(bn1)].real, MADTwiss.f1001[MADTwiss.indx(bn1)].imag, MADTwiss_ac.f1010[MADTwiss_ac.indx(bn1)].real, MADTwiss_ac.f1010[MADTwiss_ac.indx(bn1)].imag]
                #-- Output zero if the model does not have couping parameters
                except:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), (math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2)), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, 0.0, 0.0]
                tfs_file.add_table_row(list_row_entries)
        
        elif NBcpl == 2:
            if accel == "SPS" or "RHIC" in accel:
                [phasexp, Q1, MUX, bpmsx] = algorithms.helper.GetPhases(MADTwiss, PseudoListX, 'H', outputpath, beam_direction, accel, lhcphase)
                [phaseyp, Q2, MUY, bpmsy] = algorithms.helper.GetPhases(MADTwiss, PseudoListY, 'V', outputpath, beam_direction, accel, lhcphase)
                [fwqw, bpms] = algorithms.helper.GetCoupling2(MADTwiss, PseudoListX, PseudoListY, Q1, Q2, phasexp, phaseyp, beam_direction, accel)
            else:
                [fwqw, bpms] = algorithms.helper.GetCoupling2(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, phasexlist[0], phaseylist[0], beam_direction, accel)
            tfs_file = files_dict['getcouple.out']
            tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
            tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD1", "F1001R", "F1001I", "F1010W", "FWSTD2", "F1010R", "F1010I", "MDLF1001R", "MDLF1001I", "MDLF1010R", "MDLF1010I"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                try:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), (math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2)), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, math.sqrt(fwqw[bn1][0][2].real ** 2 + fwqw[bn1][0][2].imag ** 2), fwqw[bn1][0][3], fwqw[bn1][0][2].real, fwqw[bn1][0][2].imag, MADTwiss_ac.f1001[MADTwiss_ac.indx[bn1]].real, MADTwiss_ac.f1001[MADTwiss_ac.indx[bn1]].imag, MADTwiss_ac.f1010[MADTwiss_ac.indx[bn1]].real, MADTwiss_ac.f1010[MADTwiss_ac.indx[bn1]].imag]
                #-- Output zero if the model does not have couping parameters
                except:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, math.sqrt(fwqw[bn1][0][2].real ** 2 + fwqw[bn1][0][2].imag ** 2), fwqw[bn1][0][3], fwqw[bn1][0][2].real, fwqw[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0]
                tfs_file.add_table_row(list_row_entries)
            
            #-- ac to free coupling
            if with_ac_calc:
                #-- analytic eqs
                try:
                    [fwqwf, bpmsf] = algorithms.helper.GetFreeCoupling_Eq(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, Q1f, Q2f, acphasex_ac2bpmac, acphasey_ac2bpmac, beam_direction)
                    tfs_file = files_dict['getcouple_free.out']
                    tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
                    tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
                    tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD1", "F1001R", "F1001I", "F1010W", "FWSTD2", "F1010R", "F1010I", "MDLF1001R", "MDLF1001I", "MDLF1010R", "MDLF1010I"])
                    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                    for i in range(len(bpmsf)):
                        bn1 = string.upper(bpmsf[i][1])
                        bns1 = bpmsf[i][0]
                        try:
                            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), math.sqrt(fwqwf[bn1][0][0].real ** 2 + fwqwf[bn1][0][0].imag ** 2), fwqwf[bn1][0][1], fwqwf[bn1][0][0].real, fwqwf[bn1][0][0].imag, math.sqrt(fwqwf[bn1][0][2].real ** 2 + fwqwf[bn1][0][2].imag ** 2), fwqwf[bn1][0][3], fwqwf[bn1][0][2].real, fwqwf[bn1][0][2].imag, MADTwiss.f1001[MADTwiss.indx[bn1]].real, MADTwiss.f1001[MADTwiss.indx[bn1]].imag, MADTwiss.f1010[MADTwiss.indx[bn1]].real, MADTwiss.f1010[MADTwiss.indx[bn1]].imag] #-- Output zero if the model does not have couping parameters
                        except:
                            list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), math.sqrt(fwqwf[bn1][0][0].real ** 2 + fwqwf[bn1][0][0].imag ** 2), fwqwf[bn1][0][1], fwqwf[bn1][0][0].real, fwqwf[bn1][0][0].imag, math.sqrt(fwqwf[bn1][0][2].real ** 2 + fwqwf[bn1][0][2].imag ** 2), fwqwf[bn1][0][3], fwqwf[bn1][0][2].real, fwqwf[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0]
                        tfs_file.add_table_row(list_row_entries)
                
                except:
                    traceback.print_exc()
                    pass
                #-- global factor
                [fwqwf2, bpmsf2] = algorithms.helper.getFreeCoupling(Q1f, Q2f, Q1, Q2, fwqw, MADTwiss, bpms)
                tfs_file = files_dict['getcouple_free2.out']
                tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
                tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD1", "F1001R", "F1001I", "F1010W", "FWSTD2", "F1010R", "F1010I", "MDLF1001R", "MDLF1001I", "MDLF1010R", "MDLF1010I"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(len(bpmsf2)):
                    bn1 = string.upper(bpmsf2[i][1])
                    bns1 = bpmsf2[i][0]
                    try:
                        list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), math.sqrt(fwqwf2[bn1][0][0].real ** 2 + fwqwf2[bn1][0][0].imag ** 2), fwqwf2[bn1][0][1], fwqwf2[bn1][0][0].real, fwqwf2[bn1][0][0].imag, math.sqrt(fwqwf2[bn1][0][2].real ** 2 + fwqwf2[bn1][0][2].imag ** 2), fwqwf2[bn1][0][3], fwqwf2[bn1][0][2].real, fwqwf2[bn1][0][2].imag, MADTwiss.f1001[MADTwiss.indx[bn1]].real, MADTwiss.f1001[MADTwiss.indx[bn1]].imag, MADTwiss.f1010[MADTwiss.indx[bn1]].real, MADTwiss.f1010[MADTwiss.indx[bn1]].imag] #-- Output zero if the model does not have couping parameters
                    except:
                        list_row_entries = ['"' + bn1 + '"', bns1, len(ListOfZeroDPPX), math.sqrt(fwqwf2[bn1][0][0].real ** 2 + fwqwf2[bn1][0][0].imag ** 2), fwqwf2[bn1][0][1], fwqwf2[bn1][0][0].real, fwqwf2[bn1][0][0].imag, math.sqrt(fwqwf2[bn1][0][2].real ** 2 + fwqwf2[bn1][0][2].imag ** 2), fwqwf2[bn1][0][3], fwqwf2[bn1][0][2].real, fwqwf2[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0]
                    tfs_file.add_table_row(list_row_entries)
        
        else:
            raise ValueError('Number of monitors for coupling analysis should be 1 or 2 (option -n)')
        #-- Convert to C-matrix:
        if with_ac_calc and (fwqwf is not None or fwqwf2 is not None):
            try:
                [coupleterms, Qminav, Qminerr, bpms] = algorithms.helper.getCandGammaQmin(fwqwf, bpmsf, Q1f, Q2f, MADTwiss)
            except:
                [coupleterms, Qminav, Qminerr, bpms] = algorithms.helper.getCandGammaQmin(fwqwf2, bpmsf2, Q1f, Q2f, MADTwiss)
        else:
            [coupleterms, Qminav, Qminerr, bpms] = algorithms.helper.getCandGammaQmin(fwqw, bpms, Q1f, Q2f, MADTwiss)
        tfs_file = files_dict['getcoupleterms.out']
        tfs_file.add_descriptor("DQMIN", "%le", Qminav)
        tfs_file.add_descriptor("DQMINE", "%le", Qminerr)
        tfs_file.add_column_names(["NAME", "S", "DETC", "DETCE", "GAMMA", "GAMMAE", "C11", "C12", "C21", "C22"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for bpm in bpms:
            bps = bpm[0]
            bpmm = bpm[1].upper()
            list_row_entries = [bpmm, bps, coupleterms[bpmm][0], coupleterms[bpmm][1], coupleterms[bpmm][2], coupleterms[bpmm][3], coupleterms[bpmm][4], coupleterms[bpmm][5], coupleterms[bpmm][6], coupleterms[bpmm][7]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- For chromatic coupling
        fwqw['DPP'] = 0
    return Q1, Q2
# END calculate_coupling ---------------------------------------------------------------------------

def phase_and_beta_for_non_zero_dpp(outputpath, twiss_model_file, accel, lhcphase, NBcpl, BPMdictionary, beam_direction, MADTwiss, with_ac_calc, files_dict, FileOfNonZeroDPPX, FileOfNonZeroDPPY, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, with_liny_for_nonzero_dppy, with_linx_for_nonzero_dppx, Q1, Q2, Q1f, Q2f, PseudoListX, PseudoListY, phasex, phasey, phasexlist, phaseylist, bpms, betax, betay, betaxalist, betayalist):
    print "Calculating phase and Beta for non-zero DPP" 
    #TODO: what is a thingie??(vimaier)
    print "lenght of zerothingie " + str(len(ListOfNonZeroDPPX))
    print "lenght of zerothingie " + str(len(ListOfNonZeroDPPY))
    if with_linx_for_nonzero_dppx:
        plane = 'H'
        k = 0
        for j in ListOfNonZeroDPPX:
            dpop = float(j.DPP)
            list_with_single_twiss = []
            list_with_single_twiss.append(j)
            DPPTwiss = algorithms.helper.ConstructOffMomentumModel(MADTwiss, dpop, BPMdictionary)
            [phasex, Q1DPP, MUX, bpms] = algorithms.helper.GetPhases(DPPTwiss, list_with_single_twiss, Q1, plane, outputpath, beam_direction, accel, lhcphase)
            phasex['DPP'] = dpop
            phasexlist.append(phasex)
            # 'getphasex_dpp_'+str(k+1)+'.out' was inteded for debugging (vimaier)
            if DEBUG:
                filename = 'getphasex_dpp_' + str(k + 1) + '.out'
                files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
                tfs_file = files_dict[filename]
                tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPX[k])
                tfs_file.add_descriptor("DPP", "%le", str(dpop))
                tfs_file.add_descriptor("Q1", "%le", str(Q1))
                tfs_file.add_descriptor("Q2", "%le", str(Q2))
                tfs_file.add_descriptor("Q1DPP", "%le", str(Q1DPP))
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASE", "STDPH", "PHXMDL", "MUXMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                length_bpms = len(bpms)
                for i in range(0, length_bpms):
                    bn1 = string.upper(bpms[i][1])
                    bns1 = bpms[i][0]
                    index = (i + 1) % length_bpms
                    bn2 = string.upper(bpms[index][1])
                    bns2 = bpms[index][0]
                    try:
                        phmdl = phasexlist[0][bn1][4]
                    except:
                        phmdl = 0.0
                    #phmdl=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
                    list_row_entries = ['"' + bn1 + '" ', '"' + bn2 + '"', bns1, bns2, 1, phasex[bn1][0], phasex[bn1][1], phmdl, MADTwiss.MUX[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            
            betax = {}
            alfax = {}
            [betax, rmsbbx, alfax, bpms] = algorithms.helper.BetaFromPhase(MADTwiss, list_with_single_twiss, phasex, plane)
            betax['DPP'] = dpop
            betaxa = {}
            [betaxa, rmsbbx, bpms, invJx] = algorithms.helper.BetaFromAmplitude(MADTwiss, list_with_single_twiss, plane)
            betaxa['DPP'] = dpop
            betaxalist.append(betaxa)
            filename = 'getbetax_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPX[k])
            tfs_file.add_descriptor("DPP", "%le", str(dpop))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '" ', bns1, len(ListOfZeroDPPX), betax[bn1][0], betax[bn1][1], betax[bn1][2], alfax[bn1][0], alfax[bn1][1], alfax[bn1][2], MADTwiss.BETX[MADTwiss.indx[bn1]], MADTwiss.ALFX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
            
            k += 1
    
    if with_liny_for_nonzero_dppy:
        plane = 'V'
        k = 0
        for j in ListOfNonZeroDPPY:
            dpop = float(j.DPP)
            list_with_single_twiss = []
            list_with_single_twiss.append(j)
            DPPTwiss = algorithms.helper.ConstructOffMomentumModel(MADTwiss, dpop, BPMdictionary)
            [phasey, Q2DPP, MUY, bpms] = algorithms.helper.GetPhases(DPPTwiss, list_with_single_twiss, Q2, plane, outputpath, beam_direction, accel, lhcphase)
            phasey['DPP'] = dpop
            phaseylist.append(phasey)
            # 'getphasex_dpp_'+str(k+1)+'.out' was inteded for debugging (vimaier)
            if DEBUG:
                filename = 'getphasey_dpp_' + str(k + 1) + '.out'
                files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
                tfs_file = files_dict[filename]
                tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPY[k])
                tfs_file.add_descriptor("Q1", "%le", str(Q1))
                tfs_file.add_descriptor("Q2", "%le", str(Q2))
                tfs_file.add_descriptor("Q2DPP", "%le", str(Q2DPP))
                tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "COUNT", "PHASE", "STDPH", "PHYMDL", "MUYMDL"])
                tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                length_bpms = len(bpms)
                for i in range(0, length_bpms):
                    bn1 = string.upper(bpms[i][1])
                    bns1 = bpms[i][0]
                    index = (i + 1) % length_bpms
                    bn2 = string.upper(bpms[index][1])
                    bns2 = bpms[index][0]
                    try:
                        phmdl = phaseylist[0][bn1][4]
                    except:
                        phmdl = 0.0
                    #phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
                    list_row_entries = ['"' + bn1 + '" ', '"' + bn2 + '" ', bns1, bns2, 1, phasey[bn1][0], phasey[bn1][1], phmdl, MADTwiss.MUY[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            
            betay = {}
            alfay = {}
            [betay, rmsbby, alfay, bpms] = algorithms.helper.BetaFromPhase(DPPTwiss, list_with_single_twiss, phasey, plane)
            betay['DPP'] = dpop
            betaya = {}
            [betaya, rmsbby, bpms, invJy] = algorithms.helper.BetaFromAmplitude(DPPTwiss, list_with_single_twiss, plane)
            betaya['DPP'] = dpop
            betayalist.append(betaya)
            filename = 'getbetay_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPY[k])
            tfs_file.add_descriptor("DPP", "%le", str(dpop))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            #TODO: check if it should be Y instead of X in the column names since it is getbetaY_dpp...out (vimaier)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '" ', bns1, len(ListOfZeroDPPY), betay[bn1][0], betay[bn1][1], betay[bn1][2], alfay[bn1][0], alfay[bn1][1], alfay[bn1][2], MADTwiss.BETX[MADTwiss.indx[bn1]], MADTwiss.ALFX[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
            
            k += 1
    
    if with_liny_for_nonzero_dppy and with_linx_for_nonzero_dppx:
        if len(ListOfNonZeroDPPX) != len(ListOfNonZeroDPPY):
            raise ValueError("list of dppx is not equal list of dppy")
        for j in range(len(ListOfNonZeroDPPX)):
            dpop = float(ListOfNonZeroDPPX[j].DPP)
            list_with_single_twiss_x = []
            list_with_single_twiss_y = []
            list_with_single_twiss_x.append(ListOfNonZeroDPPX[j])
            list_with_single_twiss_y.append(ListOfNonZeroDPPY[j])
            ### coupling
            try:
                MADTwiss.Cmatrix()
            except:
                0.0
            if accel == "SPS" or "RHIC" in accel:
                plane = 'H'
                [phasexp, Q1, MUX, bpmsx] = algorithms.helper.GetPhases(MADTwiss, PseudoListX, plane, outputpath, beam_direction, accel, lhcphase)
                plane = 'V'
                [phaseyp, Q2, MUY, bpmsy] = algorithms.helper.GetPhases(MADTwiss, PseudoListY, plane, outputpath, beam_direction, accel, lhcphase)
                [fwqw, bpms] = algorithms.helper.GetCoupling2(MADTwiss, PseudoListX, PseudoListY, Q1, Q2, phasexp, phaseyp, beam_direction, accel)
            elif NBcpl == 1:
                [fwqw, bpms] = algorithms.helper.GetCoupling1(MADTwiss, list_with_single_twiss_x, list_with_single_twiss_y, Q1, Q2)
            elif NBcpl == 2:
                print phasexlist[j + 1]['DPP'], dpop
                [fwqw, bpms] = algorithms.helper.GetCoupling2(MADTwiss, list_with_single_twiss_x, list_with_single_twiss_y, Q1, Q2, phasexlist[j + 1], phaseylist[j + 1], beam_direction, accel)
                if with_ac_calc:
                    [fwqw, bpms] = algorithms.helper.getFreeCoupling(Q1f, Q2f, Q1, Q2, fwqw, MADTwiss, bpms)
            else:
                raise ValueError('Number of monitors for coupling analysis (option -n) should be 1 or 2.')
            fwqw['DPP'] = dpop
    
    return Q1, Q2
# END phase_and_beta_for_non_zero_dpp --------------------------------------------------------------

def calculate_getsextupoles(twiss_model_file, TBTana, MADTwiss, files_dict, ListOfZeroDPPX, Q1f, phasexlist):
    print "Calculating getsextupoles"
    # fsex1200 and 2100 is never used and ouputname is wrong(vimaier)
    #     fsex1200 = open(outputpath+'getsex1200.out','w')
    #     fsex1200.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    #     fsex2100 = open(outputpath+'getsex1200.out','w')
    #     fsex2100.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    
    #-> 1) f3000 line (-2,0)
    #-> 2) f1200 line  (2,0)
    #-> 3) f2100 line  (0,0)
    # global stuff
    # 1)
    htot, afactor, pfactor = algorithms.helper.Getsextupole(MADTwiss, ListOfZeroDPPX, phasexlist[0], Q1f, 3, 0)
    filename = 'getsex3000.out'
    files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
    tfs_file = files_dict[filename]
    tfs_file.add_descriptor("f2h_factor", "%le", str(afactor))
    tfs_file.add_descriptor("p_f2h_factor", "%le", str(pfactor))
    tfs_file.add_column_names(["NAME", "S", "AMP_20", "AMP_20std", "PHASE_20", "PHASE_20std", "f3000", "f3000std", "phase_f_3000", "phase_f_3000std", "h3000", "h3000_std", "phase_h_3000", "phase_h_3000_std"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    for bpm in htot:
        li = htot[bpm]
        list_row_entries = [li[0], li[1], li[2], li[3], li[4], li[5], li[6], li[7], li[8], li[9], li[10], li[11], li[12], li[13]]
        tfs_file.add_table_row(list_row_entries)
# END calculate_getsextupoles ----------------------------------------------------------------------


def calculate_chiterms(twiss_model_file, accel, MADTwiss, files_dict, ListOfZeroDPPX, ListOfZeroDPPY):
#-> 1) chi3000
#-> 2) chi1010
#-> 2) chi4000
    print "Calculating chiterms"
    # 1) chi3000
    filename = 'getchi3000.out'
    files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
    tfs_file = files_dict[filename]
    tfs_file.add_column_names(["NAME", "S", "S1", "S2", "X3000", "X3000i", "X3000r", "X3000RMS", "X3000PHASE", "X3000PHASERMS", "X3000M", "X3000Mi", "X3000Mr", "X3000MPHASE"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    files = [ListOfZeroDPPX, ListOfZeroDPPY]
    name = 'chi3000'
    plane = 'H'
    [dbpms, POS, XItot, XIMODEL] = algorithms.helper.getChiTerms(MADTwiss, files, plane, name, ListOfZeroDPPX, ListOfZeroDPPY)
    for i in range(0, len(dbpms) - 2):
        bn = string.upper(dbpms[i][1])
        list_row_entries = ['"' + bn + '"', POS[0][i], POS[1][i], POS[2][i], XItot[0][i], XItot[1][i], XItot[2][i], XItot[3][i], XItot[4][i], XItot[5][i], XIMODEL[0][i], XIMODEL[1][i], XIMODEL[2][i], XIMODEL[3][i]]
        tfs_file.add_table_row(list_row_entries)
    
    # 2) chi1010
    if accel != 'SPS':
        filename = 'getchi1010.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_column_names(["NAME", "S", "X1010", "X1010RMS", "X1010PHASE", "X1010PHASERMS", "X1010M", "X1010MPHASE"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        files = [ListOfZeroDPPX, ListOfZeroDPPY]
        name = 'chi1010'
        plane = 'H'
        [dbpms, XItot] = algorithms.helper.getchi1010(MADTwiss, files, plane, name, ListOfZeroDPPX, ListOfZeroDPPY)
        for i in range(len(dbpms) - 2):
            bn = string.upper(dbpms[i][1])
            bns = dbpms[i][0]
            list_row_entries = ['"' + bn + '"', bns, XItot[0][i], XItot[1][i], XItot[2][i], XItot[3][i], '0', '0']
            tfs_file.add_table_row(list_row_entries)
    
    # 1) chi4000
    filename = 'getchi4000.out'
    files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
    tfs_file = files_dict[filename]
    tfs_file.add_column_names(["NAME", "S", "S1", "S2", "X4000", "X4000i", "X4000r", "X4000RMS", "X4000PHASE", "X4000PHASERMS", "X4000M", "X4000Mi", "X4000Mr", "X4000MPHASE"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
# END calculate_chiterms ---------------------------------------------------------------------------


def calculate_kick(twiss_model_file, lhcphase, beam_direction, MADTwiss, with_ac_calc, MADTwiss_ac, files_dict, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfZeroDPPY, ListOfNonZeroDPPY, Q1, Q2, Q1f, Q2f, acphasex_ac2bpmac, acphasey_ac2bpmac, betax_ratio, betay_ratio, betaxf_ratio, betayf_ratio):
    print "Calculating kick"
    files = [ListOfZeroDPPX + ListOfNonZeroDPPX, ListOfZeroDPPY + ListOfNonZeroDPPY]
    filename = 'getkick.out'
    files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
    tfs_file = files_dict[filename]
    tfs_file.add_descriptor("RescalingFactor_for_X", "%le", str(betax_ratio))
    tfs_file.add_descriptor("RescalingFactor_for_Y", "%le", str(betay_ratio))
    tfs_file.add_column_names(["DPP", "QX", "QXRMS", "QY", "QYRMS", "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX", "2JXSTD", "2JY", "2JYSTD", "sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES", "2JXSTDRES", "2JYRES", "2JYSTDRES"])
    tfs_file.add_column_datatypes(["%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    [invarianceJx, invarianceJy, tune, tuneRMS, dpp] = algorithms.helper.getkick(files, MADTwiss)
    for i in range(0, len(dpp)):
        list_row_entries = [dpp[i], tune[0][i], tuneRMS[0][i], tune[1][i], tuneRMS[1][i], invarianceJx[i][0], invarianceJx[i][1], invarianceJy[i][0], invarianceJy[i][1], (invarianceJx[i][0] ** 2), (2 * invarianceJx[i][0] * invarianceJx[i][1]), (invarianceJy[i][0] ** 2), (2 * invarianceJy[i][0] * invarianceJy[i][1]), (invarianceJx[i][0] / math.sqrt(betax_ratio)), (invarianceJx[i][1] / math.sqrt(betax_ratio)), (invarianceJy[i][0] / math.sqrt(betay_ratio)), (invarianceJy[i][1] / math.sqrt(betay_ratio)), (invarianceJx[i][0] ** 2 / betax_ratio), (2 * invarianceJx[i][0] * invarianceJx[i][1] / betax_ratio), (invarianceJy[i][0] ** 2 / betay_ratio), (2 * invarianceJy[i][0] * invarianceJy[i][1] / betay_ratio)]
        tfs_file.add_table_row(list_row_entries)
    
    if with_ac_calc:
        files = [ListOfZeroDPPX + ListOfNonZeroDPPX, ListOfZeroDPPY + ListOfNonZeroDPPY]
        filename = 'getkickac.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_descriptor("RescalingFactor_for_X", "%le", str(betaxf_ratio))
        tfs_file.add_descriptor("RescalingFactor_for_Y", "%le", str(betayf_ratio))
        tfs_file.add_column_names(["DPP", "QX", "QXRMS", "QY", "QYRMS", "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX", "2JXSTD", "2JY", "2JYSTD", "sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES", "2JXSTDRES", "2JYRES", "2JYSTDRES"])
        tfs_file.add_column_datatypes(["%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        [invarianceJx, invarianceJy, tune, tuneRMS, dpp] = algorithms.helper.getkickac(MADTwiss_ac, files, Q1, Q2, Q1f, Q2f, acphasex_ac2bpmac, acphasey_ac2bpmac, beam_direction, lhcphase)
        for i in range(0, len(dpp)):
            list_row_entries = [dpp[i], tune[0][i], tuneRMS[0][i], tune[1][i], tuneRMS[1][i], invarianceJx[i][0], invarianceJx[i][1], invarianceJy[i][0], invarianceJy[i][1], (invarianceJx[i][0] ** 2), (2 * invarianceJx[i][0] * invarianceJx[i][1]), (invarianceJy[i][0] ** 2), (2 * invarianceJy[i][0] * invarianceJy[i][1]), (invarianceJx[i][0] / math.sqrt(betax_ratio)), (invarianceJx[i][1] / math.sqrt(betax_ratio)), (invarianceJy[i][0] / math.sqrt(betay_ratio)), (invarianceJy[i][1] / math.sqrt(betay_ratio)), (invarianceJx[i][0] ** 2 / betax_ratio), (2 * invarianceJx[i][0] * invarianceJx[i][1] / betax_ratio), (invarianceJy[i][0] ** 2 / betay_ratio), (2 * invarianceJy[i][0] * invarianceJy[i][1] / betay_ratio)]
            tfs_file.add_table_row(list_row_entries)
# END calculate_kick -------------------------------------------------------------------------------


#===================================================================================================
# main invocation
#===================================================================================================
if __name__=="__main__":
    options,args=parse_args()
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
         higher_order=options.higher)

