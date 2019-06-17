'''
@first created by R. Westenberger
@updated by F. Carlier

@version: 0.0.1

GetLLM.algorithms.resonant_driving_terms.py stores helper functions for RDT calculations for GetLLM.
This module is not intended to be executed. It stores only functions.
'''

import sys

import utils.bpm
import helper
import numpy as np

from utils import logging_tools

LogRdt = logging_tools.get_logger(__name__)

exist_curve_fit = True

try:
    from scipy.optimize import curve_fit
except ImportError:
    exist_curve_fit = False


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


RDT_LIST = ['f1001H', 'f1010H', 'f0110V', 'f1010V',  #Quadrupolar
            'f3000H', 'f1200H', 'f1020H', 'f1002H',  #Normal Sextupolar
            'f0111V', 'f1020V', 'f0120V', 'f1011V', 
            'f0030V', 'f0012V', 'f0210V', 'f2010V',  #Skew Sextupolar
            'f1101H', 'f2010H', 'f1110H', 'f2001H', 
            'f4000H', 'f1300H', 'f2002H', 'f1120H',  #Normal Octupolar
            'f1102H', 'f2020H', 'f2020V', 'f2011V', 
            'f0220V', 'f0211V', 'f0040V', 'f0013V',
            ]

def determine_lines(rdt):
    r = list(rdt)
    j, k, l, m, plane = int(r[1]), int(r[2]), int(r[3]), int(r[4]), r[5]
    if plane == 'H':
        line = (1-j+k, m-l)
    elif plane == 'V':
        line = (k-j, 1-l+m)
    return line, plane


def calculate_RDTs(mad_twiss, getllm_d, twiss_d, phase_d, tune_d, files_dict, inv_x, inv_y):
    '''
    Calculates line RDT amplitudes and phases and fills the following TfsFiles:
        f3000_line.out ...

    :Parameters:
        'mad_twiss': ? 
            MAD model
        'getllm_d': _GetllmData (In-param, values will only be read)
            lhc_phase, accel, beam_direction and num_beams_for_coupling are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances.
    '''
    print "Calculating RDTs"
    
    #if sys.flags.debug:
    #    print "DEBUG is ON"
    #    LogRdt.setLevel( logging_tools.DEBUG - 1 )
    #else:
    #    print "DEBUG is OFF"
    #    LogRdt.setLevel(logging_tools.WARNING)



    """
    The rdt_set holds all RDTs which should be investigated with the parameters to call GetRDT()
    syntax is: rdt_set = [(plane, out_file, line), ...]
    with:
        plane in ["H", "V"]
        out_file in files_dict is the out file to write the data to (must be added to GetLLM.py)
        line in (int, int) is the corresponding line to the driving term
    """
    beam = getllm_d.beam_direction

    if exist_curve_fit:
        for rdt in RDT_LIST:
            line, plane = determine_lines(rdt)
            _process_RDT(mad_twiss, phase_d, twiss_d, (plane, files_dict[rdt+'_line.out'], files_dict[rdt+'.out'], line), inv_x, inv_y, rdt, beam)
    else:
        print 'Curve fit not imported.. RDTs skipped'


def _process_RDT(mad_twiss, phase_d, twiss_d, (plane, out_file, rdt_out_file, line), inv_x, inv_y, rdt, beam):
    '''
        'mad_twiss': twiss (tfs) with the model
        'twiss_d'  : contains lists of linx and liny files grouped with zeron and non-zero dpp
                    twiss_d.zero_dpp_x
                    twiss_d.non_zero_dpp_x
                    twiss_d.zero_dpp_y
                    twiss_d.non_zero_dpp_y
        'phase_d'  : class PhaseData in phase.py, contains lists of reconstructed phases
                    fille by get_phases in phase.py
                    interesting field is ph_x, which is a dictionary indexed by name of a BPM
         
    '''
    assert plane in ["H", "V"] # check user input plane
    
    strline = plane + str(line)
    LogRdt.debug("----------------------------------------------")
    LogRdt.debug("Processing line rdt %s line %s in plane %s",rdt,strline,plane)
    LogRdt.debug("----------------------------------------------")
    LogRdt.debug("outfile: %s",out_file.get_file_name())
    
    
    # get plane corresponding phase and twiss data
    linx_data = twiss_d.zero_dpp_x  # linx for all the input files
    liny_data = twiss_d.zero_dpp_y  # liny for all the input files
    
    if plane == "H":
        phase_data = phase_d.ph_x  # for description see class PhaseData in phase.py 
        list_zero_dpp = twiss_d.zero_dpp_x # linx for all the input files
    else:
        phase_data = phase_d.ph_y
        list_zero_dpp = twiss_d.zero_dpp_y
    
    # both planes need to be there because both phases advances are needed 
    dbpms = utils.bpm.intersect(twiss_d.zero_dpp_y+twiss_d.zero_dpp_x)
    dbpms = utils.bpm.model_intersect(dbpms, mad_twiss)
    bpm_positions, bpm_names = zip(*dbpms)


    # init out file
    out_file.add_column_names(["NAME", "S", "COUNT", "AMP", "EAMP", "PHASE", "EPHASE"])
    out_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])

    line_amplitudes = [] # stored 1 entry per BPM per file, one after another 
    line_amplitudes_err = []
    line_amp_bpmidx = [] #book keeping of BPM index so the 2 above can be assigned to the good BPM when one is bad
    rdt_phases_averaged = []
    rdt_phases_averaged_std = []

    use_line = False
    use_opposite_line = False   
    
    try:
        #if DEBUG:
        #    print("Looking for normal line (%s, %s)" % (line[0],line[1]))
        LogRdt.debug("Looking for normal line (%s, %s)",line[0],line[1])
        
        _, _ = _line_to_amp_and_phase_attr(line, list_zero_dpp[0])
        use_line = True
    except AttributeError:
        #print >> sys.stderr, "Line (%s, %s) not found! Trying opposite line ... !" % line
        LogRdt.info("Line (%s, %s) not found! Trying opposite line ... !", line[0],line[1])
    
    try:
        #if DEBUG:
        #    print("Looking for oposit line (%s, %s)" % (-line[0],-line[1]))
        LogRdt.debug("Looking for oposit line (%s, %s)",-line[0],-line[1])
        
        _, _ = _line_to_amp_and_phase_attr((-line[0],-line[1]), list_zero_dpp[0])
        use_opposite_line = True
    except AttributeError:
        #print >> sys.stderr, "Opposite line (%s, %s) not found!" % (-line[0],-line[1])
        LogRdt.info("Opposite line (%s, %s) also not found!", -line[0],-line[1] )
        
        
    if use_line or use_opposite_line:
        
        # 1st BPM in the list (normally first in the measurement)
        # Its phase will be subtracted from the measured phases
        # to get rif of shot to shot phase changes
        
        bpmRef = dbpms[0][1].upper()
        LogRdt.debug("%s: First BPM is %s",strline,bpmRef)
        LogRdt.debug("%s: First BPM spos %f", strline, bpm_positions[bpm_names.index(bpmRef)] )
        

        LogRdt.debug("%s: There is %d BPMs",strline,len(dbpms))
        for i in range(len(dbpms)):
            bpm1 = dbpms[i][1].upper()
            LogRdt.debug("bpm1: %s",bpm1)
            
        
        
          
        for i in range(len(dbpms)):
            LogRdt.debug("  ")
            LogRdt.debug(" -------- ------------------------------------------------ ")
            LogRdt.debug("  ")

            bpm1 = dbpms[i][1].upper()
            LogRdt.debug("bpm1: %s",bpm1)
            bpm1spos = bpm_positions[bpm_names.index(bpm1)]
            
            try:
                # these are: best_90degrees_bpm, best_90degrees_phase, best_90degrees_phase_std
                # filled by get_phases in phase.py
                bpm_pair_data = phase_data[bpm1][7], phase_data[bpm1][8], phase_data[bpm1][9]
            except KeyError:
                LogRdt.warn("%s: Could not find a BPM pair (%s, %s)! No data for pi/2 BPM",strline, plane, bpm1)
                #print >> sys.stderr, "Could not find a BPM pair (%s, %s)!" % (plane, bpm1)
                continue
            
            # BPM name closest to 90 degreees downstream 
            bpm2 = bpm_pair_data[0]
             
            if bpm2 not in bpm_names:
                # this happens when the pi/2 BPM was cleaned in the other plane 
                LogRdt.warn("%s: Could not find a BPM pair (%s, %s)! missing data for pi/2 BPM %s",
                            strline, plane, bpm1, bpm2)
                continue
                
            bpm2spos = bpm_positions[bpm_names.index(bpm2)]
            
            rdt_phases_per_bpm = []
            
            delta, edelta = bpm_pair_data[1:]
            if beam==-1:
                delta = -delta
            
            LogRdt.debug("%s: %s <--> %s : %f",strline,bpm1, bpm2, delta)
            LogRdt.debug("  ")
            amplist = []
            amp1list = []
            amp2list = []
            phase1list = []
            phase2list = []
            linePhaseslist = []
            #Loop over the linx/liny files in list_zero_dpp list
            # list_zero_dpp is linx/liny for the current plane, the same as one of linx_data,liny_data
            
            LogRdt.debug("Starting loop !! ")
            LogRdt.debug(len(list_zero_dpp))
            LogRdt.debug("Starting loop !! ")
            
            for j in range(0,len(list_zero_dpp)):
                
                LogRdt.debug("Loop j = %d ",j)

                idxlxbpm1 = linx_data[j].indx[bpm1]
                idxlybpm1 = liny_data[j].indx[bpm1]
                LogRdt.debug("Getting mux and muy from %d and %d", idxlxbpm1, idxlybpm1)
                # phase of tune line
                ph_H10 = getattr(linx_data[j], "MUX")[idxlxbpm1]%1
                ph_V01 = getattr(liny_data[j], "MUY")[idxlybpm1]%1
                # tunes
                q1 = getattr(linx_data[j], "TUNEX")[idxlxbpm1]%1
                q2 = getattr(liny_data[j], "TUNEY")[idxlybpm1]%1
                
                #q1 = tune_d.q1f
                #q2 = tune_d.q2f
                

                idxbpm1 = list_zero_dpp[j].indx[bpm1]
                idxbpm2 = list_zero_dpp[j].indx[bpm2]
                
                # Get amplitude and phase of the line from linx/liny file
                # it calculates average for all the BPMs present
                if use_line and use_opposite_line:
                    # amp_line, phase_line: these are lists for all BPMs 
                    amp_line_l, phase_line_l = _line_to_amp_and_phase_attr(line, list_zero_dpp[j])
                    amp_line_opp_l, phase_line_opp_l = _line_to_amp_and_phase_attr((-line[0],-line[1]), list_zero_dpp[j])
                    phase_line_opp_l = -phase_line_opp_l 
                    amp_line_l = (amp_line_l + amp_line_opp_l)/2.
                    phase_line_l = (phase_line_l + phase_line_opp_l)/2.
                    
                    amp1_line_no = amp_line_l[idxbpm1]
                    amp1_line_op = amp_line_opp_l[idxbpm1]
                    amp2_line_no = amp_line_l[idxbpm2]
                    amp2_line_op = amp_line_opp_l[idxbpm2]
                    
                    amp1 = (amp1_line_no + amp1_line_op)/2.
                    amp2 = (amp2_line_no + amp2_line_op)/2.
                    

                    phase1_line_no = phase_line_l[idxbpm1]
                    phase1_line_op = phase_line_opp_l[idxbpm1]
                    phase2_line_no = phase_line_l[idxbpm2]
                    phase2_line_op = phase_line_opp_l[idxbpm2]
                    
                    phase1 = (phase1_line_no + phase1_line_op)/2.
                    phase2 = (phase2_line_no + phase2_line_op)/2.
                    
                    
                    
                elif use_line and not use_opposite_line:
                    amp_line_l, phase_line_l = _line_to_amp_and_phase_attr(line, list_zero_dpp[j])
                    amp1   = amp_line_l[idxbpm1]
                    amp2   = amp_line_l[idxbpm2]
                    phase1 = phase_line_l[idxbpm1]
                    phase2 = phase_line_l[idxbpm2]
                    
                elif use_opposite_line and not use_line:
                    amp_line_l, phase_line_l = _line_to_amp_and_phase_attr((-line[0],-line[1]), list_zero_dpp[j])
                    amp1   = amp_line_l[idxbpm1]
                    amp2   = amp_line_l[idxbpm2]
                    phase1 = -phase_line_l[idxbpm1]
                    phase2 = -phase_line_l[idxbpm2]

                
                phase1 = phase1%1
                phase2 = phase2%1
                
                #Get amplitude and pha
                
                        
                # amp1 = amp_line_l[idxbpm1]
                # amp2 = amp_line_l[idxbpm2]
                # phase1 = phase_line_l[idxbpm1]
                # phase2 = phase_line_l[idxbpm2]

                
                amp1list.append(amp1)
                amp2list.append(amp2)
                phase1list.append(phase1)
                phase2list.append(phase2)
                 
                if bpm2spos < bpm1spos:
                    LogRdt.debug("< %f %f tune %f ", phase1, phase2,q1)
                    phase2 = (phase2 + line[0]*q1 + line[1]*q2)%1
                    LogRdt.debug("> %f %f tune %f ", phase1, phase2,q1)

                #LogRdt.debug("%s:    delta Phi = %f  ",strline, (phase2-phase1)%1 )
                #LogRdt.debug("%s:    BPM1 Amp = %f Phase = %f ",strline, amp1, phase1)
                #LogRdt.debug("%s:    BPM2 Amp = %f Phase = %f ",strline, amp2, phase2)
                
                
                if amp1 == 0 or amp2 == 0:
                    line_amp, line_amp_e, line_phase, line_phase_e = 0,0,0,0
                else:
                    line_amp, line_phase, line_amp_e, line_phase_e = helper.ComplexSecondaryLineExtended(delta,edelta, amp1, amp2, phase1, phase2)
                
                LogRdt.debug("%s:   Line Amp = %f +/- %f   Phase = %f +/- %f  deltaMU = %f",
                                strline,line_amp, line_amp_e, line_phase, line_phase_e, delta)
                linePhaseslist.append(line_phase)
                
                
                LogRdt.debug("%s:   Setting _line results for BPM %s at %f",strline, bpm1, dbpms[i][0])
                out_file.add_table_row([bpm1, dbpms[i][0], len(list_zero_dpp), line_amp, line_amp_e, line_phase, line_phase_e])
                
                # these 2 lists agragate all the amps: for all BPMs and all measurements 
                line_amplitudes.append(line_amp)
                line_amplitudes_err.append(line_amp_e)
                line_amp_bpmidx.append(i)
                
                amplist.append(line_amp)    
                    

                
                if line_amp != 0:
                    rdt_phases_per_bpm.append(calculate_rdt_phases(rdt, line_phase, ph_H10, ph_V01)%1)
            
            #rdt_phases_per_bpm = [0.55, 0.55, 0.55, 0.55]
            
            
            LogRdt.debug("amps 1  for all the files")
            LogRdt.debug(amp1list)
            LogRdt.debug("amps 2  for all the files")
            LogRdt.debug(amp2list)
            LogRdt.debug("ave amlitudes  for all the files")
            LogRdt.debug(amplist)
            

            LogRdt.debug("phases 1 for all the files")
            LogRdt.debug(phase1list)
            LogRdt.debug("phases 2 for all the files")
            LogRdt.debug(phase2list)
            
            LogRdt.debug("line phases  for all the files")
            LogRdt.debug(linePhaseslist)
            
            
            np.set_printoptions(formatter={'float': '{: 0.2f}'.format},linewidth=750)

            if rdt_phases_per_bpm:
                # was simple average
                #rdt_phases_averaged.append(np.average(np.array(rdt_phases_per_bpm)))
                #rdt_phases_averaged_std.append(np.std(np.array(rdt_phases_per_bpm)))
                
                # implemented circular average
                
                phases = np.array(rdt_phases_per_bpm)*2*np.pi

                # skowron October 2018

                LogRdt.debug("phases for all the files")
                LogRdt.debug(phases)

                rdt_sinphases = np.sin( phases )
                rdt_cosphases = np.cos( phases )

                LogRdt.debug("phases sin cos")
                LogRdt.debug(rdt_sinphases)
                LogRdt.debug(rdt_cosphases)
                
                sinave = np.average(rdt_sinphases)
                cosave = np.average(rdt_cosphases)

                LogRdt.debug("phases sin cos averages %5.3f %5.3f",sinave, cosave)
                
                phase = np.arctan2(sinave,cosave)
                if phase < 0:
                    phase = phase + 2*np.pi
                
                # error propagation for atan2 if we had phase errors from Sussix/HarPy            
                # D(phi) =  D(atan(x/y)) = ( D(x)*y + D(y)*x )/ (x^2 + y^2)
                sinstd = np.std(rdt_sinphases)
                cosstd = np.std(rdt_cosphases)
                
                phaseerr =             np.abs(sinstd*cosave)
                phaseerr =  phaseerr + np.abs(cosstd*sinave)
                phaseerr =  phaseerr / (sinave*sinave + cosave*cosave)
                
                LogRdt.debug("phase = %f +/- %f [rad]",phase/(2*np.pi),phaseerr/(2*np.pi))
                
                #if DEBUG:
                #s    print("phase = %f +/- %f [rad]"%(phase/(2*np.pi),phaseerr/(2*np.pi)))
                
                rdt_phases_averaged.append(phase/(2*np.pi))
                rdt_phases_averaged_std.append(phaseerr/(2*np.pi))

                
                
            else:
                rdt_phases_averaged.append(0.0)
                rdt_phases_averaged_std.append(0.0)
            
            np.set_printoptions(formatter=None)

    else:
        #print >> sys.stderr, "Could not find line for %s !" %rdt
        LogRdt.warn("Could not find line for %s !",rdt)
        
    # init out file
    rdt_out_file.add_column_names(["NAME", "S", "COUNT", "AMP", "EAMP", "PHASE", "PHASE_STD", "REAL", "IMAG"])
    rdt_out_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le",  "%le", "%le", "%le"])

    rdt_angles = np.mod(np.array(rdt_phases_averaged), np.ones(len(rdt_phases_averaged)))
    real_part = np.cos(2*np.pi*rdt_angles)
    imag_part = np.sin(2*np.pi*rdt_angles)
    inv_x = np.array(inv_x)
    inv_y = np.array(inv_y)
    
    
    for k in range(len(line_amplitudes)/len(list_zero_dpp)):
        num_meas = len(list_zero_dpp)
        bpm_rdt_data = np.array(line_amplitudes[k*num_meas:(k+1)*num_meas])
        #bpm_idx_check = np.std(line_amp_bpmidx[k*num_meas:(k+1)*num_meas])
        # if bpm_idx_check > 1e-6: error 
        bpm_idx = line_amp_bpmidx[k*num_meas]
        bpm_name = dbpms[bpm_idx][1].upper()
        #LogRdt.debug("%s: %s in  line amps :",strline,bpm_name)
        #LogRdt.debug(bpm_rdt_data)
        res, res_err = do_fitting(bpm_rdt_data, inv_x, inv_y, rdt, plane)
        LogRdt.debug("%s: %s fit line amp : %e",strline,bpm_name,res[0])
        
        
        #print "skowron EAMP", res_err[0]
        if np.isinf(res_err[0]):
            #print "skowron EAMP is inf", res_err[0]
            res_err[0] = 0.0
        LogRdt.debug("%s:   Setting results for BPM %s at %f",strline, bpm_name, dbpms[bpm_idx][0])    
        rdt_out_file.add_table_row([bpm_name, dbpms[bpm_idx][0], len(list_zero_dpp), res[0], res_err[0], rdt_angles[k], rdt_phases_averaged_std[k], res[0]*real_part[k], res[0]*imag_part[k]])


def calculate_rdt_phases(rdt, line_phase, ph_H10, ph_V01):
    r = list(rdt)
    j, k, l, m, plane = int(r[1]), int(r[2]), int(r[3]), int(r[4]), r[5]
    if plane == 'H':
        rdt_phase = line_phase - (k-j+1)*ph_H10 - (m-l)*ph_V01 + 0.25
    elif plane == 'V':
        rdt_phase = line_phase - (k-j)*ph_H10 - (m-l+1)*ph_V01 + 0.25
    return rdt_phase 


def rdt_function_gen(rdt, plane):
    '''
    Note that the factor 2 in 2*j*f_jklm*.... is absent due to the normalization with the main line. 
    The main line has an amplitude of sqrt(2J*beta)/2
    '''
    r = list(rdt)
    j, k, l, m, plane = int(r[1]), int(r[2]), int(r[3]), int(r[4]), r[5]
    if plane == 'H':
        def rdt_function(x, f):
            return 2 * j * f * x[0]**((j+k-2)/2.) * x[1]**((l+m)/2.)
    elif plane == 'V':
        def rdt_function(x, f):
            return 2 * l * f * x[0]**((j+k)/2.) * x[1]**((l+m-2)/2.)
    return rdt_function


def do_fitting(bpm_rdt_data, kick_x, kick_y, rdt, plane):
    mask = bpm_rdt_data!=0
    
    if sum(mask)!=0:
        func = rdt_function_gen(rdt, plane)
        kick_data = np.vstack((np.transpose(kick_x[mask])[0]**2, np.transpose(kick_y[mask])[0]**2))
        popt, pcov = curve_fit(func, kick_data, bpm_rdt_data[mask])
        perr = np.sqrt(np.diag(pcov))
    else:
        popt, perr = [0], [0]
    return popt, perr


def _line_to_amp_and_phase_attr(line, zero_dpp):
    '''
    Returns for the line amplitude and phase columns (for all BPMS)
    'line': list of 2 integers defining a line is spectra
    'zero_dpp': is a linx or liny TFS file 
    '''
    #To turn input line (-1,2) to (zero_dpp.AMP_12, zero_dpp.PHASE_12).
    line = (str(line[0])+str(line[1])).replace("-", "_")
    return (getattr(zero_dpp, "AMP"+line), getattr(zero_dpp, "PHASE"+line))

def _adjustTo01range(ph):
    '''
    Returns phase in 0-1 range
    'ph': a phase in 2pi units 
    '''
    
    while ph < 0:
        ph += 1
    
    while ph > 1:
        ph -= 1
    
    return ph
