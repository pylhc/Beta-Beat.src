r'''
.. module: drive.phase_correction

Created on 08/08/16

:author: Lukas Malina  

Phase_correction corrects phases in DRIVE output disturbed by the effect of oscillation damping (W. Guo @ IPAC16). 
It rewrites _linx and _liny files unless the suffix is specified (for tests). 
The correction performs better (much better), when the average tune is imposed in GetLLM. 
As an option the dpp value can be added, instead of manually written through GUI.
'''

import __init__  # @UnusedImport
from optparse import OptionParser
from Python_Classes4MAD import metaclass
import time
import numpy
import multiprocessing
from Utilities import tfs_file_writer

HEADERS_X = ["NAME", "S", "BINDEX", "SLABEL", "TUNEX", "NOISE", "PK2PK", "CO", "CORMS", "AMPX", "MUX", "AVG_MUX", "AMP01", "PHASE01", "AMP_20", "PHASE_20", "AMP02", "PHASE02", "AMP_30", "PHASE_30", "AMP_1_1", "PHASE_1_1", "AMP2_2", "PHASE2_2", "AMP0_2", "PHASE0_2", "AMP1_2", "PHASE1_2", "AMP_13", "PHASE_13", "AMP12", "PHASE12", "AMP_21", "PHASE_21", "AMP11", "PHASE11", "AMP20", "PHASE20", "AMP_1_2", "PHASE_1_2", "AMP30", "PHASE30", "NATTUNEX", "NATAMPX", "MUXCOR","AVG_MUXCOR"]
HEADERS_Y = ["NAME", "S", "BINDEX", "SLABEL", "TUNEY", "NOISE", "PK2PK", "CO", "CORMS", "AMPY", "MUY", "AVG_MUY", "AMP10", "PHASE10", "AMP_1_1", "PHASE_1_1", "AMP_20", "PHASE_20", "AMP1_1", "PHASE1_1", "AMP0_2", "PHASE0_2", "AMP0_3", "PHASE0_3", "AMP_11", "PHASE_11", "AMP21", "PHASE21", "AMP_13", "PHASE_13", "AMP11", "PHASE11", "AMP_12", "PHASE_12", "NATTUNEY", "NATAMPY", "MUYCOR","AVG_MUYCOR"]
MAX_NO_PROC = 16

def _parse_args():
    parser = OptionParser()
    parser.add_option("-f", "--files",
                    help="List of comma-separated sdds files from analysis",
                    metavar="FILES", dest="files")
    parser.add_option("-s", "--suffix",
                    help="Suffix to be put after linx/liny: creates a new file not used by GetLLM",
                    metavar="SUFFIX", default="", dest="suffix")
    parser.add_option("-p", "--dpp",
                    help="Edited dpoverp which is then used by GetLLM",
                    metavar="DPP", default=0.0, dest="dpp")
    options, _ = parser.parse_args()

    return options.files, options.suffix, options.dpp

def phase_correction(files, suffix,dpp):
    file_list = [file_name.strip() for file_name in files.strip("\"").split(",")]
    nProc=numpy.min([len(file_list),MAX_NO_PROC])
    process_pool = multiprocessing.Pool(processes=nProc)
    for file_name in file_list:
        process_pool.apply_async(_one_file_correction,(file_name,suffix,dpp))
    process_pool.close()
    process_pool.join()    

def _one_file_correction(file_name,suffix,dpp):
    linx,liny = _parse_twiss_files(file_name)
    _single_file_phasex(linx)
    _single_file_phasey(liny)
    linx_file = _create_linx_file(linx, suffix, dpp)
    liny_file = _create_liny_file(liny, suffix, dpp)
    for bpm in linx.NAME:
        _write_single_line_x(linx_file,linx,bpm)
    linx_file.write_to_file()
    for bpm in liny.NAME:
        _write_single_line_y(liny_file,liny,bpm)
    liny_file.write_to_file()

def _parse_twiss_files(file_name):
    startturn, endturn = _get_turns(file_name)
    twx = metaclass.twiss(file_name + "_linx")
    twy = metaclass.twiss(file_name + "_liny")
    twx.keys.append("SDDS")
    twy.keys.append("SDDS")
    empty_list = numpy.empty([len(twx.NAME),endturn-startturn])
    setattr(twx, "SDDS", empty_list)
    empty_list = numpy.empty([len(twy.NAME),endturn-startturn])
    setattr(twy, "SDDS", empty_list)
    update_twisses_with_sdds(file_name, twx, twy)
    return twx,twy

def update_twisses_with_sdds(file_path, twx, twy):
    column = [getattr(twx, "SDDS"), getattr(twy, "SDDS")]
    startturn, endturn = _get_turns(file_path)
    with open(file_path, "r") as filesdds:
        for line in filesdds:
            if line.startswith("#"):
                continue
            list_splitted_line_values = line.split()
            number_of_columns = len(list_splitted_line_values)
            if (number_of_columns < 3):
                print "not a valid data line: " + str(line)
                continue
            plane = int(list_splitted_line_values[0])
            bpm_name = list_splitted_line_values[1]
            if plane == 0:
                bpm_index = twx.indx[bpm_name]
            else:
                bpm_index = twy.indx[bpm_name]
            column[plane][bpm_index] = numpy.array(\
                    list_splitted_line_values[3 + startturn:3 + endturn], \
                    dtype=numpy.float64) - numpy.average(numpy.array(\
                    list_splitted_line_values[3:], \
                    dtype=numpy.float64))
        setattr(twx, "SDDS", column[0])
        setattr(twy, "SDDS", column[1])

def _get_turns(file_path):
    f = open(file_path[:(file_path.rindex("/")+1)] + "DrivingTerms", "r")
    for line in f:
        l=line.split()
    return int(l[1]),int(l[2])
     
def _single_file_phasex(tw):
    damping=_get_daming_time(tw.SDDS)
    damp=numpy.average(damping)
    print "Damping factor X: {0:2.2e} +- {1:2.2e}".format(damp,numpy.std(damping))
    q1=tw.Q1 * 2 * numpy.pi
    muxcor = []
    avgmuxcor = []
    for bpm_name in tw.NAME:
        ind=tw.indx[bpm_name]
        muxcor.append(_get_phase_correction(tw.PK2PK[ind]/2, damp, tw.SDDS[ind], tw.TUNEX[ind] * 2 * numpy.pi, tw.MUX[ind] * 2 * numpy.pi))
        avgmuxcor.append(_get_phase_correction(tw.PK2PK[ind]/2, damp, tw.SDDS[ind], q1, tw.AVG_MUX[ind] * 2 * numpy.pi))
    tw.keys.append("MUXCOR")
    mux = getattr(tw, "MUX")
    tw.keys.append("AVG_MUXCOR")
    avgmux = getattr(tw, "AVG_MUX")
    setattr(tw, "MUXCOR", muxcor)
    setattr(tw, "MUX", mux+muxcor) 
    setattr(tw, "AVG_MUXCOR", avgmuxcor)
    setattr(tw, "AVG_MUX", avgmux+avgmuxcor)   


def _single_file_phasey(tw):
    damping=_get_daming_time(tw.SDDS)
    damp=numpy.average(damping)
    print "Damping factor Y: {0:2.2e} +- {1:2.2e}".format(damp,numpy.std(damping))
    q2=tw.Q2 * 2 * numpy.pi
    muycor = []
    avgmuycor = []
    for bpm_name in tw.NAME:
        ind=tw.indx[bpm_name]
        muycor.append(_get_phase_correction(tw.PK2PK[ind]/2, damp, tw.SDDS[ind], tw.TUNEY[ind] * 2 * numpy.pi, tw.MUY[ind] * 2 * numpy.pi))
        avgmuycor.append(_get_phase_correction(tw.PK2PK[ind]/2, damp, tw.SDDS[ind], q2, tw.AVG_MUY[ind] * 2 * numpy.pi))
    tw.keys.append("MUYCOR")
    muy = getattr(tw, "MUY")
    tw.keys.append("AVG_MUYCOR")
    avgmuy = getattr(tw, "AVG_MUY")
    setattr(tw, "MUYCOR", muycor)
    setattr(tw, "MUY", muy+muycor) 
    setattr(tw, "AVG_MUYCOR", avgmuycor)
    setattr(tw, "AVG_MUY", avgmuy+avgmuycor)   
   
def _get_daming_time(sdds):         #assumes average around zero
    damp = []
    for a in sdds:
        m,b=numpy.polyfit(numpy.linspace(1.0,float(len(a)),num=len(a)),numpy.maximum.accumulate(numpy.log(numpy.abs(a[::-1]))),1)
        damp.append(-m)
    return damp  

def _get_phase_correction(amp,damp,x,tune,phase):
    i=numpy.linspace(1.0,float(len(x)),num=len(x))
    e1=numpy.sum(numpy.exp(2*damp*i)*numpy.sin(2*(i*tune + phase)))* amp/2
    e2=numpy.sum(x*numpy.exp(damp*i)*numpy.sin(i*tune + phase))
    e3=numpy.sum(x*numpy.exp(damp*i)*numpy.cos(i*tune + phase))
    e4=numpy.sum(numpy.exp(2*damp*i)*numpy.cos(2*(i*tune + phase)))* amp
    mux_cor = (e1-e2) / ((e3-e4) * 2 * numpy.pi)
    return mux_cor
        
def _create_linx_file(tw, suffix,dpp):
    lin_outfile = tfs_file_writer.TfsFileWriter(tw.filename + suffix)
    lin_outfile.add_float_descriptor("DPP",dpp)
    lin_outfile.add_float_descriptor("Q1",tw.Q1)
    lin_outfile.add_float_descriptor("Q1RMS",tw.Q1RMS)
    if "NATQ1" in tw.keys:
        lin_outfile.add_float_descriptor("NATQ1",tw.NATQ1)  
    if "NATQ1RMS" in tw.keys:    
        lin_outfile.add_float_descriptor("NATQ1RMS",tw.NATQ1RMS)   
    lin_outfile.add_column_names(HEADERS_X)
    lin_outfile.add_column_datatypes(["%s"] + ["%le"] * (len(HEADERS_X) - 1))
    return lin_outfile

def _create_liny_file(tw, suffix, dpp):
    lin_outfile = tfs_file_writer.TfsFileWriter(tw.filename + suffix)
    lin_outfile.add_float_descriptor("DPP",dpp)
    lin_outfile.add_float_descriptor("Q2",tw.Q2)
    lin_outfile.add_float_descriptor("Q2RMS",tw.Q2RMS)
    if "NATQ2" in tw.keys:
        lin_outfile.add_float_descriptor("NATQ2",tw.NATQ2)  
    if "NATQ2RMS" in tw.keys:    
        lin_outfile.add_float_descriptor("NATQ2RMS",tw.NATQ2RMS)   
    lin_outfile.add_column_names(HEADERS_Y)
    lin_outfile.add_column_datatypes(["%s"] + ["%le"] * (len(HEADERS_Y) - 1))
    return lin_outfile

def _write_single_line_x(lin_outfile,tw,name):
    i = tw.indx[name]
    row = [tw.NAME[i], tw.S[i], tw.BINDEX[i], tw.SLABEL[i], tw.TUNEX[i], tw.NOISE[i], tw.PK2PK[i], tw.CO[i], tw.CORMS[i], tw.AMPX[i], tw.MUX[i], tw.AVG_MUX[i], tw.AMP01[i], tw.PHASE01[i], tw.AMP_20[i], tw.PHASE_20[i], tw.AMP02[i], tw.PHASE02[i], tw.AMP_30[i], tw.PHASE_30[i], tw.AMP_1_1[i], tw.PHASE_1_1[i], tw.AMP2_2[i], tw.PHASE2_2[i], tw.AMP0_2[i], tw.PHASE0_2[i], tw.AMP1_2[i], tw.PHASE1_2[i], tw.AMP_13[i], tw.PHASE_13[i], tw.AMP12[i], tw.PHASE12[i], tw.AMP_21[i], tw.PHASE_21[i], tw.AMP11[i], tw.PHASE11[i], tw.AMP20[i], tw.PHASE20[i], tw.AMP_1_2[i], tw.PHASE_1_2[i], tw.AMP30[i], tw.PHASE30[i], tw.NATTUNEX[i], tw.NATAMPX[i], tw.MUXCOR[i], tw.AVG_MUXCOR[i]]
    lin_outfile.add_table_row(row)

def _write_single_line_y(lin_outfile,tw,name):
    i = tw.indx[name]
    row = [tw.NAME[i], tw.S[i], tw.BINDEX[i], tw.SLABEL[i], tw.TUNEY[i], tw.NOISE[i], tw.PK2PK[i], tw.CO[i], tw.CORMS[i], tw.AMPY[i], tw.MUY[i], tw.AVG_MUY[i], tw.AMP10[i], tw.PHASE10[i], tw.AMP_1_1[i], tw.PHASE_1_1[i], tw.AMP_20[i], tw.PHASE_20[i], tw.AMP1_1[i], tw.PHASE1_1[i], tw.AMP0_2[i], tw.PHASE0_2[i], tw.AMP0_3[i], tw.PHASE0_3[i], tw.AMP_11[i], tw.PHASE_11[i], tw.AMP21[i], tw.PHASE21[i], tw.AMP_13[i], tw.PHASE_13[i], tw.AMP11[i], tw.PHASE11[i], tw.AMP_12[i], tw.PHASE_12[i], tw.NATTUNEY[i], tw.NATAMPY[i], tw.MUYCOR[i], tw.AVG_MUYCOR[i]]
    lin_outfile.add_table_row(row)

   
if __name__ == "__main__":
    timeStartGlobal = time.time()
    _files, _suffix, _dpp = _parse_args()
    phase_correction(_files, _suffix, float(_dpp))
    timeGlobal = time.time() - timeStartGlobal
    print "Duration:", timeGlobal, "s"
