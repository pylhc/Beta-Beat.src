'''
Created sometime 2009-2011

:maintainer: Yngve Inntjore Levinsen

:author: Glenn Vanbavinckhove, Yngve Inntjore Levinsen

:version: 3.1.3

Getsuper calculates with the help of GetLLM chromatic beta functions.
The Montague W and Phi functions are output.

What getsuper essentially does is run GetLLM on files with different dp/p and then afterwards
interpolate the results to see how the functions vary with dp/p.

To run getsuper.py you need several source files(sdds files) with different DPP(delta_p/p).
At least one of the source files must have DPP=0.0 .
Hint: You can change DPP in the GUI application in the Analysis panel. Change the entries in the
column 'dp/d' in the table at the top.

Further you need AC dipole in your model. If your twissfile is 'Twiss.dat' you need to provide the
file 'Twiss_ac.dat' in the same directory.

Usage cmd line::

    Usage: getsuper.py [options] sdds-file1 [sdds-file2 ...]

    Options:
      -h, --help            show this help message and exit
      -f TwissFile, --files=TwissFile
                            Files from analysis, separated by comma
      --madxbin=<path>      Path to mad-x binary
      --twissfile=/path/to/twiss.dat
                            Twiss file to use
      -o <path>, --output=<path>
                            Output path, where to store the results
      -b <path>, --beta=<path>
                            Path to Beat-Beat.src folder
      -t ALGORITHM, --algorithm=ALGORITHM
                            Which algorithm to use (SUSSIX/SVD)
      -a ACCEL, --accel=ACCEL
                            Which accelerator: LHCB1 LHCB2 SPS RHIC
      -d <deltapScalingFactor>, --deltapScalingFactor=<deltapScalingFactor>
                            Scaling factor for deltap, remember final value must
                            be in MAD units

Usage in another Python module::

    import GetLLM.getsuper
    ...
    GetLLM.getsuper.main(my_src_files_list, path_to_twiss_file)

'''

import optparse
import os
import sys
import shutil
import subprocess
import math
import re

import __init__  # @UnusedImport init will include paths
import Python_Classes4MAD.metaclass as metaclass
import Python_Classes4MAD.linreg as linreg
import Utilities.bpm
import Utilities.iotools
import superutils

#===================================================================================================
# _parse_args()-function
#===================================================================================================
def _parse_args():
    ''' Parses arguments from command line. '''
    usage = "usage: %prog [options] sdds-file1 [sdds-file2 ...]"
    parser = optparse.OptionParser(usage)
    # general
    parser.add_option("-f", "--files",
        help="Files from analysis, separated by comma",
        metavar="TwissFile", default="", dest="files")
    parser.add_option("--madxbin",
            help="Path to mad-x binary",
            metavar="<path>", default="madx", dest="madx")
    parser.add_option("--twissfile",
            help="Twiss file to use",
            metavar="/path/to/twiss.dat", default="", dest="twissfile")
    parser.add_option("-o", "--output",
            help="Output path, where to store the results",
            metavar="<path>", default="./", dest="output")
    # By default we take the path from where getsuper.py is ran from..
    parser.add_option("-b", "--beta",
            help="Path to Beat-Beat.src folder",
            metavar="<path>", dest="brc",
            default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) )
    parser.add_option("-t", "--algorithm",
            help="Which algorithm to use (SUSSIX/SVD)",
            metavar="ALGORITHM", default="SUSSIX", dest="technique")
    parser.add_option("-a", "--accel",
            help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
            metavar="ACCEL", default="LHCB1",dest="accel")

    parser.add_option("-d", "--deltapScalingFactor",
            help="Scaling factor for deltap, remember final value must be in MAD units",
            metavar="<deltapScalingFactor>", default=1.0, type=float,dest="deltap_scaling_factor")

    return parser.parse_args()


#===================================================================================================
# main()-function
#===================================================================================================
def main(
        source_files,
        twissfile,
        path_to_madx = "madx",
        output_path = "./",
        path_to_beta_beat = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
        technique = "SUSSIX",
        accel = "LHCB1",
        deltap_scaling_factor = 1.0
        ):
    """ getsuper main function 
    
    :param list source_files: list of strings with file_paths to TFS files
    :param string twissfile: path to TFS model file
    :param string path_to_madx: path to madx binary
    :param string output_path: Path to store created files
    :param string path_to_beta_beat: root of Beta-Beat.src
    :param string technique: Used Turn-by-turn data analysis algorithm: 'SUSSIX' or 'SVD'
    :param string accel: Accelerator(LHCB1 LHCB2 SPS or RHIC)
    :param float deltap_scaling_factor: Scaling factor for deltap, remember final value must be in MAD units
    """
    _InputData.static_init(
            source_files, path_to_madx, twissfile, output_path, path_to_beta_beat, technique, accel, deltap_scaling_factor
    )


    files_dict = {} # dpp --> files with corresponding dpp

    for f_name in _InputData.source_files:

        if f_name.endswith('.gz'):
            datax = metaclass.twiss(f_name.replace(".gz", "_linx.gz"))
            datay = metaclass.twiss(f_name.replace(".gz", "_liny.gz"))
        else:
            datax = metaclass.twiss(f_name+"_linx")
            datay = metaclass.twiss(f_name+"_liny")

        dppx = datax.DPP * _InputData.deltap_scaling_factor      # Quick hack to be able to use old files with bad dpp input
        dppy = datay.DPP * _InputData.deltap_scaling_factor

        if dppx != dppy:
            raise ValueError("Discrepancy between horizontal and vertical => "+str(dppx)+" "+str(dppy))
        else:
            dpp = dppx/1.0

        if dpp not in files_dict:
            print "Adding dpp", dpp
            files_dict[dpp] = [f_name]
        else:
            files_dict[dpp].append(f_name)

    if 0 not in files_dict:
        raise ValueError("NO DPP=0.0. Provide at least one source file with DPP=0.0.\n"+
                         "In GUI you can change DPP in the table at the top of the analyses panel."+
                         "Click in the column 'delta p /p' on a value to change it.")

    _madcreator(files_dict.keys(), files_dict)
    print "All models are created"
    for dpp in files_dict:
        files = files_dict[dpp]
        _rungetllm(_join_with_output_path("twiss_"+str(dpp)+".dat"), files, dpp)
    # The GUI wants the default files to have the names without _0.0
    _copy_default_outfiles()

    ##adding data
    betalistx = {}
    betalisty = {}
    couplelist = {}
    betalistxf = {}
    betalistyf = {}
    couplelistf = {}

    listx = []
    listxf = []
    listy = []
    listyf = []
    listc = []
    listcf = []

    try:
        metaclass.twiss(_join_with_output_path('getbetax_free_'+str(dpp)+'.out'))
        freeswitch = 1
    except IOError:
        print "WARNING: Could not open", _join_with_output_path('getbetax_free_'+str(dpp)+'.out')
        freeswitch = 0

    for dpp in files_dict.keys():
        print "Loading driven data for ", dpp
        betx = metaclass.twiss(_join_with_output_path('getbetax_'+str(dpp)+'.out'))
        bety = metaclass.twiss(_join_with_output_path('getbetay_'+str(dpp)+'.out'))
        couple = metaclass.twiss(_join_with_output_path('getcouple_'+str(dpp)+'.out'))
        #couple=twiss(options.output+'/getbetay_'+str(dpp)+'.out')
        betalistx[dpp] = betx
        betalisty[dpp] = bety
        couplelist[dpp] = couple

        if float(dpp)==0.0:
            zerobx = betx
            zeroby = bety

        listx.append(betx)
        listy.append(bety)
        listc.append(couple)

        modeld = metaclass.twiss(_InputData.twissfile)

        #try:
        if freeswitch == 1:
            print "Loading free data"
            print 'getbetax_free_'+str(dpp)+'.out'
            betxf = metaclass.twiss(_join_with_output_path('getbetax_free_'+str(dpp)+'.out'))
            betyf = metaclass.twiss(_join_with_output_path('getbetay_free_'+str(dpp)+'.out'))
            couplef = metaclass.twiss(_join_with_output_path('getcouple_free_'+str(dpp)+'.out'))
            betalistxf[dpp] = betxf
            betalistyf[dpp] = betyf
            couplelistf[dpp] = couplef
            listxf.append(betxf)
            listyf.append(betyf)
            listcf.append(couplef)
            modelf = modeld
            path_ac_file = _InputData.twissfile.replace(".dat","_ac.dat")
            if not os.path.isfile(path_ac_file):
                print >> sys.stderr, "Ac file does not exist:", path_ac_file,"\nIn GUI check 'Ac dipole' box to create a model with ac dipole."
                sys.exit(1)
            modeld = metaclass.twiss(path_ac_file)
            if float(dpp)==0.0:
                zerobxf = betalistxf[dpp]
                zerobyf = betalistyf[dpp]

        #except:
        #       print "No free data"

    #
    # driven beta
    #

    print "Driven beta"

    #H
    fileobj = _chromFileWriter('beta', _join_with_output_path("chrombetax.out"), 'H')

    bpms = Utilities.bpm.intersect(listx)
    bpms = Utilities.bpm.model_intersect(bpms, modeld)
    _do_lin_reg_bet(fileobj, files_dict.keys(), betalistx, bpms, "H", zerobx, modeld)
    del fileobj

    #V
    fileobj = _chromFileWriter('beta', _join_with_output_path("chrombetay.out"), 'V')

    bpms = Utilities.bpm.intersect(listy)
    bpms = Utilities.bpm.model_intersect(bpms, modeld)
    _do_lin_reg_bet(fileobj, files_dict.keys(), betalisty, bpms, "V", zeroby, modeld)
    del fileobj

    print "Driven beta finished"

    #
    # driven coupling
    #
    print "Driven coupling"

    fileobj = _chromFileWriter('coupling', _join_with_output_path("chromcoupling.out"), '')

    bpms = Utilities.bpm.intersect(listc)
    bpms = Utilities.bpm.model_intersect(bpms, modeld)

    _do_linreg_coupling(couplelist, bpms, files_dict.keys(), fileobj)
    del fileobj

    print "Driven coupling finished"

    if freeswitch == 1:
        #
        # free beta
        #
        print "Free beta"
        #H
        fileobj = _chromFileWriter('beta', _join_with_output_path("chrombetax_free.out"), 'H')

        bpms = Utilities.bpm.intersect(listxf)
        bpms = Utilities.bpm.model_intersect(bpms, modelf)
        _do_lin_reg_bet(fileobj, files_dict.keys(), betalistxf, bpms, "H", zerobxf, modelf)

        #V
        fileobj = _chromFileWriter('beta', _join_with_output_path("chrombetay_free.out"), 'V')

        bpms = Utilities.bpm.intersect(listyf)
        bpms = Utilities.bpm.model_intersect(bpms, modelf)
        _do_lin_reg_bet(fileobj, files_dict.keys(), betalistyf, bpms, "V", zerobyf, modelf)

        print "Free beta finished"

        #
        # free coupling
        #
        print "Free coupling"

        fileobj = _chromFileWriter('coupling', _join_with_output_path("chromcoupling_free.out"), '')

        bpms = Utilities.bpm.intersect(listcf)
        bpms = Utilities.bpm.model_intersect(bpms, modelf)

        _do_linreg_coupling(couplelistf, bpms, files_dict.keys(), fileobj)

        print "Free coupling finished"


#===================================================================================================
# helper-functions
#===================================================================================================
def _madcreator(dpps, files_dict):
    '''
    :param dpps: list of dp/p to create model for
    :param dict files_dict: dpp_value --> corresponding_filenames
    '''

    # creating the DPP
    dppstring = ''
    dppstring_ac = ''
    for dpp in dpps:
        path_to_twiss_dpp = _join_with_output_path("twiss_"+str(dpp)+".dat")
        if not os.path.exists(path_to_twiss_dpp):
            dppstring = dppstring+'twiss, chrom,sequence='+_InputData.accel+', deltap='+str(dpp)+', file="'+path_to_twiss_dpp+'";\n'
            # if the model has twiss_ac.dat:
            if os.path.exists(_InputData.twissfile.replace(".dat","_ac.dat")): # this is only correct as long as the filenames are <filename>_ac.dat and <filename>.dat!
                path_to_twiss_ac_dpp = _join_with_output_path("twiss_"+str(dpp)+"_ac.dat")
                dppstring_ac = dppstring_ac+'twiss, chrom,sequence='+_InputData.accel+', deltap='+str(dpp)+', file="'+path_to_twiss_ac_dpp+'";\n'
            else: # do not create ac file if we don't have ac in our original model..
                dppstring_ac = ''

    if not dppstring:
        print "No need to run madx"
        return 0
    
    dict_for_replacing = {}
    dict_for_replacing["DPP"] = dppstring
    dict_for_replacing["DP_AC_P"] = dppstring_ac
    dict_for_replacing["ACCEL"] = _InputData.accel
    if _InputData.accel == 'LHCB1':
        dict_for_replacing["BEAM"] = "B1"
    elif _InputData.accel == 'LHCB2':
        dict_for_replacing["BEAM"] = "B2"
    else:
        print "WARNING: Could not decide what BEAM should be"

    (qx, qy, qdx, qdy) = _get_tunes(files_dict)

    dict_for_replacing["QX"] = qx
    dict_for_replacing["QY"] = qy
    dict_for_replacing["QDX"] = qdx
    dict_for_replacing["QDY"] = qdy
    dict_for_replacing["QMX"] = int(qx*1000000)
    dict_for_replacing["QMY"] = int(qy*1000000)
    dict_for_replacing["STOP"] = "!"

    for testpath in [_InputData.output_path, os.path.dirname(_InputData.twissfile)]:
        _tmpmod = os.path.join(testpath, 'modifiers.madx')
        if os.path.isfile(_tmpmod):
            print "INFO: Using", _tmpmod
            dict_for_replacing["MODIFIERS"] = _tmpmod
            break

    print "Creating madx"
    path_to_job_chrom_madx = _join_with_output_path("job.chrom.madx")
    Utilities.iotools.replace_keywords_in_textfile(
                                                   os.path.join(_InputData.path_to_beta_beat, "MODEL", "LHCB", "model", "job.twiss_chrom.madx.macro"), 
                                                   dict_for_replacing, 
                                                   path_to_job_chrom_madx
                                                   )

    print "Running madx"
    process = subprocess.Popen(_InputData.path_to_madx+' < '+path_to_job_chrom_madx,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)
    # wait for the process to terminate
    (out, err) = process.communicate()
    errcode = process.returncode

    if 0 != errcode:
        print "Mad-X failed. Printing output:-------------------------"
        print out
        print >> sys.stderr, "Mad-X failed. Printing error output:-------------------"
        print >> sys.stderr, err
        raise ValueError("Mad-X failed")

def _get_filenames_in_output_with_dpp(dpp=''):
    '''
    Returns list of available file names
    for the given dpp.

    Example: _get_filenames_in_output_with_dpp('0.0')

    :param dpp: [str] dpp appendix to list of files
    '''
    ret = []
    for fname in os.listdir(_InputData.output_path):
        if dpp:
            extra = '_'+dpp
        else:
            extra = ''
        # thank you for making this easy..
        if re.search("get[A-Za-z]*[0-9]*[xy]*(_free[2]*)*"+extra+".out", fname):
            ret.append(fname)
    return ret


def _rungetllm(twiss_filename, files, dpp):
    '''
    Running GetLLM...
    '''
    import GetLLM

    print "Will run getllm for ", dpp

    GetLLM.main(outputpath=_InputData.output_path,
            files_to_analyse=','.join(files),
            model_filename=twiss_filename,
            accel=_InputData.accel,
            TBTana=_InputData.technique)
    print "GetLLM finished"

    for fname in _get_filenames_in_output_with_dpp():
        src_path = _join_with_output_path(fname)
        dst_path = _join_with_output_path(fname.replace(".out", "_"+str(dpp)+".out"))
        shutil.move(src_path, dst_path)


def _copy_default_outfiles():
    for fname in _get_filenames_in_output_with_dpp('0.0'):
        src_path = _join_with_output_path(fname)
        dst_path = _join_with_output_path(fname.replace("_0.0.out", ".out"))
        shutil.copy(src_path, dst_path)


##### for chromatic

def _do_lin_reg_bet(fileobj, listx, listy, bpms, plane, zero, twiss):
    '''
    Calculates stuff and writes to the file in a table
    Closes the file afterwards

    :param fileobj: _chromFileWriter for output table
    :param listx: List of variables...
    :param listy: List of variables...
    :param bpms: List of BPMs
    :param plane: Which plane (H/V)
    :param zero: Twiss for dp/p = 0
    :param twiss: Twiss
    '''
    for bpm in bpms:
        name = bpm[1]
        sloc = bpm[0]
        indx = []
        b = []
        a = []
        bm = []
        am = []
        if "H" in plane:
            beta0 = zero.BETX[zero.indx[name]]
            alfa0 = zero.ALFX[zero.indx[name]]
            alfa0err = zero.STDALFX[zero.indx[name]]

            beta0m = twiss.BETX[twiss.indx[name]]
            alfa0m = twiss.ALFX[twiss.indx[name]]

            wmo = twiss.WX[twiss.indx[name]]
            pmo = twiss.PHIX[twiss.indx[name]]
        else:

            beta0 = zero.BETY[zero.indx[name]]
            alfa0 = zero.ALFY[zero.indx[name]]
            alfa0err = zero.STDALFY[zero.indx[name]]

            beta0m = twiss.BETY[twiss.indx[name]]
            alfa0m = twiss.ALFY[twiss.indx[name]]

            wmo = twiss.WY[twiss.indx[name]]
            pmo = twiss.PHIY[twiss.indx[name]]
        for dpp in listx:
            _file = listy[dpp]
            ix = _file.indx[name]
            indx.append(ix)
            if "H" in plane:
                b.append(_file.BETX[ix])
                a.append(_file.ALFX[ix])

                bm.append(_file.BETXMDL[_file.indx[name]])
                am.append(_file.ALFXMDL[_file.indx[name]])
            else:
                b.append(_file.BETY[ix])
                a.append(_file.ALFY[ix])

                bm.append(_file.BETYMDL[_file.indx[name]])
                am.append(_file.ALFYMDL[_file.indx[name]])

        bfit = linreg.linreg(listx, b)
        afit = linreg.linreg(listx, a)

        bfitm = linreg.linreg(listx, bm)
        afitm = linreg.linreg(listx, am)

        # measurement
        dbb = bfit[0]/beta0
        dbberr = bfit[3]/beta0
        da = afit[0]
        daerr = afit[3]
        A = dbb
        Aerr = dbberr
        B = da-alfa0*dbb
        Berr = math.sqrt(daerr**2 + (alfa0err*dbb)**2 + (alfa0*dbberr)**2)
        w = 0.5*math.sqrt(A**2+B**2)
        werr = 0.5*math.sqrt( (Aerr*A/w)**2 + (Berr*B/w)**2  )
        phi = math.atan2(B,A)/2./math.pi
        phierr = 1./(1.+(A/B)**2)*math.sqrt( (Aerr/B)**2 + (A/B**2*Berr)**2)/2./math.pi

        #model
        dbbm = bfitm[0]/beta0m
        dbberrm = bfitm[3]/beta0m
        dam = afitm[0]
        daerrm = afitm[3]
        Am = dbbm
        Aerrm = dbberrm
        Bm = dam-alfa0m*dbbm
        Berrm = math.sqrt(daerrm**2 + (alfa0m*dbberrm)**2)
        wm = 0.5*math.sqrt(Am**2+Bm**2)
        werrm = 0.5*math.sqrt( (Aerrm*Am/wm)**2 + (Berrm*Bm/wm)**2  )
        phim = math.atan2(Bm,Am)/2./math.pi
        phierrm = 1./(1.+(Am/Bm)**2)*math.sqrt( (Aerrm/Bm)**2 + (Am/Bm**2*Berrm)**2)/2./math.pi

        fileobj.writeLine(locals().copy())


def _get_f( couplelist, dpplist, bpm_name, value):
    '''
    calculates the linear regression of 'value' for each
    dpp in dpplist

    :param couplelist: list of getcouple files (for each dpp)
    :param dpplist: list of all dpp values available
    :param bpm_name: name of bpm
    :param value: name of column (e.g. F1001R)
    '''
    lst = []
    x = []
    for dpp in dpplist:
        x.append(dpp)
        couplefile = couplelist[dpp]
        lst.append( getattr(couplefile, value)[ couplefile.indx[ bpm_name]])
    lreg = linreg.linreg(x, lst)
    return lreg[0], lreg[3]


def _do_linreg_coupling(couplelist, bpms, dpplist, fileobj):
    '''
    linreg for chromatic coupling

    Writes to fileobj the chromatic coupling.
    f1001, f1010 derivatives wrt dp/p, and errors.
    '''
    for bpm in bpms:
        name = bpm[1]
        sloc = bpm[0]

        chr_f1001r, chr_err_f1001r = _get_f( couplelist, dpplist, name, 'F1001R')
        chr_f1001i, chr_err_f1001i = _get_f( couplelist, dpplist, name, 'F1001I')
        chr_f1010r, chr_err_f1010r = _get_f( couplelist, dpplist, name, 'F1010R')
        chr_f1010i, chr_err_f1010i = _get_f( couplelist, dpplist, name, 'F1010I')

        mdl_chr_f1001r, mdl_chr_err_f1001r = _get_f( couplelist, dpplist, name, 'MDLF1001R')
        mdl_chr_f1001i, mdl_chr_err_f1001i = _get_f( couplelist, dpplist, name, 'MDLF1001I')
        mdl_chr_f1010r, mdl_chr_err_f1010r = _get_f( couplelist, dpplist, name, 'MDLF1010R')
        mdl_chr_f1010i, mdl_chr_err_f1010i = _get_f( couplelist, dpplist, name, 'MDLF1010I')

        fileobj.writeLine(locals().copy())


def _get_tunes(fileslist):
    '''
    Reads in the driven tunes from the
    file with dpp=0
    Reads in the model tunes from the
    twiss model (twiss.dat)

    :param fileslist: dictionary of files, dpp used as key
    :returns: (
    :raise ValueError: If fileslist[0] does not exist
    '''
    if fileslist[0][0].endswith(".gz"):
        fname = fileslist[0][0].replace(".gz", "")
        end = '.gz'
    else:
        fname = fileslist[0][0]
        end = ''
    tw_x = metaclass.twiss(fname+'_linx'+end)
    tw_y = metaclass.twiss(fname+'_liny'+end)
    tw = metaclass.twiss(_InputData.twissfile)

    qdx = tw_x.TUNEX[0]
    qdy = tw_y.TUNEY[0]
    qx = tw.Q1 % 1
    qy = tw.Q2 % 1

    return (qx, qy, qdx, qdy)


def _join_with_output_path(*path_tokens):
    return os.path.join(_InputData.output_path, *path_tokens)


class _chromFileWriter:
    def __init__(self, ftype, fname, plane, overwrite=True):
        '''
        :param ftype: string, 'beta' or 'coupling'
        :param fname: string, name of file
        :param plane: "H" or "V"
        :param overwrite: Overwrite file if it already exist
        '''
        self.fstream = file(fname,'w')
        self._clen = 19 # column length..

        if plane.upper()=="H":
            plane = "X"
        elif plane.upper()=="V":
            plane = "Y"

        if os.path.isfile(fname) and not overwrite:
            raise ValueError("Cannot overwrite file "+fname)

        betacolumns = ['name', 'sloc',  'dbb', 'dbberr', 'da', 'daerr',
                    'w', 'werr','wm', 'werrm', 'wmo', 'phi', 'phierr','phim','phierrm','pmo',
                    'A','Aerr','Am','Aerrm','B','Berr','Bm','Berrm', 'dbbm','dbberrm',
                    'dam','daerrm']
        couplecolumns = ['name', 'sloc',
            'chr_f1001r', 'chr_err_f1001r', 'chr_f1001i', 'chr_err_f1001i',
            'chr_f1010r', 'chr_err_f1010r', 'chr_f1010i', 'chr_err_f1010i',
            'mdl_chr_f1001r', 'mdl_chr_err_f1001r', 'mdl_chr_f1001i', 'mdl_chr_err_f1001i',
            'mdl_chr_f1010r', 'mdl_chr_err_f1010r', 'mdl_chr_f1010i', 'mdl_chr_err_f1010i',
          ]

        headnames = {
            'name':'NAME',
            'sloc':'S',
            'dbb':'dbb',
            'dbberr':'dbberr',
            'da':'dalfa',
            'daerr':'daerr',
            'w':'W%(plane)s',
            'werr':'W%(plane)sERR',
            'wmo':'WMO',
            'phi':'PHI%(plane)s',
            'phierr':'PHI%(plane)sERR',
            'phim':'PHI%(plane)sM',
            'phierrm':'PHI%(plane)sERRM',
            'pmo':'PHIZERO',
            'dbbm':'dbbM',
            'dbberrm':'dbberrM',
            'dam':'dalfaM',
            'daerrm':'daerrM',
            'wm':'W%(plane)sM',
            'werrm':'W%(plane)sERRM',
            'A':'CHROM_A%(plane)s',
            'Aerr':'CHROM_Aerr%(plane)s',
            'Am':'CHROM_AM%(plane)s',
            'Aerrm':'CHROM_AERRM%(plane)s',
            'B':'CHROM_B%(plane)s',
            'Berr':'CHROM_Berr%(plane)s',
            'Bm':'CHROM_BM%(plane)s',
            'Berrm':'CHROM_BERRM%(plane)s',
            'chr_f1001r':'Cf1001r',
            'chr_err_f1001r':'Cf1001rERR',
            'chr_f1001i':'Cf1001i',
            'chr_err_f1001i':'Cf1001iERR',
            'chr_f1010r':'Cf1010r',
            'chr_err_f1010r':'Cf1010rERR',
            'chr_f1010i':'Cf1010i',
            'chr_err_f1010i':'Cf1010iERR',
            'mdl_chr_f1001r':'Cf1001r_MDL',
            'mdl_chr_err_f1001r':'Cf1001rERR_MDL',
            'mdl_chr_f1001i':'Cf1001i_MDL',
            'mdl_chr_err_f1001i':'Cf1001iERR_MDL',
            'mdl_chr_f1010r':'Cf1010r_MDL',
            'mdl_chr_err_f1010r':'Cf1010rERR_MDL',
            'mdl_chr_f1010i':'Cf1010i_MDL',
            'mdl_chr_err_f1010i':'Cf1010iERR_MDL'}
        headtypes = {'name':'%s'}
        for key in headnames:
            # for all others we use '%le'...
            if key not in headtypes:
                headtypes[key] = '%le'

        if ftype.lower().strip()=='beta':
            self.ftype = 'beta'
            self.columns = betacolumns
        elif ftype.lower().strip()=='coupling':
            self.ftype = 'coupling'
            self.columns = couplecolumns
        else:
            raise ValueError("ftype %s not understood" % (ftype))

        self.head = [headnames[c] % locals() for c in self.columns]
        headcount = [i+1 for i in xrange(len(self.head))]
        self.types = [headtypes[c] for c in self.columns]

        self.head[0] = '* '+(self.head[0].rjust(self._clen-3))
        headcount[0] = '# '+(str(headcount[0]).rjust(self._clen-3))
        self.types[0] = '$ '+(self.types[0].rjust(self._clen-3))

        self._write_list(self.head)
        # According to SL-CO-Note-91-32, this is allowed...
        self._write_list(headcount)
        self._write_list(self.types)

    def writeLine(self, data):
        '''
        Write one line.
        :param data: Dictionary of data columns to write
        '''
        tmp_line = [data[c] for c in self.columns]
        self._write_list(tmp_line)

    def _write_list(self, tmp_list):
        line = ''
        for entry in tmp_list:
            line += str(entry).rjust(self._clen-1)+' '
        self.fstream.write(line[:-1]+'\n')


class _InputData(object):
    """This static class holds all input variables for getsuper """
    source_files = ""
    path_to_madx = ""
    twissfile = ""
    output_path = ""
    path_to_beta_beat = ""
    technique = ""
    accel = ""
    deltap_scaling_factor = 0.0

    @staticmethod
    def static_init(source_files, path_to_madx, twissfile, output_path, path_to_beta_beat, technique, accel, deltap_scaling_factor):
        _InputData.source_files = source_files
        _InputData.__check_src_files()
        _InputData.path_to_madx = path_to_madx
        _InputData.twissfile = superutils.get_twissfile(twissfile)
        if not os.path.isfile(_InputData.twissfile):
            raise ValueError("Twissfile does not exist: "+_InputData.twissfile)
        _InputData.output_path = output_path
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        _InputData.path_to_beta_beat = path_to_beta_beat
        _InputData.technique = technique
        if technique not in ("SUSSIX", "SVD"):
            raise ValueError("Chosen  technique is unknown: "+technique)
        _InputData.accel = accel
        if _InputData.accel not in ("LHCB1", "LHCB2", "SPS", "RHIC"):
            raise ValueError("Provided accelerator is not valid: "+_InputData.accel)
        _InputData.deltap_scaling_factor = deltap_scaling_factor


    @staticmethod
    def __check_src_files():
        if 2 > len(_InputData.source_files):
            raise ValueError("Provide at least two files. Files: "+str(_InputData.source_files))
        for f_name in _InputData.source_files:
            if not os.path.isfile(f_name) and not os.path.isfile(f_name+'.gz'):
                raise ValueError(f_name+' does not exist')

    def __init__(self):
        raise NotImplementedError("static class _InputData cannot be instantiated")


#===================================================================================================
# main invocation
#===================================================================================================
def _start():
    """ Starter function to not pollute the global space with variables options and args """
    options,args = _parse_args()
    main(
        source_files=superutils.get_filelist(options, args),
        path_to_madx=options.madx,
        twissfile=options.twissfile,
        output_path=options.output,
        path_to_beta_beat=options.brc,
        technique=options.technique,
        accel=options.accel,
        deltap_scaling_factor=options.deltap_scaling_factor
        )

if __name__ == "__main__":
    _start()
