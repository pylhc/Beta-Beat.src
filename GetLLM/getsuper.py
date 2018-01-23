"""
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

"""

import argparse
import os
import sys
import shutil
import math
import re
import numpy as np

import __init__  # @UnusedImport init will include paths
import Python_Classes4MAD.metaclass as metaclass
from Utilities import bpm as bpm_util
from Utilities import logging_tools as logtools
from Utilities.dict_tools import DotDict
from Utilities import tfs_remove_nan
from model import manager, creator
from model.accelerators.lhc import LhcExcitationMode

LOG = logtools.get_logger(__name__)


# ==================================================================================================
# Argument Handling
# ==================================================================================================

ALGO_CHOICES = ["SUSSIX", "SVD", "HA"]


def _parse_args(args=None):
    """ Parses arguments from command line. """
    parser = argparse.ArgumentParser()

    # general
    parser.add_argument("-f", "--files",
            help="Files from analysis, separated by comma",
            metavar="TwissFile", dest="files", required=True)
    parser.add_argument("--twissfile",
            help="Twiss file to use",
            metavar="/path/to/twiss.dat", dest="twissfile", required=True)
    parser.add_argument("-o", "--output",
            help="Output path, where to store the results",
            metavar="<path>", default="./", dest="output_path")
    parser.add_argument("-t", "--algorithm",
            help="Which algorithm to use {:s}".format(ALGO_CHOICES),
            metavar="ALGORITHM", default=ALGO_CHOICES[0], dest="algorithm",
                        choices=ALGO_CHOICES)
    parser.add_argument("-d", "--deltapScalingFactor",
            help="Scaling factor for deltap, remember final value must be in MAD units",
            metavar="<deltapScalingFactor>", default=1.0, type=float, dest="deltap_scaling_factor")

    # parse arguments
    accel_cls, remain_args = manager.get_accel_class_from_args(args)
    options = parser.parse_args(remain_args)
    source_files = [f.strip() for f in options.files.split(',')]

    # put all arguments into one dict
    options_dict = {
        "accel_cls": accel_cls,
        "source_files": source_files,
    }
    options_dict.update(options.__dict__)

    options_dict.pop("files")  # is "source_files" now

    return options_dict


def check_input(opt):
    """ Approves the input and sets default """

    # files
    if "source_files" not in opt or len(opt.source_files) < 2:
        raise ValueError("Provide at least two source files!")
    for f_name in opt.source_files:
        if not os.path.isfile(f_name) and not os.path.isfile(f_name + '.gz'):
            raise ValueError(f_name + ' does not exist')

    if "twissfile" not in opt:
        raise ValueError("Twissfile needed for execution.")

    if not os.path.isfile(opt.twissfile):
        raise ValueError("Twissfile does not exist: " + opt.twissfile)

    # set defaults
    opt.deltap_scaling_factor = opt.get("deltap_scaling_factor", 1.0)

    opt.output_path = opt.get("output_path", "./")
    if not os.path.isdir(opt.output_path):
        os.makedirs(opt.output_path)

    opt.algorithm = opt.get("algorithm", ALGO_CHOICES[0])
    if opt.algorithm not in ALGO_CHOICES:
        raise ValueError("Algorithm needs to be either one of  '" + ALGO_CHOICES + "'")

    return opt


# ==================================================================================================
# main()-function
# ==================================================================================================

def main(**kwargs):
    """ getsuper main function

    Keyword Args:
        source_files (list): list of strings with file_paths to TFS files
        twissfile (str): path to TFS model file
        output_path (str): Path to store created files
        algorithm (str): Used Turn-by-turn data analysis algorithm: 'SUSSIX', 'SVD' or 'HA'
        accel_cls (accelerator): Accelerator class object
        deltap_scaling_factor (float): Scaling factor for deltap,
                                       remember final value must be in MAD units
    """

    options = check_input(DotDict(kwargs))

    files_dict = {}  # dpp --> files with corresponding dpp

    for f_name in options.source_files:

        datax, datay = _load_from_file(f_name)

        dppx = datax.DPP * options.deltap_scaling_factor  # scaling factor hack for old files
        dppy = datay.DPP * options.deltap_scaling_factor

        if dppx != dppy:
            raise ValueError("Discrepancy between horizontal"
                             "{:f} and vertical {:f} dpp".format(dppx, dppy))
        else:
            dpp = float(dppx)

        if dpp not in files_dict:
            LOG.debug("Adding dpp {:f}".format(dpp))
            files_dict[dpp] = [f_name]
        else:
            files_dict[dpp].append(f_name)

    if len(files_dict.keys()) < 2:
        raise ValueError("Less than two DPP-Values found. Cannot do W-Analysis.")

    if 0 not in files_dict:
        raise ValueError("NO DPP=0.0. Provide at least one source file with DPP=0.0.")

    accel_inst = _create_accel_instance(options.accel_cls, files_dict.keys(),
                                        files_dict, options.output_path, options.twissfile)
    LOG.debug("All models are created")

    for dpp in files_dict:
        files = files_dict[dpp]
        twiss_dpp_path = _join_with_output_path(options.output_path, "twiss_{:f}.dat".format(dpp))
        _rungetllm(twiss_dpp_path, files, dpp,
                   options.output_path, accel_inst, options.algorithm)

    # The GUI wants the default files to have the names without _0.0
    _copy_default_outfiles(options.output_path)

    #TODO: HOPE THAT GETLLM DOES A BETTER JOB
    LOG.warn("Cleaning files of NAN! This should not be necessary!")
    all_files = os.listdir(options.output_path)
    tfs_remove_nan.clean_files(
        [os.path.join(options.output_path, f) for f in all_files], replace=True)

    # adding data
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

    for dpp in files_dict.keys():
        LOG.debug("Loading driven data for dpp {:f}".format(dpp))
        betx_path = _join_with_output_path(options.output_path, 'getbetax{ext:s}'.format(ext=_ext(dpp)))
        bety_path = _join_with_output_path(options.output_path, 'getbetay{ext:s}'.format(ext=_ext(dpp)))
        couple_path = _join_with_output_path(options.output_path, 'getcouple{ext:s}'.format(ext=_ext(dpp)))

        betx = metaclass.twiss(betx_path)
        bety = metaclass.twiss(bety_path)
        couple = metaclass.twiss(couple_path)

        betalistx[dpp] = betx
        betalisty[dpp] = bety
        couplelist[dpp] = couple

        if float(dpp) == 0.0:
            zerobx = betx
            zeroby = bety

        listx.append(betx)
        listy.append(bety)
        listc.append(couple)

        modeld = metaclass.twiss(options.twissfile)

        try:
            betaxf_path = _join_with_output_path(options.output_path, 'getbetax_free{ext:s}'.format(ext=_ext(dpp)))
            betxf = metaclass.twiss(betaxf_path)
            LOG.debug("Loaded betax free data from '{:s}".format(betaxf_path))


            betayf_path = _join_with_output_path(options.output_path, 'getbetay_free{ext:s}'.format(ext=_ext(dpp)))
            betyf = metaclass.twiss(betayf_path)
            LOG.debug("Loaded betay free data from '{:s}".format(betayf_path))

            couplef_path = _join_with_output_path(options.output_path, 'getcouple_free{ext:s}'.format(ext=_ext(dpp)))
            couplef = metaclass.twiss(couplef_path)
            LOG.debug("Loaded coupling free data from '{:s}".format(couplef_path))
        except IOError:
            use_free = False
            LOG.warn("WARNING: Could not open all of the free data files.")
        else:
            use_free = True
            betalistxf[dpp] = betxf
            betalistyf[dpp] = betyf
            couplelistf[dpp] = couplef

            listxf.append(betxf)
            listyf.append(betyf)
            listcf.append(couplef)

            modelf = modeld

            path_ac_file = options.twissfile.replace(".dat", "_ac.dat")
            if not os.path.isfile(path_ac_file):
                LOG.error("Ac file '{:s}' does not exist.".format(path_ac_file))
                LOG.error("  -> In GUI check 'Ac dipole' box to create a model with ac dipole.")
                sys.exit(1)
            modeld = metaclass.twiss(path_ac_file)
            if float(dpp) == 0.0:
                zerobxf = betalistxf[dpp]
                zerobyf = betalistyf[dpp]

    LOG.debug("Getting Driven beta")
    # H
    fileobj = _chromFileWriter('beta', _join_with_output_path(options.output_path, "chrombetax" + _ext()), 'H')
    bpms = bpm_util.intersect(listx)
    bpms = bpm_util.model_intersect(bpms, modeld)
    _do_lin_reg_bet(fileobj, files_dict.keys(), betalistx, bpms, "H", zerobx, modeld)
    del fileobj

    # V
    fileobj = _chromFileWriter('beta', _join_with_output_path(options.output_path, "chrombetay" + _ext()), 'V')
    bpms = bpm_util.intersect(listy)
    bpms = bpm_util.model_intersect(bpms, modeld)
    _do_lin_reg_bet(fileobj, files_dict.keys(), betalisty, bpms, "V", zeroby, modeld)
    del fileobj
    LOG.debug("Driven beta finished")


    LOG.debug("Getting Driven coupling")
    fileobj = _chromFileWriter('coupling', _join_with_output_path(options.output_path, "chromcoupling" + _ext()), '')
    bpms = bpm_util.intersect(listc)
    bpms = bpm_util.model_intersect(bpms, modeld)
    _do_linreg_coupling(couplelist, bpms, files_dict.keys(), fileobj)
    del fileobj
    LOG.debug("Driven coupling finished")

    if use_free:
        #
        # free beta
        #
        LOG.debug("Getting Free beta")
        # H
        fileobj = _chromFileWriter('beta',
                                   _join_with_output_path(options.output_path, "chrombetax_free" + _ext()), 'H')
        bpms = bpm_util.intersect(listxf)
        bpms = bpm_util.model_intersect(bpms, modelf)
        _do_lin_reg_bet(fileobj, files_dict.keys(), betalistxf, bpms, "H", zerobxf, modelf)

        # V
        fileobj = _chromFileWriter('beta',
                                   _join_with_output_path(options.output_path, "chrombetay_free" + _ext()), 'V')
        bpms = bpm_util.intersect(listyf)
        bpms = bpm_util.model_intersect(bpms, modelf)
        _do_lin_reg_bet(fileobj, files_dict.keys(), betalistyf, bpms, "V", zerobyf, modelf)
        LOG.debug("Free beta finished")


        LOG.debug("GettingFree coupling")
        fileobj = _chromFileWriter('coupling',
                                   _join_with_output_path(options.output_path, "chromcoupling_free" + _ext()), '')
        bpms = bpm_util.intersect(listcf)
        bpms = bpm_util.model_intersect(bpms, modelf)
        _do_linreg_coupling(couplelistf, bpms, files_dict.keys(), fileobj)
        LOG.debug("Free coupling finished")

# ===================================================================================================
# helper-functions
# ===================================================================================================


def _load_from_file(filepath):
    if filepath.endswith('.gz'):
        try:
            datax = metaclass.twiss(filepath.replace(".gz", ".linx.gz"))
            datay = metaclass.twiss(filepath.replace(".gz", ".liny.gz"))
        except IOError:
            datax = metaclass.twiss(filepath.replace(".gz", "_linx.gz"))
            datay = metaclass.twiss(filepath.replace(".gz", "_liny.gz"))
    else:
        try:
            datax = metaclass.twiss(filepath + ".linx")
            datay = metaclass.twiss(filepath + ".liny")
        except IOError:
            datax = metaclass.twiss(filepath + "_linx")
            datay = metaclass.twiss(filepath + "_liny")
    return datax, datay


def _create_accel_instance(accel_cls, dpps, files_dict, output_path, twissfile):
    """
    Args:
        dpps: list of dp/p to create model for
        dict files_dict: dpp_value --> corresponding_filenames
    """

    (exp_qx, exp_qy, mdl_qx, mdl_qy) = _get_tunes(twissfile, files_dict)

    for testpath in [output_path, os.path.dirname(twissfile)]:
        modifiers = os.path.join(testpath, 'modifiers.madx')
        if os.path.isfile(modifiers):
            LOG.debug("Using file '{:s}'".format(modifiers))
            break

    accel_inst = accel_cls()
    accel_inst.optics_file = modifiers
    accel_inst.nat_tune_x = mdl_qx
    accel_inst.nat_tune_y = mdl_qy
    accel_inst.drv_tune_x = exp_qx
    accel_inst.drv_tune_y = exp_qy
    accel_inst.dpp = dpps
    accel_inst.xing = False


    accel_inst.excitation = LhcExcitationMode.FREE
    if os.path.exists(twissfile.replace(".dat", "_ac.dat")):
        accel_inst.excitation = LhcExcitationMode.ACD
    elif os.path.exists(twissfile.replace(".dat", "_adt.dat")):
        accel_inst.excitation = LhcExcitationMode.ADT

    creator.create_model(accel_inst, "nominal", output_path)

    return accel_inst


def _get_output_filenames(output_path, dpp=None):
    """
    Returns list of available file names
    for the given dpp.

    Example: _get_output_filenames(0.)

    Args:
        dpp: dpp appendix to list of files
    """
    ret = []
    for fname in os.listdir(output_path):
        ext = _ext(dpp)
        if re.match(r"get[^_]+[_free\d?]?" + ext, fname):
            ret.append(fname)
    return ret


def _rungetllm(twiss_filename, files, dpp, output_path, accel_inst, algorithm):
    """
    Running GetLLM...
    """
    import GetLLM

    LOG.debug("Will run getllm for dpp {:f}".format(dpp))

    if "lhc" == accel_inst.NAME:
        lhcphase = "1"
        accel_name = "LHCB" + str(accel_inst.get_beam())
    else:
        lhcphase = "0"
        accel_name = accel_inst.NAME.upper()  #TODO: TEST that! Should work for ESRF at least.

    GetLLM.main(outputpath=output_path,
            files_to_analyse=','.join(files),
            model_filename=twiss_filename,
            accel=accel_name,
            tbtana=algorithm,
            lhcphase=lhcphase)
    LOG.debug("GetLLM finished")

    for fname in _get_output_filenames(output_path):
        src_path = _join_with_output_path(output_path, fname)
        dst_path = _join_with_output_path(output_path, fname.replace(_ext(), _ext(dpp)))
        shutil.move(src_path, dst_path)


def _copy_default_outfiles(output_path):
    for fname in _get_output_filenames(output_path, dpp=0.):
        src_path = _join_with_output_path(output_path, fname)
        dst_path = _join_with_output_path(output_path, fname.replace(_ext(0.), _ext()))
        shutil.copy(src_path, dst_path)


def _ext(dpp=None):
    if dpp is None:
        return '.out'
    else:
        return "_{:f}.out".format(dpp)


# for chromatic

def _do_lin_reg_bet(fileobj, listx, listy, bpms, plane, zero, twiss):
    """
    Calculates stuff and writes to the file in a table
    Closes the file afterwards

    Args:
        fileobj: _chromFileWriter for output table
        listx: List of variables...
        listy: List of variables...
        bpms: List of BPMs
        plane: Which plane (H/V)
        zero: Twiss for dp/p = 0
        twiss: Twiss
    """
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
            try:
                alfa0err = zero.STDALFX[zero.indx[name]]
            except AttributeError:
                alfa0err = zero.ERRALFX[zero.indx[name]]


            beta0m = twiss.BETX[twiss.indx[name]]
            alfa0m = twiss.ALFX[twiss.indx[name]]

            wmo = twiss.WX[twiss.indx[name]]
            pmo = twiss.PHIX[twiss.indx[name]]
        else:

            beta0 = zero.BETY[zero.indx[name]]
            alfa0 = zero.ALFY[zero.indx[name]]
            try:
                alfa0err = zero.STDALFY[zero.indx[name]]
            except AttributeError:
                alfa0err = zero.ERRALFY[zero.indx[name]]

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

        bfit = linreg(listx, b)
        afit = linreg(listx, a)

        bfitm = linreg(listx, bm)
        afitm = linreg(listx, am)

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


def _get_f(couplelist, dpplist, bpm_name, value):
    """
    calculates the linear regression of 'value' for each
    dpp in dpplist

    Args:
        couplelist: list of getcouple files (for each dpp)
        dpplist: list of all dpp values available
        bpm_name: name of bpm
        value: name of column (e.g. F1001R)
    """
    lst = []
    x = []
    for dpp in dpplist:
        x.append(dpp)
        couplefile = couplelist[dpp]
        lst.append(getattr(couplefile, value)[couplefile.indx[bpm_name]])

    lreg = linreg(x, lst)

    return lreg[0], lreg[3]


def _do_linreg_coupling(couplelist, bpms, dpplist, fileobj):
    """
    linreg for chromatic coupling

    Writes to fileobj the chromatic coupling.
    f1001, f1010 derivatives wrt dp/p, and errors.
    """
    for bpm in bpms:
        name = bpm[1]
        sloc = bpm[0]

        chr_f1001r, chr_err_f1001r = _get_f(couplelist, dpplist, name, 'F1001R')
        chr_f1001i, chr_err_f1001i = _get_f(couplelist, dpplist, name, 'F1001I')
        chr_f1010r, chr_err_f1010r = _get_f(couplelist, dpplist, name, 'F1010R')
        chr_f1010i, chr_err_f1010i = _get_f(couplelist, dpplist, name, 'F1010I')

        mdl_chr_f1001r, mdl_chr_err_f1001r = _get_f(couplelist, dpplist, name, 'MDLF1001R')
        mdl_chr_f1001i, mdl_chr_err_f1001i = _get_f(couplelist, dpplist, name, 'MDLF1001I')
        mdl_chr_f1010r, mdl_chr_err_f1010r = _get_f(couplelist, dpplist, name, 'MDLF1010R')
        mdl_chr_f1010i, mdl_chr_err_f1010i = _get_f(couplelist, dpplist, name, 'MDLF1010I')

        fileobj.writeLine(locals().copy())


def _get_tunes(model_file, fileslist):
    """
    Reads in the driven tunes from the
    file with dpp=0
    Reads in the model tunes from the
    twiss model (twiss.dat)

    Args:
        fileslist: dictionary of files, dpp used as key

    Returns:
        (qx, qy, qdx, qdy, qmx, qmy)

    Raises:
        ValueError: If fileslist[0] does not exist
    """
    tw_x, tw_y = _load_from_file(fileslist[0][0])
    tw = metaclass.twiss(model_file)

    exp_qx = tw_x.Q1
    exp_qy = tw_y.Q2
    mdl_qx = tw.Q1
    mdl_qy = tw.Q2

    return (exp_qx, exp_qy, mdl_qx, mdl_qy)


def linreg(X, Y):
    """
    Summary
        Linear regression of y = ax + b
    Usage
        real, real, real = linreg(list, list)
    Returns coefficients to the regression line "y=ax+b" from x[] and y[], and R^2 Value
    """
    if len(X) != len(Y):  raise ValueError, 'unequal length'
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in map(None, X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
    meanerror = residual = 0.0
    for x, y in map(None, X, Y):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2
    if residual == 0 and meanerror == 0:
        RR = 1.0
    else:
        RR = 1 - residual/meanerror
    if N>2:
        ss = residual / (N-2)
    else:
        ss = 0
    Var_a, Var_b = ss * N / det, ss * Sxx / det
    #print "y=ax+b"
    #print "N= %d" % N
    #print "a= %g \pm t_{%d;\alpha/2} %g" % (a, N-2, sqrt(Var_a))
    #print "b= %g \pm t_{%d;\alpha/2} %g" % (b, N-2, sqrt(Var_b))
    #print "R^2= %g" % RR
    #print "s^2= %g" % ss
    return a, b, RR, np.sqrt(Var_a), np.sqrt(Var_b)


def _join_with_output_path(output_path, *path_tokens):
    return os.path.join(output_path, *path_tokens)


class _chromFileWriter:
    def __init__(self, ftype, fname, plane, overwrite=True):
        """
        Args:
            ftype: string, 'beta' or 'coupling'
            fname: string, name of file
            plane: "H" or "V"
            overwrite: Overwrite file if it already exist
        """
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
            'name': 'NAME',
            'sloc': 'S',
            'dbb': 'dbb',
            'dbberr': 'dbberr',
            'da': 'dalfa',
            'daerr': 'daerr',
            'w': 'W%(plane)s',
            'werr': 'W%(plane)sERR',
            'wmo': 'WMO',
            'phi': 'PHI%(plane)s',
            'phierr': 'PHI%(plane)sERR',
            'phim': 'PHI%(plane)sM',
            'phierrm': 'PHI%(plane)sERRM',
            'pmo': 'PHIZERO',
            'dbbm': 'dbbM',
            'dbberrm': 'dbberrM',
            'dam': 'dalfaM',
            'daerrm': 'daerrM',
            'wm': 'W%(plane)sM',
            'werrm': 'W%(plane)sERRM',
            'A': 'CHROM_A%(plane)s',
            'Aerr': 'CHROM_Aerr%(plane)s',
            'Am': 'CHROM_AM%(plane)s',
            'Aerrm': 'CHROM_AERRM%(plane)s',
            'B': 'CHROM_B%(plane)s',
            'Berr': 'CHROM_Berr%(plane)s',
            'Bm': 'CHROM_BM%(plane)s',
            'Berrm': 'CHROM_BERRM%(plane)s',
            'chr_f1001r': 'Cf1001r',
            'chr_err_f1001r': 'Cf1001rERR',
            'chr_f1001i': 'Cf1001i',
            'chr_err_f1001i': 'Cf1001iERR',
            'chr_f1010r': 'Cf1010r',
            'chr_err_f1010r': 'Cf1010rERR',
            'chr_f1010i': 'Cf1010i',
            'chr_err_f1010i': 'Cf1010iERR',
            'mdl_chr_f1001r': 'Cf1001r_MDL',
            'mdl_chr_err_f1001r': 'Cf1001rERR_MDL',
            'mdl_chr_f1001i': 'Cf1001i_MDL',
            'mdl_chr_err_f1001i': 'Cf1001iERR_MDL',
            'mdl_chr_f1010r': 'Cf1010r_MDL',
            'mdl_chr_err_f1010r': 'Cf1010rERR_MDL',
            'mdl_chr_f1010i': 'Cf1010i_MDL',
            'mdl_chr_err_f1010i': 'Cf1010iERR_MDL'}
        headtypes = {'name': '%s'}
        for key in headnames:
            # for all others we use '%le'...
            if key not in headtypes:
                headtypes[key] = '%le'

        if ftype.lower().strip() == 'beta':
            self.ftype = 'beta'
            self.columns = betacolumns
        elif ftype.lower().strip() == 'coupling':
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
        """
        Write one line.

        Args:
            data: Dictionary of data columns to write
        """
        tmp_line = [data[c] for c in self.columns]
        self._write_list(tmp_line)

    def _write_list(self, tmp_list):
        line = ''
        for entry in tmp_list:
            line += str(entry).rjust(self._clen-1)+' '
        self.fstream.write(line[:-1]+'\n')


# ==================================================================================================
# main invocation
# ==================================================================================================

def _start(args=None):
    """ Starter function to not pollute the global space with variables options and args """
    options = _parse_args(args)
    main(**options)

if __name__ == "__main__":
    _start()
