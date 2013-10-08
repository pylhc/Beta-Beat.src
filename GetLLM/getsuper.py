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

Use argument --help for further information

Changelog:

 - 3.1 ylevinse 2012: 
    Cleaned macro writer in madcreator
    modifiers.madx should be in options.output
 - 3.1.1 vimaier 28th May 2013: 
    Added module docstring
    Removed option 'twiss'. See github issue #15
 - 3.1.2 vimaier 31th May 2013:
    Insterted checks for preconditions:
        more than 1 file needed
        at least one file with DPP=0.0 (adapted)
        ac file
    Extracted functions modelIntersect and intersect to Utilities.bpm
    Cleaned import section
    Restructured into parse_args, main, helper and main invocation
 - 3.1.3 ylevinse 6th of August 2013:
    Fix when analysing non-AC data
    Moved out some helper functions to superutils.py
    Some minor cleaning..
'''



import optparse
import os
import sys
import shutil
import subprocess
import math
import re

import __init__  # @UnusedImport init will include paths
import metaclass
import linreg
import Utilities.bpm

import superutils

#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
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
            metavar="<deltapScalingFactor>", default=1.0, type=float,dest="deltapScalingFactor")

    return parser.parse_args()


#===================================================================================================
# main()-function
#===================================================================================================
def main(options, args):
    files = superutils.get_filelist(options, args)
    
    if 2 > len(files):
        print >> sys.stderr, "Provide at least two files. Files:", str(files)
        sys.exit(1)

    if not os.path.isdir(options.output):
        os.makedirs(options.output)
    accel = options.accel
    technique = options.technique

    fileslist = {}

    for f_name in files:

        if(f_name[-3:]=='.gz'):
            datax = metaclass.twiss(f_name[:-3]+"_linx.gz")
            datay = metaclass.twiss(f_name[:-3]+"_liny.gz")
        else:
            datax = metaclass.twiss(f_name+"_linx")
            datay = metaclass.twiss(f_name+"_liny")

        dppx = datax.DPP*options.deltapScalingFactor      # Quick hack to be able to use old files with bad dpp input
        dppy = datay.DPP*options.deltapScalingFactor

        if dppx != dppy:
            raise ValueError("Discrepancy between horizontal and vertical => "+str(dppx)+" "+str(dppy))
        else:
            dpp = dppx/1.0

        if dpp not in fileslist:
            print "Adding dpp", dpp
            fileslist[dpp] = [f_name]
        else:
            fileslist[dpp].append(f_name)

    if 0 not in fileslist:
        raise ValueError("NO DPP=0.0. Provide at least one source file with DPP=0.0.\n"+
                         "In GUI you can change DPP in the table at the top of the analyses panel."+
                         "Click in the column 'delta p /p' on a value to change it.")

    getTunes(options, fileslist)

    madcreator(fileslist.keys(), options)
    print "All models are created"
    for dpp in fileslist:
        files = fileslist[dpp]
        rungetllm(options.output+"/twiss_"+str(dpp)+".dat", accel, technique, files, options, dpp)
    # The GUI wants the default files to have the names without _0.0
    copy_default_outfiles(options)

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
        metaclass.twiss(options.output+'/getbetax_free_'+str(dpp)+'.out')
        freeswitch = 1
    except IOError:
        print "WARNING: Could not open", options.output+'/getbetax_free_'+str(dpp)+'.out'
        freeswitch = 0

    for dpp in fileslist.keys():
        print "Loading driven data for ", dpp
        betx = metaclass.twiss(options.output+'/getbetax_'+str(dpp)+'.out')
        bety = metaclass.twiss(options.output+'/getbetay_'+str(dpp)+'.out')
        couple = metaclass.twiss(options.output+'/getcouple_'+str(dpp)+'.out')
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
        
        path_twissfile = superutils.get_twissfile(options)
        if not os.path.isfile(path_twissfile):
            print >> sys.stderr, "Twissfile does not exist:", path_twissfile
            sys.exit(1)
        modeld = metaclass.twiss(path_twissfile)

        #try:
        if freeswitch == 1:
            print "Loading free data"
            freeswitch = 1 #TODO: Actually senseless. Previous condition asserts ==1(vimaier)
            print 'getbetax_free_'+str(dpp)+'.out'
            betxf = metaclass.twiss(options.output+'/getbetax_free_'+str(dpp)+'.out')
            betyf = metaclass.twiss(options.output+'/getbetay_free_'+str(dpp)+'.out')
            couplef = metaclass.twiss(options.output+'/getcouple_free_'+str(dpp)+'.out')
            betalistxf[dpp] = betxf
            betalistyf[dpp] = betyf
            couplelistf[dpp] = couplef
            listxf.append(betxf)
            listyf.append(betyf)
            listcf.append(couplef)
            modelf = modeld
            path_ac_file = path_twissfile.replace(".dat","_ac.dat")
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
    fileobj = chromFileWriter('beta', options.output+"/chrombetax.out", 'H')

    bpms = Utilities.bpm.intersect(listx)
    bpms = Utilities.bpm.model_intersect(bpms, modeld)
    dolinregbet(fileobj, fileslist.keys(), betalistx, bpms, "H", zerobx, modeld)
    del fileobj

    #V
    fileobj = chromFileWriter('beta', options.output+"/chrombetay.out", 'V')

    bpms = Utilities.bpm.intersect(listy)
    bpms = Utilities.bpm.model_intersect(bpms, modeld)
    dolinregbet(fileobj, fileslist.keys(), betalisty, bpms, "V", zeroby, modeld)
    del fileobj

    print "Driven beta finished"

    #
    # driven coupling
    #
    print "Driven coupling"
    
    fileobj = chromFileWriter('coupling', options.output+"/chromcoupling.out", '')

    bpms = Utilities.bpm.intersect(listc)
    bpms = Utilities.bpm.model_intersect(bpms, modeld)

    dolinregCoupling(couplelist, bpms, fileslist.keys(), fileobj)
    del fileobj

    print "Driven coupling finished"

    if freeswitch == 1:
        #
        # free beta
        #
        print "Free beta"
        #H
        fileobj = chromFileWriter('beta', options.output+"/chrombetax_free.out", 'H')

        bpms = Utilities.bpm.intersect(listxf)
        bpms = Utilities.bpm.model_intersect(bpms, modelf)
        dolinregbet(fileobj, fileslist.keys(), betalistxf, bpms, "H", zerobxf, modelf)

        #V
        fileobj = chromFileWriter('beta', options.output+"/chrombetay_free.out", 'V')

        bpms = Utilities.bpm.intersect(listyf)
        bpms = Utilities.bpm.model_intersect(bpms, modelf)
        dolinregbet(fileobj, fileslist.keys(), betalistyf, bpms, "V", zerobyf, modelf)

        print "Free beta finished"

        #
        # free coupling
        #
        print "Free coupling"

        fileobj = chromFileWriter('coupling', options.output+"/chromcoupling_free.out", '')

        bpms = Utilities.bpm.intersect(listcf)
        bpms = Utilities.bpm.model_intersect(bpms, modelf)

        dolinregCoupling(couplelistf, bpms, fileslist.keys(), fileobj)

        print "Free coupling finished"


#===================================================================================================
# helper-functions
#===================================================================================================
def check_input(options, args):
    files = superutils.get_filelist(options, args)
    if len(files) == 0:
        raise SyntaxError("You need to define at least one file input")
    for f_name in files:
        if not os.path.isfile(f_name) and not os.path.isfile(f_name+'.gz'):
            raise ValueError(f_name+' does not exist')


def madcreator(dpps, options):
    '''
    :param dpps: list of dp/p to create model for
    :param options: dictionary of options from input arguments
    '''
    madfile = options.brc+"/MODEL/LHCB/model/"

    linesmad = open(madfile+"/job.twiss_chrom.madx.macro","r").read()

    # creating the DPP
    dppstring = ''
    dppstring_ac = ''
    for dpp in dpps:
        if not os.path.exists(options.output+'/twiss_'+str(dpp)+'.dat'):
            dppstring = dppstring+'twiss, chrom,sequence='+options.accel+', deltap='+str(dpp)+', file="'+options.output+'/twiss_'+str(dpp)+'.dat";\n'
            # if the model has twiss_ac.dat:
            if os.path.exists(superutils.get_twissfile(options)[:-4]+"_ac.dat"): # this is only correct as long as the filenames are <filename>_ac.dat and <filename>.dat!
                dppstring_ac = dppstring_ac+'twiss, chrom,sequence='+options.accel+', deltap='+str(dpp)+', file="'+options.output+'/twiss_'+str(dpp)+'_ac.dat";\n'
            else: # do not create ac file if we don't have ac in our original model..
                dppstring_ac = ''

    if not dppstring:
        print "No need to run madx"
        return 0

    DPP = dppstring
    DP_AC_P = dppstring_ac
    ACCEL=options.accel
    if options.accel == 'LHCB1':
        BEAM = 'B1'
    elif options.accel == 'LHCB2':
        BEAM = 'B2'
    else:
        print "WARNING: Could not decide what BEAM should be"
    QX = options.qx
    QY = options.qy
    QDX = options.qdx
    QDY = options.qdy
    QMX = int(options.qx*1000000)
    QMY = int(options.qy*1000000)
    STOP = '!'

    for testpath in [options.output, os.path.dirname(options.twissfile)]:
        _tmpmod = os.path.join(testpath, 'modifiers.madx')
        if os.path.isfile(_tmpmod):
            print "INFO: Using", _tmpmod
            MODIFIERS = _tmpmod
            break

    print "Creating madx"
    filetoprint = open(options.output+"/job.chrom.madx", "w")

    #changing variables
    filetoprint.write(linesmad % locals())

    filetoprint.close()
    
    print "Running madx"
    process = subprocess.Popen(options.madx+' < '+options.output+'/job.chrom.madx',
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

def filenames(dpp=''):
    '''
    Returns list of available file names
    for the given dpp.

    Example: varnames('0.0')

    :param dpp: [str] dpp appendix to list of files
    '''
    ret = []
    for fname in os.listdir(options.output):
        if dpp:
            extra = '_'+dpp
        else:
            extra = ''
        # thank you for making this easy..
        if re.search("get[A-Za-z]*[0-9]*[xy]*(_free[2]*)*"+extra+".out", fname): 
            ret.append(fname)
    return ret


def rungetllm(twiss_filename, accel, technique, files, options, dpp):
    '''
    Running GetLLM...
    '''
    import GetLLM

    print "Will run getllm for ", dpp #, command

    GetLLM.main(outputpath=options.output,
            files_to_analyse=','.join(files),
            model_filename=twiss_filename,
            accel=accel,
            TBTana=technique)
    print "GetLLM finished"

    for fname in filenames():
        f_identifier = fname.strip('.out')
        shutil.move(options.output+'/'+fname, options.output+'/'+f_identifier+'_'+str(dpp)+'.out')


def copy_default_outfiles(options):
    for fname in filenames('0.0'):
        f_identifier = fname.strip('_0.0.out')
        shutil.copy(options.output+'/'+fname, options.output+'/'+f_identifier+'.out')


##### for chromatic

def dolinregbet(fileobj, listx, listy, bpms, plane, zero, twiss):
    '''
    Calculates stuff and writes to the file in a table
    Closes the file afterwards

    :param fileobj: chromFileWriter for output table
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


def get_f( couplelist, dpplist, bpm_name, value):
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


def dolinregCoupling(couplelist, bpms, dpplist, fileobj):
    '''
    linreg for chromatic coupling

    Writes to fileobj the chromatic coupling.
    f1001, f1010 derivatives wrt dp/p, and errors.
    '''
    for bpm in bpms:
        name = bpm[1]
        sloc = bpm[0]

        chr_f1001r, chr_err_f1001r = get_f( couplelist, dpplist, name, 'F1001R')
        chr_f1001i, chr_err_f1001i = get_f( couplelist, dpplist, name, 'F1001I')
        chr_f1010r, chr_err_f1010r = get_f( couplelist, dpplist, name, 'F1010R')
        chr_f1010i, chr_err_f1010i = get_f( couplelist, dpplist, name, 'F1010I')

        mdl_chr_f1001r, mdl_chr_err_f1001r = get_f( couplelist, dpplist, name, 'MDLF1001R')
        mdl_chr_f1001i, mdl_chr_err_f1001i = get_f( couplelist, dpplist, name, 'MDLF1001I')
        mdl_chr_f1010r, mdl_chr_err_f1010r = get_f( couplelist, dpplist, name, 'MDLF1010R')
        mdl_chr_f1010i, mdl_chr_err_f1010i = get_f( couplelist, dpplist, name, 'MDLF1010I')

        fileobj.writeLine(locals().copy())
        

def getTunes(options,fileslist):
    '''
    Reads in the driven tunes from the
    file with dpp=0
    Reads in the model tunes from the 
    twiss model (twiss.dat)
    Appends the attributes to options.
    
    :param options: options from parse_args
    :param fileslist: dictionary of files, dpp used as key
    :raise ValueError: If fileslist[0] does not exist
    '''
    if fileslist[0][0][-3:] == '.gz':
        fname=fileslist[0][0][:-3]
        end='.gz'
    else:
        fname=fileslist[0][0]
        end=''
    tw_x = metaclass.twiss(fname+'_linx'+end)
    tw_y = metaclass.twiss(fname+'_liny'+end)
    tw = metaclass.twiss(superutils.get_twissfile(options))

    qdx,qdy = tw_x.TUNEX[0],tw_y.TUNEY[0]
    qx,qy = tw.Q1%1,tw.Q2%1

    setattr(options, "qx", qx)
    setattr(options, "qy", qy)
    setattr(options, "qdx", qdx)
    setattr(options, "qdy", qdy)


class chromFileWriter:
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


#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":

    options,args = parse_args()
    check_input(options, args)
    main(options, args)
