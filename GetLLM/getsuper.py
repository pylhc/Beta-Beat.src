'''
Created on ??/??/??

@author: ?

@version: 0.0.3

TODO: Description
What getsuper essentially does is run GetLLM on files with different dp/p and then afterwards 
interpolate the results to see how the functions vary with dp/p.

To run getsuper.py you need several source files(sdds files which stated in a comma separated string
in option -f/--files) with different DPP(delta_p/p). 
At least one of the source files must have DPP=0.0 .
Hint: You can change DPP in the GUI application in the Analysis panel. Change the entries in the
column 'dp/d' in the table at the top.

Further you need AC dipole in your model. If your twissfile is 'Twiss.dat' you need to provide the 
file 'Twiss_ac.dat' in the same directory. 


Change history:

 - 0.0.2 vimaier 28th May 2013: 
    Added module docstring
    Removed option 'twiss'. See github issue #15
 - 0.0.3 vimaier 31th May 2013:
    Insterted checks for preconditions:
        more than 1 file needed
        at least one file with DPP=0.0 (adapted)
        ac file
'''



###### imports
from optparse import OptionParser
import os,sys,shutil,subprocess,re
if "/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/" not in sys.path: # add internal path for python scripts to current environment (tbach)
    sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')
import math
from metaclass import twiss
import linreg

##
# YIL changes v 3.1:
#  - Cleaned macro writer in madcreator
#  - modifiers.madx should be in options.output
#


class chromFileWriter:
    def __init__(self,ftype,fname,plane,overwrite=True):
        '''
        :param ftype: string, 'beta' or 'coupling'
        :param fname: string, name of file
        :param plane: "H" or "V"
        :param overwrite: Overwrite file if it already exist
        '''
        self.fstream=file(fname,'w')
        self._clen=19 # column length..

        if plane.upper()=="H":
            plane="X"
        elif plane.upper()=="V":
            plane="Y"

        if os.path.isfile(fname) and not overwrite:
            raise ValueError("Cannot overwrite file "+fname)

        betacolumns= ['name', 'sloc',  'dbb', 'dbberr', 'da', 'daerr', 
                    'w', 'werr','wm', 'werrm', 'wmo', 'phi', 'phierr','phim','phierrm','pmo',
                    'A','Aerr','Am','Aerrm','B','Berr','Bm','Berrm', 'dbbm','dbberrm',
                    'dam','daerrm']
        couplecolumns=['name', 'sloc', 'chr_f1001r', 'chr_err_f1001r', 'chr_f1001i',
            'chr_f1001r', 'chr_err_f1001r', 'chr_f1001i', 'chr_err_f1001i',
            'chr_f1010r', 'chr_err_f1010r', 'chr_f1010i', 'chr_err_f1010i',
            'mdl_chr_f1001r', 'mdl_chr_err_f1001r', 'mdl_chr_f1001i', 'mdl_chr_err_f1001i',
            'mdl_chr_f1010r', 'mdl_chr_err_f1010r', 'mdl_chr_f1010i', 'mdl_chr_err_f1010i',
          ]

        headnames={
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
        headtypes={'name':'%s'}
        for key in headnames:
            # for all others we use '%le'...
            if key not in headtypes:
                headtypes[key]='%le'

        if ftype.lower().strip()=='beta':
            self.ftype='beta'
            self.columns=betacolumns
        elif ftype.lower().strip()=='coupling':
            self.ftype='coupling'
            self.columns=couplecolumns
        else:
            raise ValueError("ftype %s not understood" % (ftype))
        
        self.head=[headnames[c] % locals() for c in self.columns]
        headcount=[i+1 for i in xrange(len(self.head))]
        self.types=[headtypes[c] for c in self.columns]
        
        self.head[0]='* '+(self.head[0].rjust(self._clen-3))
        headcount[0]='# '+(str(headcount[0]).rjust(self._clen-3))
        self.types[0]='$ '+(self.types[0].rjust(self._clen-3))
        
        self._write_list(self.head)
        # According to SL-CO-Note-91-32, this is allowed...
        self._write_list(headcount)
        self._write_list(self.types)

    def writeLine(self,data):
        '''
        Write one line. 
        :param data: Dictionary of data columns to write
        '''
        tmp_line=[data[c] for c in self.columns]
        self._write_list(tmp_line)

    def _write_list(self,tmp_list):
        line=''
        for l in tmp_list:
            line+=str(l).rjust(self._clen-1)+' '
        self.fstream.write(line[:-1]+'\n')

def parse_args():
    ###### optionparser
    usage = "usage: %prog [options] sdds-file1 [sdds-file2 ...]"
    parser = OptionParser(usage)
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
    # By default we take the path from where getsuper_dev.py is ran from..
    parser.add_option("-b", "--beta",
            help="Path to Beat-Beat.src folder",
            metavar="<path>", default=os.path.dirname(__file__)+'/../', dest="brc")
    parser.add_option("-t", "--algorithm",
            help="Which algorithm to use (SUSSIX/SVD)",
            metavar="ALGORITHM", default="SUSSIX", dest="technique")
    parser.add_option("-a", "--accel",
            help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
            metavar="ACCEL", default="LHCB1",dest="accel")
    parser.add_option("", "--llm_version",
            help="Run with specific version of GetLLM.py",
            metavar="<version>", default=None,dest="llm_version")

    parser.add_option("-d", "--deltapScalingFactor",
            help="Scaling factor for deltap, remember final value must be in MAD units",
            metavar="<deltapScalingFactor>", default=1.0, type=float,dest="deltapScalingFactor")

    return parser.parse_args()



def check_input(options,args):
    files=get_filelist(options,args)
    if len(files)==0:
        raise SyntaxError("You need to define at least one file input")
    for f in files:
        if not os.path.isfile(f) and not os.path.isfile(f+'.gz'):
            raise ValueError(f+' does not exist')


def get_filelist(options,args):
    '''
    Returns list of files to be analysed
    '''
    if options.files:
        files=[f.strip() for f in options.files.split(',')]
        files.extend(args)
        return files
    return args

def get_twissfile(options):
    '''
    Returns the full path to
    the twiss file
    '''
    if options.twissfile:
        return options.twissfile
    if os.path.isfile(options.twissfile+'/twiss.dat'):
        return options.twissfile+'/twiss.dat'
    if os.path.isfile(options.twissfile+'/twiss.dat.gz'):
        return options.twissfile+'/twiss.dat.gz'
    # did not find any file..
    raise ValueError("Could not find twissfile! "+options.twissfile)

## ############
#functions
## ############


#####
def madcreator(dpps,options):
    '''
    :param dpps: list of dp/p to create model for
    :param options: dictionary of options from input arguments
    '''

    madfile=options.brc+"/MODEL/LHCB/model/"

    linesmad=open(madfile+"/job.twiss_chrom.madx.macro","r").read()


    # creating the DPP
    dppstring=''
    dppstring_ac=''
    for dpp in dpps:
        if not os.path.exists(options.output+'/twiss_'+str(dpp)+'.dat'):
            dppstring=dppstring+'twiss, chrom,sequence='+options.accel+', deltap='+str(dpp)+', file="'+options.output+'/twiss_'+str(dpp)+'.dat";\n'
            dppstring_ac=dppstring_ac+'twiss, chrom,sequence='+options.accel+', deltap='+str(dpp)+', file="'+options.output+'/twiss_'+str(dpp)+'_ac.dat";\n'

    if not dppstring:
        print "No need to run madx"
        return 0

    DPP=dppstring
    DP_AC_P=dppstring_ac
    ACCEL=options.accel
    if options.accel=='LHCB1':
        BEAM='B1'
    elif options.accel=='LHCB2':
        BEAM='B2'
    else:
        print "WARNING: Could not decide what BEAM should be"
    QX=options.qx
    QY=options.qy
    QDX=options.qdx
    QDY=options.qdy
    QMX=int(options.qx*1000000)
    QMY=int(options.qy*1000000)
    STOP='!'

    for testpath in [options.output,os.path.dirname(options.twissfile)]:
        _tmpmod=os.path.join(testpath,'modifiers.madx')
        if os.path.isfile(_tmpmod):
            print "INFO: Using",_tmpmod
            MODIFIERS=_tmpmod
            break

    print "Creating madx"
    filetoprint=open(options.output+"/job.chrom.madx","w")


    #changing variables
    filetoprint.write(linesmad % locals())

    filetoprint.close()
    print "Running madx"
    if subprocess.call(options.madx+' < '+options.output+'/job.chrom.madx',shell=True):
        raise ValueError("Mad-X failed")

def append(files):
    return ','.join(files)

def filenames(dpp=''):
    '''
    Returns list of available file names
    for the given dpp.

    Example: varnames('0.0')

    :param dpp: [str] dpp appendix to list of files
    '''
    ret=[]
    for fname in os.listdir(options.output):
        if dpp:
            extra='_'+dpp
        else:
            extra=''
        # thank you for making this easy..
        if re.search("get[A-Za-z]*[0-9]*[xy]*(_free[2]*)*"+extra+".out",fname): 
            ret.append(fname)
    return ret

def rungetllm(twissfile,accel,technique,files,options,dpp):
    '''
    Running GetLLM...
    '''
    if options.llm_version:
        print "GetLLM_V"+options.llm_version+" as GetLLM"
        exec("import GetLLM_V"+options.llm_version+" as GetLLM")
    else:
        import GetLLM

    print "Will run getllm for ",dpp #, command

    GetLLM.main(outputpath=options.output,
            files_to_analyse=append(files),
            twiss_model_file=twissfile,
            accel=accel,
            TBTana=technique)
    print "GetLLM finished"

    for fname in filenames():
            v=fname.strip('.out')
            shutil.move(options.output+'/'+fname,options.output+'/'+v+'_'+str(dpp)+'.out')

def copy_default_outfiles(options):
    for fname in filenames('0.0'):
        v=fname.strip('_0.0.out')
        shutil.copy(options.output+'/'+fname,options.output+'/'+v+'.out')

##### for chromatic
# model intersect
def modelIntersect(expbpms, model):

    bpmsin=[]
    for bpm in expbpms:
        try:
            model.indx[bpm[1].upper()]
            bpmsin.append(bpm)
        except:
            print bpm, "Not in Model"
    if len(bpmsin)==0:
        print "Zero intersection of Exp and Model"
        print "Please, provide a good Dictionary"
        print "Now we better leave!"
        sys.exit()
    return bpmsin


#intersect
def intersect(ListOfFile):
    '''Pure intersection of all bpm names in all files '''
    if len(ListOfFile)==0:
        print "Nothing to intersect!!!!"
        sys.exit()
    z=ListOfFile[0].NAME
    for b in ListOfFile:
        z=filter(lambda x: x in z   , b.NAME)
    #SORT by S
    result=[]
    x0=ListOfFile[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))

    result.sort()
    return result

def dolinregbet(fileobj,listx,listy,bpms,plane,zero,twiss):
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
        name=bpm[1]
        sloc=bpm[0]
        indx=[]
        b=[]
        a=[]
        bm=[]
        am=[]
        if "H" in plane:
            beta0=zero.BETX[zero.indx[name]]
            alfa0=zero.ALFX[zero.indx[name]]
            alfa0err=zero.STDALFX[zero.indx[name]]

            beta0m=twiss.BETX[twiss.indx[name]]
            alfa0m=twiss.ALFX[twiss.indx[name]]

            wmo=twiss.WX[twiss.indx[name]]
            pmo=twiss.PHIX[twiss.indx[name]]
        else:

            beta0=zero.BETY[zero.indx[name]]
            alfa0=zero.ALFY[zero.indx[name]]
            alfa0err=zero.STDALFY[zero.indx[name]]

            beta0m=twiss.BETY[twiss.indx[name]]
            alfa0m=twiss.ALFY[twiss.indx[name]]

            wmo=twiss.WY[twiss.indx[name]]
            pmo=twiss.PHIY[twiss.indx[name]]
        for dpp in listx:
            _file=listy[dpp]
            ix=_file.indx[name]
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

        bfit=linreg.linreg(listx, b)
        afit=linreg.linreg(listx, a)

        bfitm=linreg.linreg(listx, bm)
        afitm=linreg.linreg(listx, am)

        # measurement
        dbb=bfit[0]/beta0
        dbberr=bfit[3]/beta0
        da=afit[0]
        daerr=afit[3]
        A=dbb
        Aerr=dbberr
        B=da-alfa0*dbb
        Berr=math.sqrt(daerr**2 + (alfa0err*dbb)**2 + (alfa0*dbberr)**2)
        w=0.5*math.sqrt(A**2+B**2)
        werr=0.5*math.sqrt( (Aerr*A/w)**2 + (Berr*B/w)**2  )
        phi=math.atan2(B,A)/2./math.pi
        phierr=1./(1.+(A/B)**2)*math.sqrt( (Aerr/B)**2 + (A/B**2*Berr)**2)/2./math.pi

        #model
        dbbm=bfitm[0]/beta0m
        dbberrm=bfitm[3]/beta0m
        dam=afitm[0]
        daerrm=afitm[3]
        Am=dbbm
        Aerrm=dbberrm
        Bm=dam-alfa0m*dbbm
        Berrm=math.sqrt(daerrm**2 + (alfa0m*dbberrm)**2)
        wm=0.5*math.sqrt(Am**2+Bm**2)
        werrm=0.5*math.sqrt( (Aerrm*Am/wm)**2 + (Berrm*Bm/wm)**2  )
        phim=math.atan2(Bm,Am)/2./math.pi
        phierrm=1./(1.+(Am/Bm)**2)*math.sqrt( (Aerrm/Bm)**2 + (Am/Bm**2*Berrm)**2)/2./math.pi

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
    lst=[]
    x=[]
    for dpp in dpplist:
        x.append(dpp)
        couplefile=couplelist[dpp]
        lst.append( getattr(couplefile, value)[ couplefile.indx[ bpm_name]])
    lreg=linreg.linreg(x,lst)
    return lreg[0],lreg[3]

def dolinregCoupling(couplelist,bpms,dpplist,fileobj):
    '''
    linreg for chromatic coupling

    Writes to fileobj the chromatic coupling.
    f1001, f1010 derivatives wrt dp/p, and errors.
    '''
    for bpm in bpms:

        name=bpm[1]
        sloc=bpm[0]


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
    if fileslist[0][0][-3:]=='.gz':
        fname=fileslist[0][0][:-3]
        end='.gz'
    else:
        fname=fileslist[0][0]
        end=''
    tw_x=twiss(fname+'_linx'+end)
    tw_y=twiss(fname+'_liny'+end)
    tw=twiss(get_twissfile(options))

    qdx,qdy=tw_x.TUNEX[0],tw_y.TUNEY[0]
    qx,qy=tw.Q1%1,tw.Q2%1

    setattr(options,"qx",qx)
    setattr(options,"qy",qy)
    setattr(options,"qdx",qdx)
    setattr(options,"qdy",qdy)


def main(options,args):

    ## ##############
    #   main
    ## ##############

    files=get_filelist(options,args)
    
    if 2 > len(files):
        print >> sys.stderr,"Provide at least two files. Files:",str(files)
        sys.exit(1)

    if not os.path.isdir(options.output):
        os.makedirs(options.output)
    accel=options.accel
    technique=options.technique

    fileslist={}

    for f in files:

        if(f[-3:]=='.gz'):
            datax=twiss(f[:-3]+"_linx.gz")
            datay=twiss(f[:-3]+"_liny.gz")
        else:
            datax=twiss(f+"_linx")
            datay=twiss(f+"_liny")

        dppx=datax.DPP*options.deltapScalingFactor      # Quick hack to be able to use old files with bad dpp input
        dppy=datay.DPP*options.deltapScalingFactor

        if dppx!=dppy:
            raise ValueError("Discrepancy between horizontal and vertical => "+str(dppx)+" "+str(dppy))
        else:
            dpp=dppx/1.0

        if dpp not in fileslist:
            print "Adding dpp",dpp
            fileslist[dpp]=[f]
        else:
            fileslist[dpp].append(f)


    if 0 not in fileslist:
        raise ValueError("NO DPP=0.0. Provide at least one source file with DPP=0.0")

    getTunes(options,fileslist)

    madcreator(fileslist.keys(),options)
    print "All models are created"
    for dpp in fileslist:
        files=fileslist[dpp]
        rungetllm(options.output+"/twiss_"+str(dpp)+".dat",accel,technique,files,options,dpp)
    # The GUI wants the default files to have the names without _0.0
    copy_default_outfiles(options)


    ##adding data
    betalistx={}
    betalisty={}
    couplelist={}
    betalistxf={}
    betalistyf={}
    couplelistf={}

    listx=[]
    listxf=[]
    listy=[]
    listyf=[]
    listc=[]
    listcf=[]

    try:
        twiss(options.output+'/getbetax_free_'+str(dpp)+'.out')
        freeswitch=1
    except:
        print "WARNING: Could not open",options.output+'/getbetax_free_'+str(dpp)+'.out'
        freeswitch=0


    for dpp in fileslist.keys():
        print "Loading driven data for ",dpp
        betx=twiss(options.output+'/getbetax_'+str(dpp)+'.out')
        bety=twiss(options.output+'/getbetay_'+str(dpp)+'.out')
        couple=twiss(options.output+'/getcouple_'+str(dpp)+'.out')
        #couple=twiss(options.output+'/getbetay_'+str(dpp)+'.out')
        betalistx[dpp]=betx
        betalisty[dpp]=bety
        couplelist[dpp]=couple

        if float(dpp)==0.0:
            zerobx=betx
            zeroby=bety

        listx.append(betx)
        listy.append(bety)
        listc.append(couple)
        
        path_twissfile = get_twissfile(options)
        if not os.path.isfile(path_twissfile):
            print >> sys.stderr, "Twissfile does not exist:",path_twissfile
            sys.exit(1)
        modeld = twiss(path_twissfile)

        #try:
        if freeswitch==1:
            print "Loading free data"
            freeswitch=1
            print 'getbetax_free_'+str(dpp)+'.out'
            betxf=twiss(options.output+'/getbetax_free_'+str(dpp)+'.out')
            betyf=twiss(options.output+'/getbetay_free_'+str(dpp)+'.out')
            couplef=twiss(options.output+'/getcouple_free_'+str(dpp)+'.out')
            betalistxf[dpp]=betxf
            betalistyf[dpp]=betyf
            couplelistf[dpp]=couplef
            listxf.append(betxf)
            listyf.append(betyf)
            listcf.append(couplef)
            modelf = modeld
            path_ac_file = path_twissfile.replace(".dat","_ac.dat")
            if not os.path.isfile(path_ac_file):
                print >> sys.stderr, "Ac file does not exist:",path_ac_file,"\nIn GUI check 'Ac dipole' box to create a model with ac dipole."
                sys.exit(1)
            modeld=twiss(path_ac_file)
            if float(dpp)==0.0:
                zerobxf=betalistxf[dpp]
                zerobyf=betalistyf[dpp]


        #except:
        #       print "No free data"

    #
    # driven beta
    #

    print "Driven beta"

    #H
    fileobj=chromFileWriter('beta',options.output+"/chrombetax.out",'H')

    bpms=intersect(listx)
    bpms=modelIntersect(bpms,modeld)
    dolinregbet(fileobj,fileslist.keys(),betalistx,bpms,"H",zerobx,modeld)
    del fileobj

    #V
    fileobj=chromFileWriter('beta',options.output+"/chrombetay.out",'V')

    bpms=intersect(listy)
    bpms=modelIntersect(bpms,modeld)
    dolinregbet(fileobj,fileslist.keys(),betalisty,bpms,"V",zeroby,modeld)
    del fileobj

    print "Driven beta finished"

    #
    # driven coupling
    #
    print "Driven coupling"
    
    fileobj=chromFileWriter('coupling',options.output+"/chromcoupling.out",'')

    bpms=intersect(listc)
    bpms=modelIntersect(bpms,modeld)

    dolinregCoupling(couplelist,bpms,fileslist.keys(),fileobj)
    del fileobj

    print "Driven coupling finished"

    if freeswitch==1:
        #
        # free beta
        #
        print "Free beta"
        #H
        fileobj=chromFileWriter('beta',options.output+"/chrombetax_free.out",'H')

        bpms=intersect(listxf)
        bpms=modelIntersect(bpms,modelf)
        dolinregbet(fileobj,fileslist.keys(),betalistxf,bpms,"H",zerobxf,modelf)

        #V
        fileobj=chromFileWriter('beta',options.output+"/chrombetay_free.out",'V')

        bpms=intersect(listyf)
        bpms=modelIntersect(bpms,modelf)
        dolinregbet(fileobj,fileslist.keys(),betalistyf,bpms,"V",zerobyf,modelf)

        print "Free beta finished"

        #
        # free coupling
        #
        print "Free coupling"

        fileobj=chromFileWriter('coupling',options.output+"/chromcoupling_free.out",'')

        bpms=intersect(listcf)
        bpms=modelIntersect(bpms,modelf)

        dolinregCoupling(couplelistf,bpms,fileslist.keys(),fileobj)


        print "Free coupling finished"


if __name__=="__main__":

    options,args = parse_args()
    check_input(options,args)
    main(options,args)
