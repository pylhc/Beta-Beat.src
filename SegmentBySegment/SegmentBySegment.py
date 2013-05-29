'''
Created on 11/09/09

@author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)

@version: 0.34

TODO: Description

Change history:
#  !=> SegementBySegment_0.0.py : - Construction of the program (11/09/09)
#  !=> SegementBySegment_0.1.py : - Adding the dispersion to the program (23/09/09)
#                                 - Make distuinsh between BPMs and instuments (24/09/09)
#                                 - Adding output tables (sbs...out) (21/10/09)
#                                 - Adding coupling propogation (21/10/09)
#                                 - changing output when segment is start and end (24/11/09)
#  !=> SegementBySegment_0.22.py :- Fixing some small bugs
#                                 - Added Instruments
#
#  !=> SegementBySegment_0.23.py :- Adding summary list for instruments
#                                 - Applying better cuts for beta
#  !=> SegementBySegment_0.24.py :- Including one/two beam matching (two beam for IP) (7/12/09)
#    
#  !=> SegementBySegment_0.25.py :- Adding phase to the propagation (10/12/09)
#
#  !=> SegementBySegment_0.26.py :- General clean-up (19/01/10)
#                                 - Adding  Total phase to the propagation (19/01/10)
#                                 - Coupling propogation (20/01/10)
#
#  !=> SegementBySegment_0.27.py :- Fixing mistake of total phase (1/03/10)
#                                 - Include real and imaginary part in the output of coupling (1/03/10)
#
#  !=> SegementBySegment_0.28.py : -Coupling initial conditions added to MADX segment  24 March 2010
#                                  -New file getfterms.py required to convert MADX C matrix to observable f terms
#        4 April 2010              -Introducing the modifiers.madx     
#                                   file, it should be in twisspath or an 
#                                   empty file will be created.             
#                                   The MAD mask file has also been 
#                                   updated
#  !=> SegementBySegment_0.29.py : - Update for IP2  beam1 and IP8 beam 2... problem with S coordinate
#  !=>  SegementBySegment_0.30.py : - New filtering... 
#  !=>  SegementBySegment_0.31.py : - New and easier filtering
#  !=>  SegementBySegment_0.31.py : - Adjusting filtering
#                                   - Adding C-matrix,angle,waist calculation,...
#                                   - Adding IP calculation from getllm... if BPM is not found you take another.
#   Added system path to the beginning to avoid having to define it in the shell
#
# 0.33 20120820 (tbach):
# - changed savepath = options.SAVE to savepath = options.SAVE + "/"
# - fixed errors and warnings from static code analysis (unused variables, not defined variables, variable names overriding bultin reserved key words, broken indentation) 
# - fixed mixed indentation (hint: if you want to a diff which ignores the changed indentation, use sdiff -W -s
# - deleted code which is commented out
 - 0.34 vimaier 28th May 2013: 
    Restructured code into parse_args(),main(),helper-functions and main invocation
    Adapted docstrings.
    Removed unused parameters from getIPfrompara
    Removed unused function getTwiss
    Previous version can be found in git commit eca47b2d06ab7eab911c830f30dc222a8bcd91a0

'''


import os
import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import optparse
from math import sqrt,cos,sin,pi, tan

try:
    from metaclass import twiss
except:
    from metaclass25 import twiss



#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
    ''' Parses the arguments, checks for valid input and returns tupel '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                    help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                    metavar="ACCEL", default="LHCB1",dest="accel")
    parser.add_option("-f", "--path", # assumes that output is same as input
                    help="Path to measurement files",
                    metavar="PATH", default="./", dest="path")
    parser.add_option("-i", "--path2", # assumes that output is same as input
                    help="Path to second measurement files",
                    metavar="PATH2", default="./", dest="path2")
    parser.add_option("-s", "--start",
                    help="give start,endbpm,name (multiple allowed) eg: start1,end1,name1,start2,end2,name2,...",
                    metavar="SEGF", default="./", dest="segf")
    parser.add_option("-t", "--twiss",
                    help="basic twiss file, the modifiers.madx is assumed to be in the same direcory",
                    metavar="TWISS", default="./", dest="twiss")
    parser.add_option("-r", "--response",
                    help="switch calculate corrections 0 no, 1 yes",
                    metavar="COR", default="0", dest="cor")
    parser.add_option("-g", "--graphic",
                    help="choice between graphic or table output 0=graphic and 1=table output",
                    metavar="GRA", default="1", dest="gra")
    parser.add_option("-p", "--save",  
                    help="Output path",
                    metavar="SAVE", default="./", dest="SAVE")
    parser.add_option("-m", "--mad", # assumes that output is same as input
                    help="mad link",
              metavar="mad", default="", dest="mad")
    parser.add_option("-b", "--bbsource", # assumes that output is same as input
                    help="beta beat source",
                    metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="bb")
    parser.add_option("-x", "--take", # take or create mad input, default should be 0 for creating
                    help="take or create madx 0/1",
                    metavar="mad", default="0", dest="madpass")
    parser.add_option("-c", "--cuts", 
                    help="cut on error of beta in percentage",
                    metavar="cuts", default="10", dest="cuts")
    parser.add_option("-w", "--w", # Path to Chromaticity functions
                        help="Path to  chromaticity functions, by default this is skiped",
                        metavar="wpath", default="0", dest="wpath")
    
    (options, args) = parser.parse_args()
    
    return options,args


#===================================================================================================
# main()-function
#===================================================================================================
def main(options):
    '''
    :Parameters:
        'options': Values
            Values instance with all options from OptionParser
    :Return: int
        0 if execution was successful otherwise !=0
    '''
    path=options.path
    ##########
    #
    # LOADING DATA
    #
    ##########
    try:
        filedatax=twiss(path+"/getbetax_free.out")
    except:
        filedatax=twiss(path+"/getbetax.out")
    betxtwiss=filedatax
    QXX=filedatax.Q1
    QYY=filedatax.Q2
    
    try:
        filedatay=twiss(path+"/getbetay_free.out")
    except:
        filedatay=twiss(path+"/getbetay.out")    
    
    betytwiss=filedatay
    filedataxA=twiss(path+"/getampbetax.out")
    filedatayA=twiss(path+"/getampbetay.out")
    try:
        filephasex=twiss(path+"/getphasex_free.out")
    except:
        filephasex=twiss(path+"/getphasex.out")    
    
    try:
        filephasey=twiss(path+"/getphasey.out")
    except:
        filephasey=twiss(path+"/getphasey_free.out")
    try:
        filephasextot=twiss(path+"/getphasetotx_free.out")
    except:
        filephasextot=twiss(path+"/getphasetotx.out")
    
    try:
        filephaseytot=twiss(path+"/getphasetoty_free.out")
    except:
        filephaseytot=twiss(path+"/getphasetoty.out")
            
    ### check if dispersion exist
    try:
        filedx=twiss(path+"/getDx.out")
        filendx=twiss(path+"/getNDx.out")
        filedy=twiss(path+"/getDy.out")
        disp=1
        print "Disp files OK", filedx.DX[0], filedx.NAME[0]
    except:
        print "no dispersion files... will continue without taking into account dispersion"
        filendx=[]
        filedx=[]
        filedy=[]
        disp=0
    
    ### check if coupling exist
    try:
        filecouple=twiss(path+"/getcouple.out")
        filecoupleC=twiss(path+"/getcoupleterms.out")
        coupleswitch=1
        method="driven"
    except:
        print "no coupling file... will continue without taking into account coupling"
        filecouple=[]
        filecoupleC=[]
        coupleswitch=0
        fs=[0,0,0,0]
    
        
    try:
        filecouple=twiss(path+"/getcouple_free.out")
        method="_free"
        print "Free coupling found"
    
    except:
        print 
        #method=""
        #%method_label="driven"
        #filecouple=twiss(path+"/getcouple.out")
        #print "WARN: Free coupling not found!!"
    
    
    # checking if end element is BPM or instrument (IP,collimators)
    elementswitch=-1
    savepath = options.SAVE + "/"
    if not os.path.isdir(savepath):
        os.makedirs(savepath)
    list2run=options.segf.split(',')
    errorcut=float(options.cuts)
    accel=options.accel
    
    start={}
    end={}
    names=[]
    
    
    twissfile=options.twiss
    twisspath=os.path.dirname(twissfile)+'/'
    if twissfile=="./":
        twissfile=options.bb+"/MODEL/"+options.accel+"/nominal.opt/twiss_elements.dat"
    
    print "twissfile ",twissfile
    
    
    twisstwiss=twiss(twissfile)
    
    step=0# for selecting segments and elements
    for run in range(len(list2run)):
        run=run+step
    
        if run==(len(list2run)):
            break
    
        print list2run[run],run
        
        if "BPM" in list2run[run]: # segment
            
            start[list2run[run+2]]=list2run[run]
            end[list2run[run+2]]=list2run[run+1]
            names.append(list2run[run+2])
            step=2+step
        else: #element
            step=0+step
            names.append(list2run[run])
    
    
    ####### Difference between segment or instrument #######
    # For segment => Startbpm,endbpm,NameOfSegment
    # For instrument => instrument
    ####### 
    
    for namename in names:
    
        print namename
        ##
        # Getting the BPM's
        ##
        if start.has_key(namename) and end.has_key(namename):
            print "Segment has been choosen"
            elementswitch=0
            segment=[start[namename],end[namename]]
            startbpm,endbpm=filterandfind(betxtwiss,betytwiss,"null",segment,twisstwiss,errorcut)
    
        elif start.has_key(namename) or end.has_key(namename):
    
            print "Something strange ....Did you properly define the input ?"
            print "Like: BPM1,BPM2,ARC12,IP1"
            print "Structure must be for Segment => BPML,BPMR,NAME"
            print "For Instrument just name"
            sys.exit()
        else:
            print "Element has been choosen"
            elementswitch=1
            startbpm,endbpm=filterandfind(betxtwiss,betytwiss,namename,[],twisstwiss,errorcut)
    
        passs=1 # why?
    
        #gathering data
        hor=[betxtwiss.BETX[betxtwiss.indx[startbpm]],
             sqrt(betxtwiss.STDBETX[betxtwiss.indx[startbpm]]**2+betxtwiss.ERRBETX[betxtwiss.indx[startbpm]]**2),
             betxtwiss.ALFX[betxtwiss.indx[startbpm]],
             sqrt(betxtwiss.ERRALFX[betxtwiss.indx[startbpm]]**2+betxtwiss.STDALFX[betxtwiss.indx[startbpm]]**2)]
        ver=[betytwiss.BETY[betytwiss.indx[startbpm]],
             sqrt(betytwiss.STDBETY[betytwiss.indx[startbpm]]**2+betytwiss.ERRBETY[betytwiss.indx[startbpm]]**2),
             betytwiss.ALFY[betytwiss.indx[startbpm]],
             sqrt(betytwiss.ERRALFY[betytwiss.indx[startbpm]]**2+betytwiss.STDALFY[betytwiss.indx[startbpm]]**2)]
        hore=[betxtwiss.BETX[betxtwiss.indx[endbpm]],
             sqrt(betxtwiss.STDBETX[betxtwiss.indx[endbpm]]**2+betxtwiss.ERRBETX[betxtwiss.indx[endbpm]]**2),
             betxtwiss.ALFX[betxtwiss.indx[endbpm]],
             sqrt(betxtwiss.ERRALFX[betxtwiss.indx[endbpm]]**2+betxtwiss.STDALFX[betxtwiss.indx[endbpm]]**2)]
        vere=[betytwiss.BETY[betytwiss.indx[endbpm]],
             sqrt(betytwiss.STDBETY[betytwiss.indx[endbpm]]**2+betytwiss.ERRBETY[betytwiss.indx[endbpm]]**2),
             betytwiss.ALFY[betytwiss.indx[endbpm]],
             sqrt(betytwiss.ERRALFY[betytwiss.indx[endbpm]]**2+betytwiss.STDALFY[betytwiss.indx[endbpm]]**2)]
    
        beta4plot=hor[2]
    
        if disp==1:
    
            try: #horizontal
                filedx.indx[startbpm]
                dxpass=1
            except:
                print "startbpm for horizontal dispersion not found"
                dxpass=0
    
            try: #vertical
                filedy.indx[startbpm]
                dypass=1
            except:
                print "startbpm for  vertical dispersion not found"
                dypass=0
    
            try: #horizontal
                filedx.indx[endbpm]
                dxpasse=1
            except:
                print "endbpm for horizontal dispersion not found"
                dxpasse=0
    
            try: #vertical
                filedy.indx[endbpm]
                dypasse=1
            except:
                print "endbpm for  vertical dispersion not found"
                dypasse=0
    
            if dxpass==1:
                dxx=filedx.DX[filedx.indx[startbpm]]
                dxp=filedx.DPX[filedx.indx[startbpm]]
                dxe=filedx.STDDX[filedx.indx[startbpm]]
            else:
                dxx=0
                dxp=0
                dxe=0
    
            if dypass==1:
                dyy=filedy.DY[filedy.indx[startbpm]]
                dyp=filedy.DPY[filedy.indx[startbpm]]
                dye=filedy.STDDY[filedy.indx[startbpm]]
            else:
                dyy=0
                dyp=0
                dye=0            
            
            dp=[dxx,
                dxp,
                dyy,
                dyp,
                dxe,
                dye]
    
        
            if dxpasse==1:
                dxx=filedx.DX[filedx.indx[endbpm]]
                dxp=filedx.DPX[filedx.indx[endbpm]]
                dxe=filedx.STDDX[filedx.indx[endbpm]]
            else:
                dxx=0
                dxp=0
                dxe=0
    
            if dypasse==1:
                dyy=filedy.DY[filedy.indx[endbpm]]
                dyp=filedy.DPY[filedy.indx[endbpm]]
                dye=filedy.STDDY[filedy.indx[endbpm]]
            else:
                dyy=0
                dyp=0
                dye=0    
                
            dpe=[dxx,
                dxp,
                dyy,
                dyp,
                dxe,
                dye]
        else:
            dp=[0,0,0,0,0,0]
            dpe=[0,0,0,0,0,0]
            print "No dispersion"
    
    
        if passs==1:#TODO: Pretty useless. passs will always be initialized with 1(vimaier)
    
    
            if coupleswitch==1:
                print "coupling,"
                try:
                    filecouple.indx[startbpm]
                    cpass=1
                except:
                    cpass=0
                    print "startbpm ",startbpm," not found in coupling measurement => values=0"
                    
                if cpass==1:
                    f1001r=filecouple.F1001R[filecouple.indx[startbpm]]
                    f1001i=filecouple.F1001I[filecouple.indx[startbpm]]
                    f1010r=filecouple.F1010R[filecouple.indx[startbpm]]
                    f1010i=filecouple.F1010I[filecouple.indx[startbpm]]
                    
                    f1001std=filecouple.FWSTD1[filecouple.indx[startbpm]]
                    f1010std=filecouple.FWSTD2[filecouple.indx[startbpm]]
                
                    ### add error
                else:
                    f1001r=0
                    f1001i=0
                    f1010r=0
                    f1010i=0
                    f1001std=0
                    f1010std=0
                    
                fs=[f1001r,f1001i,f1010r,f1010i,f1001std,f1010std]
    
            
    
            else:
                fs=[0,0,0,0,0,0]
                print "No coupling"
    
            #sys.exit()
                print "madpass", options.madpass        
            if str(options.madpass)=="0":
                print "Going to run4mad"
                run4mad(savepath,hor,ver,hore,vere,dp,dpe,startbpm,endbpm,namename, fs, options.path+"/",twisspath, method)
    
            else:
                runmad(savepath,namename)
                print "Just rerunning mad"
            reversetable(savepath,namename)
    
    
                    ###############################################
            #
            # =>>>>>> including new writer
            #
                    ###############################################
    
            # loading twiss
            errbetamin=twiss(savepath+'/twiss.b-.dat')
            
            errbetamax=twiss(savepath+'/twiss.b+.dat')
            
            erralfmin=twiss(savepath+'/twiss.a-.dat')
            erralfminb=twiss(savepath+'/twiss.a-_back.dat')
            
            erralfmax=twiss(savepath+'/twiss.a+.dat')
            erralfmaxb=twiss(savepath+'/twiss.a+_back.dat')
            
            errdmin=twiss(savepath+'/twiss.d-.dat')
            errdminb=twiss(savepath+'/twiss.d-_back.dat')
    
            errdmax=twiss(savepath+'/twiss.d+.dat')
            errdmaxb=twiss(savepath+'/twiss.d+_back.dat')
    
            errcmax=twiss(savepath+'/twiss_c_max.dat')
            errcmin=twiss(savepath+'/twiss_c_min.dat')
    
            modelcor=twiss(savepath+'/twiss_'+namename+'_cor.dat')
    
            normal_pro=twiss(savepath+'/twiss_'+namename+'.dat')
            back_pro=twiss(savepath+'/twiss_'+namename+'_back_rev.dat')
    
    
            normal_pro.Cmatrix()
    
#             print normal_pro.C
#             sys.exit()
            # writing data in list
            phases=[filephasex,filephasey,filephasextot,filephaseytot]
            betah=[filedatax,errbetamin,errbetamax,errbetamin,errbetamax,erralfmin,erralfmax,erralfminb,erralfmaxb,filedataxA]
            betav=[filedatay,errbetamin,errbetamax,errbetamin,errbetamax,erralfmin,erralfmax,erralfminb,erralfmaxb,filedatayA]
            if disp==1:
                disph=[filedx,filendx,errdmin,errdmax,errdminb,errdmaxb]
                dispv=[filedy,errdmin,errdmax,errdminb,errdmaxb]
            else:
                disph=[]
                dispv=[]
            couple=[filecouple,filecoupleC,errcmax,errcmin]
            chromatic=[]
    
            # calling function to collect and write data
            print "Writing data for function ",namename
            getAndWriteData(namename,phases,betah,betav,disph,dispv,couple,chromatic,twisstwiss,modelcor,normal_pro,back_pro,savepath,elementswitch,accel)
    
    
    
            # gnuplot
            if elementswitch==0:
                
                startpos=twisstwiss.S[twisstwiss.indx[startbpm]]
                endpos=twisstwiss.S[twisstwiss.indx[endbpm]]
                print startpos,endpos
                run4plot(savepath,startpos,endpos,beta4plot,options.bb,path,namename,QXX,QYY,options.accel,method)
                
    return 0

# END main() ---------------------------------------------------------------------------------------

#===================================================================================================
# helper-functions
#===================================================================================================
def modelIntersect(expbpms, model):
    bpmsin=[]
    for bpm in expbpms:
            
        try:
            model.indx[bpm.replace("-","_").upper()]
            bpmsin.append(bpm)
        except:
            print bpm, "Not in Model"
    if len(bpmsin)==0:
        print "Zero intersection of Exp and Model"
        print "Please, provide a good Dictionary"
        print "Now we better leave!"
        sys.exit()            
    return bpmsin


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


def filterandfind(betaxx,betayy,element,segment,model,errorcut):
    '''
    Automatic BPM filter and BPM finder
    
    :Parameters:
        'betaxx': twiss
            twiss for horizontal beta function
        'betayy': twiss
            twiss for vertcal beta function
        'element': 
            instrument,collimator or IP . Specify as "null" if it as a segment
        'segment': 
            leftbpm and rightbpm where you want to compute the segment
        'model': twiss
            twiss model
        'errorcut': float
            errorcut for beta in percentage (err/bet)
    :Return: list
        Left and right BPM in a list
    :Exception:
        If n BPMS is  < 3 => System exit        
    '''
    # initial points
    translate={}
    location=[]

    if element!="null":
        seg=0
        try:
            sele=model.S[model.indx[element]]
        except:
            print element, " Not found in model=> System exit"
            sys.exit()
        location.append(sele)
        translate[sele]=[1,element]
        print "You selected an element"
    else:
        seg=1
        bpml=segment[0]
        bpmr=segment[1]
        try:
            sbpml=model.S[model.indx[bpml]]
            sbpmr=model.S[model.indx[bpmr]]
        except:
            print bpml,bpmr, " Not found in model=> System exit"
            sys.exit()            
        location.append(sbpml)
        location.append(sbpmr)
        try:
            betay=betayy.BETY[betayy.indx[bpml]]
            betax=betaxx.BETX[betaxx.indx[bpml]]
            betay=betayy.BETY[betayy.indx[bpmr]]
            betax=betaxx.BETX[betaxx.indx[bpmr]]
            translate[sbpml]=[1,bpml]
            translate[sbpmr]=[1,bpmr]
        except:
            translate[sbpml]=[0,bpml]
            translate[sbpmr]=[0,bpmr]
            
        print "You selected a segment"


    names=model.NAME
    goodones=0

    #filtering
    for name in names:

        if "BPM" in name:
            s=model.S[model.indx[name]]
            flag=0
            try:
                betay=betayy.BETY[betayy.indx[name]]
                betax=betaxx.BETX[betaxx.indx[name]]
                flag=1
            except:
                flag=0

                
            if flag==1:
                errbetax=betaxx.ERRBETX[betaxx.indx[name]]
                stdbetax=betaxx.STDBETX[betaxx.indx[name]]
                errbetay=betayy.ERRBETY[betayy.indx[name]]
                stdbetay=betayy.STDBETY[betayy.indx[name]]

                errx=sqrt(errbetax**2+stdbetax**2)
                erry=sqrt(errbetay**2+stdbetay**2)

                #print (betax>0),(betay>0),(errx<betax),(erry<betay),(errx>0),(erry>0),(((errx/betax)*100)<errorcut),(((erry/betay)*100)<errorcut)

                
                if (betax>0) and (betay>0) and (errx<betax) and (erry<betay) and (errx>0) and (erry>0) and (((errx/betax)*100)<errorcut) and(((erry/betay)*100)<errorcut):
                    translate[s]=[1,name]
                    goodones=goodones+1
                    location.append(s)

                elif seg==1 and (bpml in name or bpmr in name):
                    translate[s]=[0,name] # exclude giving bpms for segment
                    # if they dont pass the cut
            
                    
    
    location.sort()
    if goodones<3:
        print "Not enough BPMs! System exit"
        sys.exit()

    # finding the BPMs
    if seg==0:
        indin=location.index(sele)

        if indin==0:
            bpleft=translate[location[len(location)-2]][1]
            bpright=translate[location[1]][1]
        elif indin==(len(location)-1):
            bpleft=translate[location[indin-1]][1]
            bpright=translate[location[0]][1]
        else:
            bpleft=translate[location[indin-1]][1]
            bpright=translate[location[indin+1]][1]

    else:
        fll=translate[sbpml][0]
        flr=translate[sbpmr][0]
        print fll,flr

        if fll==1:
            bpleft=translate[sbpml][1]
        else:
            indin=location.index(sbpml)
            if indin==0:
                bpleft=translate[location[len(location)-2]][1]
            else:
                bpleft=translate[location[indin-1]][1]

        if flr==1:
            bpright=translate[sbpmr][1]
        else:
            indin=location.index(sbpmr)
            if indin==(len(location)-1):
                bpright=translate[location[0]][1]
            else:
                bpright=translate[location[indin+1]][1]

    print bpleft,bpright," Will be used for the propogation"

#    sys.exit()

    return [bpleft,bpright]

def TransverseDampers(twissp,twissb,element,model,savepath,phasex,phasey,errors):
    '''
    Function for getting the phase advance between the dampers and pick-ups
    
    :Parameters:
        'twissp': twiss?
            containing propagation from BPM to damper
        'twissb': twiss?
            containing back propagation from BPM to damper
        'element': string?
            element name
        'model': twiss
            twiss model
        'savepath': string
            where to save file
        'phasex': 
            measured phase advances
        'phasey': 
            measured phase advances
        'errors': 
            twissmin,twissminb,twissmax,twissmaxb
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    # initial
    dampermap={}
    bpmmap={}
    path=savepath+"phases.out"
    accel=model.SEQUENCE
    if os.path.exists(path):
        writer=open(path,"a")
        firsttime=0
    else:
        writer=open(path,"w")
        firsttime=1
        print >> writer,"@ BEAM %s",accel
        print >> writer,"@ NIM %s Not In Measurement"
        print >> writer,"* NAME1 NAME2 S1 S2 PHX PHXe PHXM PHY PHYe PHYM"
        print >> writer,"$ %s %s %le %le %le %le %le %le %le %le "
        
    
    #beam1
    dampermap['ADTKH.D5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.C5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.B5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.A5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']

    dampermap['ADTKV.A5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.B5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']    
    dampermap['ADTKV.C5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.D5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']

    bpmmap['LHCB1_1']=['BPMC.9L4.B1','BPMC.7L4.B1','BPMC.8L4.B1'] #H
    bpmmap['LHCB1_2']=['BPMCA.7R4.B1','BPMC.9R4.B1','BPMC.8R4.B1'] #V

    #beam2
    dampermap['ADTKV.D5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.C5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.B5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.A5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']

    dampermap['ADTKH.A5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.B5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.C5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.D5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']


    bpmmap['LHCB2_1']=['BPMCA.7R4.B2','BPMC.9R4.B2','BPMC.8R4.B2'] #H 
    bpmmap['LHCB2_2']=['BPMC.9L4.B2','BPMC.7L4.B2','BPMC.8L4.B2'] #V 

    ## main
    if firsttime==1:

        #
        counter=2
        #

        for count in range(0,counter):

            count=count+1

            print "counter for loop ", count

            bpms=bpmmap[accel+"_"+str(count)]
            refBpm=bpms[0]
            eightbpm=bpms[2]
            endbpm=bpms[1]

            # hor bpms's
            try: 
                in1x=phasex.indx[refBpm]
                
            except:
                in1x=-1

            if (in1x<len(phasex.indx)-1) and (in1x!=-1): # when BPM is not last
                if phasex.NAME[phasex.indx[refBpm]+1]==endbpm:
                    in2x=1
                    eightinphasex=0
                elif phasex.NAME[phasex.indx[refBpm]+2]==endbpm:
                    in2x=1
                    eightinphasex=1
                else: # not in measurement
                    in2x=-1
                    eightinphasex=0
            elif (in1x!=-1): # when BPM is last
                if phasex.NAME[0]==endbpm:
                    in2x=1
                    eightinphasex=0
                elif phasex.NAME[1]==endbpm:
                    in2x=1
                    eightinphasex=1
                else: # not in measurement
                    in2x=-1
                    eightinphasex=0
            else:
                in2x=-1    
                eightinphasex=0
                    

            # ver bpm's
            try:
                in1y=phasey.indx[refBpm]                
            except:
                in1y=-1

            if (in1y<len(phasey.indx)-1) and (in1y!=-1): # when BPM is not last
                if phasey.NAME[phasey.indx[refBpm]+1]==endbpm:
                    in2y=1
                    eightinphasey=0
                elif phasey.NAME[phasey.indx[refBpm]+2]==endbpm:
                    in2y=1
                    eightinphasey=1
                else: # not in measurement
                    in2y=-1
                    eightinphasey=0
            elif in1y!=-1: # when BPM is last
                if phasey.NAME[0]==endbpm:
                    in2y=1
                    eightinphasey=0
                elif phasey.NAME[1]==endbpm:
                    in2y=1
                    eightinphasey=1
                else: # not in measurement
                    in2y=-1
                    eightinphasey=0
            else:
                in2y=-1
                eightinphasey=0

            
        
            

            ###### H plane
            if (in1x!=-1) and (in2x!=-1):
                if eightinphasex==1:
                    phaseh="%.5f" %float(phasex.PHASEX[phasex.indx[refBpm]]+phasex.PHASEX[phasex.indx[eightbpm]])
                    errphaseh="%.5f" %float(phasex.STDPHX[phasex.indx[refBpm]]+phasex.STDPHX[phasex.indx[eightbpm]])
                    phasemodelh="%.5f" %float(phasex.PHXMDL[phasex.indx[refBpm]]+phasex.PHXMDL[phasex.indx[eightbpm]])
                else:
                    phaseh="%.5f" %float(phasex.PHASEX[phasex.indx[refBpm]])
                    errphaseh="%.5f" %float(phasex.STDPHX[phasex.indx[refBpm]])
                    phasemodelh="%.5f" %float(phasex.PHXMDL[phasex.indx[refBpm]])
            else:
                print "Horizontal plane not found for transverse dampers pick-ups"
                phaseh='NIM'
                errphaseh='NIM'
                phasemodelh='NIM'                        
                    
            ###### V plane
            print in1y,in2y,eightinphasey
            if (in1y!=-1) and (in2y!=-1):
                if eightinphasey==1:
                    phasev="%.5f" %float(phasey.PHASEY[phasey.indx[refBpm]]+phasey.PHASEY[phasey.indx[eightbpm]])
                    errphasev="%.5f" %float(phasey.STDPHY[phasey.indx[refBpm]]+phasey.STDPHY[phasey.indx[eightbpm]])
                    phasemodelv="%.5f" %float(phasey.PHYMDL[phasey.indx[refBpm]]+phasey.PHYMDL[phasey.indx[eightbpm]])
                else:
                    
                    phasev="%.5f" %float(phasey.PHASEY[phasey.indx[refBpm]])
                    errphasev="%.5f" %float(phasey.STDPHY[phasey.indx[refBpm]])
                    phasemodelv="%.5f" %float(phasey.PHYMDL[phasey.indx[refBpm]])
                                    
            else:
                print "Vertical plane not found for transverse dampers pick-ups"
                phasev='NIM'
                errphasev='NIM'
                phasemodelv='NIM'

                
            print >> writer,refBpm,endbpm,model.S[model.indx[refBpm]],model.S[model.indx[endbpm]],phaseh,errphaseh,phasemodelh,phasev,errphasev,phasemodelv
                                          
    if element in dampermap:
        bpms=dampermap[element]
        passs=1
    else:
        print "WARN: Element ",element," Not found in dampermap"
        passs=0

    if passs==1:
        muxbpm1=twissp.MUX[twissp.indx[bpms[0]]]
        muxbpm2=twissp.MUX[twissp.indx[bpms[1]]]
        mux=twissp.MUX[twissp.indx[element]]
        muybpm1=twissp.MUY[twissp.indx[bpms[0]]]
        muybpm2=twissp.MUY[twissp.indx[bpms[1]]]
        muy=twissp.MUY[twissp.indx[element]]

        muxbpm1m=model.MUX[model.indx[bpms[0]]]
        muxbpm2m=model.MUX[model.indx[bpms[1]]]
        muxm=model.MUX[model.indx[element]]

        muybpm1m=model.MUY[model.indx[bpms[0]]]
        muybpm2m=model.MUY[model.indx[bpms[1]]]
        muym=model.MUY[model.indx[element]]

        pha1xp=abs(mux-muxbpm1)
        pha2xp=abs(muxbpm2-mux)
        pha1yp=abs(muy-muybpm1)
        pha2yp=abs(muybpm2-muy)
        
        muxbpm1=twissb.MUX[twissb.indx[bpms[0]]]
        muxbpm2=twissb.MUX[twissb.indx[bpms[1]]]
        mux=twissb.MUX[twissb.indx[element]]
        muybpm1=twissb.MUY[twissb.indx[bpms[0]]]
        muybpm2=twissb.MUY[twissb.indx[bpms[1]]]
        muy=twissb.MUY[twissb.indx[element]]
        

        pha1xb=abs(muxbpm1-mux)
        pha2xb=abs(mux-muxbpm2)
        pha1yb=abs(muybpm1-muy)
        pha2yb=abs(muy-muybpm2)


        pha1xm=abs(muxm-muxbpm1m)
        pha2xm=abs(muxbpm2m-muxm)
        pha1ym=abs(muym-muybpm1m)
        pha2ym=abs(muybpm2m-muym)




        #mean
        pha1x=(pha1xp+pha1xb)/2
        pha1y=(pha1yp+pha1yb)/2
        pha2x=(pha2xp+pha2xb)/2
        pha2y=(pha2yp+pha2yb)/2

        #error
        twissmin=errors[0]
        twissminb=errors[1]
        twissmax=errors[2]
        twissmaxb=errors[3]
        
        minn=abs(abs(twissmin.MUX[twissmin.indx[element]]-twissmin.MUX[twissmin.indx[bpms[0]]])-abs(twissminb.MUX[twissminb.indx[element]]-twissminb.MUX[twissminb.indx[bpms[0]]]))
        maxx=abs(abs(twissmax.MUX[twissmax.indx[element]]-twissmax.MUX[twissmax.indx[bpms[0]]])-abs(twissmaxb.MUX[twissmaxb.indx[element]]-twissmaxb.MUX[twissmaxb.indx[bpms[0]]]))
        pha1xe=minn+maxx


        minn=abs(abs(twissmin.MUX[twissmin.indx[bpms[1]]]-twissmin.MUX[twissmin.indx[element]])-abs(twissminb.MUX[twissminb.indx[bpms[1]]]-twissminb.MUX[twissminb.indx[element]]))
        maxx=abs(abs(twissmax.MUX[twissmax.indx[bpms[1]]]-twissmax.MUX[twissmax.indx[element]])-abs(twissmaxb.MUX[twissmaxb.indx[bpms[1]]]-twissmaxb.MUX[twissmaxb.indx[element]]))
        pha2xe=minn+maxx

        minn=abs(abs(twissmin.MUY[twissmin.indx[element]]-twissmin.MUY[twissmin.indx[bpms[0]]])-abs(twissminb.MUY[twissminb.indx[element]]-twissminb.MUY[twissminb.indx[bpms[0]]]))
        maxx=abs(abs(twissmax.MUY[twissmax.indx[element]]-twissmax.MUY[twissmax.indx[bpms[0]]])-abs(twissmaxb.MUY[twissmaxb.indx[element]]-twissmaxb.MUY[twissmaxb.indx[bpms[0]]]))
        pha1ye=minn+maxx

        minn=abs(abs(twissmin.MUY[twissmin.indx[bpms[1]]]-twissmin.MUY[twissmin.indx[element]])-abs(twissminb.MUY[twissminb.indx[bpms[1]]]-twissminb.MUY[twissminb.indx[element]]))
        maxx=abs(abs(twissmax.MUY[twissmax.indx[bpms[1]]]-twissmax.MUY[twissmax.indx[element]])-abs(twissmaxb.MUY[twissmaxb.indx[bpms[1]]]-twissmaxb.MUY[twissmaxb.indx[element]]))
        pha2ye=minn+maxx
        
        print >> writer,bpms[0],element,model.S[model.indx[bpms[0]]],model.S[model.indx[element]],"%.5f" %float(pha1x),"%.5f" %float(pha1xe),"%.5f" %float(pha1xm),"%.5f" %float(pha1y),"%.5f" %float(pha1ye),"%.5f" %float(pha1ym)
        print >> writer,element,bpms[1],model.S[model.indx[element]],model.S[model.indx[bpms[1]]],"%.5f" %float(pha2x),"%.5f" %float(pha2xe),"%.5f" %float(pha2xm),"%.5f" %float(pha2y),"%.5f" %float(pha2ye),"%.5f" %float(pha2ym)

    writer.close()

def getIP(betameA,basetwiss,betatwiss,alfatwiss,model,phasex,phasey,name,accel,path):
    '''
    Function calculating the optics parameters at the IP
    
    :Parameters:
        'betameA': list
            contaning horizontal and vertical measured amplitude functions
        'basetwiss': list
            list with propogated and back propogated twiss
        'betatwiss': list
            list of twisses for beta errors
        'alfatwiss': list
            list of twisses for alfa errors
        'model': twiss
            twiss model
        'phasex': 
            measured phase advances
        'phasey': 
            measured phase advances
        'name': string
            name of the IP
        'accel': string
            name of the accelerator
        'path': string
            where to save file
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    if os.path.isfile(path+"/IP_para.out"):
        ipfilepara=open(path+"/IP_para.out","a")
    else:
        ipfilepara=open(path+"/IP_para.out","w")
        print >> ipfilepara,"* NAME S BETX EBETX X EX BETY EBETY Y EY BETX_ph EBETX_ph BETY_ph EBETY_ph"
        print >> ipfilepara,"$ %s %le %le %le %le %le %le %le %le %le %le %le %le %le"


    if os.path.isfile(path+"/IP_pro_x.out"):
        ipfileprox=open(path+"/IP_pro_x.out","a")
    else:
        ipfileprox=open(path+"/IP_pro_x.out","w")
        print >> ipfileprox,"* NAME S BETSTARX EBETSTARX X[cm] EX[cm] BETIPX EBETIPX ALFIPX EALFIPX"
        print >> ipfileprox,"$ %s %le %le %le %le %le %le %le %le %le"

    if os.path.isfile(path+"/IP_pro_y.out"):
        ipfileproy=open(path+"/IP_pro_y.out","a")
    else:
        ipfileproy=open(path+"/IP_pro_y.out","w")
        print >> ipfileproy,"* NAME S BETY EBETY Y[cm] EY[cm] BETIPY EBETIPY ALFIPY EALFIPY"
        print >> ipfileproy,"$ %s %le %le %le %le %le %le %le %le %le"

    ## constructor
    if "B1" in accel:
        accelb="B1"
    else:
        accelb="B2"
    twissp=basetwiss[0]
    twissb=basetwiss[1]

    twissbmip=betatwiss[0]
    twissbmib=betatwiss[1]
    twissbmap=betatwiss[2]
    twissbmab=betatwiss[3]

    twissamip=alfatwiss[0]
    twissamib=alfatwiss[1]
    twissamap=alfatwiss[2]
    twissamab=alfatwiss[3]

    ampbetx=betameA[0]
    ampbety=betameA[1]

    ## beta from parabola calculation
    # Map bpm's
    bpmmap={}
    bpmmap['IP1']=["BPMSW.1L1.","BPMSW.1R1."]
    bpmmap['IP2']=["BPMSW.1L2.","BPMSW.1R2."]
    bpmmap['IP5']=["BPMSW.1L5.","BPMSW.1R5."]
    bpmmap['IP8']=["BPMSW.1L8.","BPMSW.1R8."]
    
    if(name in bpmmap):
        # Removed unusedphix, ephix, phiy, ephiy. 
        # If needed --> git commit eca47b2d06ab7eab911c830f30dc222a8bcd91a0
        # --vimaier
        betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_p,ebetastar_h_p,betastar_v_p,ebetastar_v_p=getIPfrompara(bpmmap[name][0]+accelb,bpmmap[name][1]+accelb,ampbetx,ampbety,phasex,phasey)
    else:
        betastar_h_p=ebetastar_h_p=betastar_v_p=ebetastar_v_p=betastar_h=errbx=location_h=errsx=betastar_v=errby=location_v=errsy=0    

    print >> ipfilepara,name,model.S[model.indx[name]],betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_p,ebetastar_h_p,betastar_v_p,ebetastar_v_p
    ipfilepara.close()

    ## beta from propogations
    #gathering info
    #error alfa

    dalfaxp=1/abs(twissamap.ALFX[twissamap.indx[name]]-twissamip.ALFX[twissamip.indx[name]])
    dalfaxb=1/abs(twissamab.ALFX[twissamab.indx[name]]-twissamib.ALFX[twissamib.indx[name]])

    dalfayp=1/abs(twissamap.ALFY[twissamap.indx[name]]-twissamip.ALFY[twissamip.indx[name]])
    dalfayb=1/abs(twissamab.ALFY[twissamab.indx[name]]-twissamib.ALFY[twissamib.indx[name]])

    #error beta's
    dbxp=1/abs(twissbmap.BETX[twissbmap.indx[name]]-twissbmip.BETX[twissbmip.indx[name]])
    dbxb=1/abs(twissbmab.BETX[twissbmab.indx[name]]-twissbmib.BETX[twissbmib.indx[name]])

    dbyp=1/abs(twissbmap.BETY[twissbmap.indx[name]]-twissbmip.BETY[twissbmip.indx[name]])
    dbyb=1/abs(twissbmab.BETY[twissbmab.indx[name]]-twissbmib.BETY[twissbmib.indx[name]])

    bb=twissb.ALFX[twissb.indx[name]]
    if bb<0.0:
        Balfx=abs(bb)
    else:
        Balfx=-bb

    bb=twissb.ALFY[twissb.indx[name]]
    if bb<0.0:
        Balfy=abs(bb)
    else:
        Balfy=-bb
    
    alfax=(1/(dalfaxp+dalfaxb))*(dalfaxp*twissp.ALFX[twissp.indx[name]]+Balfx*dalfaxb)

    ealfax=(sqrt((1/dalfaxp)**2+(1/dalfaxb)**2))/2
    
    alfay=(1/(dalfayp+dalfayb))*(dalfayp*twissp.ALFY[twissp.indx[name]]+Balfy*dalfayb)
    ealfay=(sqrt((1/dalfayp)**2+(1/dalfayb)**2))/2

    betxip=(1/(dbxp+dbxb))*(dbxp*twissp.BETX[twissp.indx[name]]+dbxb*twissb.BETX[twissb.indx[name]])
    
    betyip=(1/(dbyp+dbyb))*(dbyp*twissp.BETY[twissp.indx[name]]+dbyb*twissb.BETY[twissb.indx[name]])
    errbetxip=sqrt((1/dbxp)**2+(1/dbxb)**2)/2
    errbetyip=sqrt((1/dbyp)**2+(1/dbyb)**2)/2

    #betstar and waist
    #x
    betastar,ebetastar,waist,ewaist=getIPfromProp(betxip,errbetxip,alfax,ealfax)
    print >> ipfileprox,name,model.S[model.indx[name]],betastar,ebetastar,waist,ewaist,round(betxip,3),round(errbetxip,3),round(alfax,4),round(ealfax,4)
    ipfileprox.close()

    #y
    betastar,ebetastar,waist,ewaist=getIPfromProp(betyip,errbetyip,alfay,ealfay)    
    print >> ipfileproy,name,model.S[model.indx[name]],betastar,ebetastar,waist,ewaist,round(betyip,3),round(errbetyip,3),round(alfay,4),round(ealfay,4)
    ipfileproy.close()

    ##dispersion

def getIPfromProp(betaip,errbetaip,alfaip,ealfaip):

    #values
    betastar=betaip/(1+alfaip**2)
    waist=alfaip*betaip #(sign flip)

    #errors
    ewaist=((ealfaip/abs(alfaip))+(errbetaip/abs(betaip)))*abs(waist)
    ebetastar=sqrt(errbetaip**2+(((-2*alfaip)/(1+alfaip**2)**2)*alfaip)**2 )

    waist=waist*100 # transferring to CM!!!!
    ewaist=ewaist*100 # transferring to CM!!!!


    return round(betastar,3),round(ebetastar,3),round(waist,3),round(ewaist,3)


def getIPfrompara(bpmleft,bpmright,betax,betay,phasex,phasey):
    '''
    Function calculating beta at IP (based on Rioichy thesis)
    b(s)=b*+(s^2/b*)
    
    :Parameters:
        'bpmleft,bpmright': <dtype?>
            bpm's used to caculate the beta at the IP
        'ampbetx,ampbety': twiss
            twiss containing the horizontal and vertical beta from amp
        'phasex,phasey': <dtype?>
            measured phases
    :Return: <dtype?>
        beta at waist and offset
    '''
    #left horizontal
    try:
        blx=betax.BETX[betax.indx[bpmleft]]
        betax.BETXSTD[betax.indx[bpmleft]]
        hlp=1
    except:
        hlp=0
            
    #right horizontal
    try:
        brx=betax.BETX[betax.indx[bpmright]]
        betax.BETXSTD[betax.indx[bpmright]]
        hrp=1
    except:
        hrp=0
        

    #left vertical
    try:
        bly=betay.BETY[betay.indx[bpmleft]]
        betay.BETYSTD[betay.indx[bpmleft]]
        vlp=1
    except:
        vlp=0

    #right vertical
    try:
        bry=betay.BETY[betay.indx[bpmright]]
        betay.BETYSTD[betay.indx[bpmright]]
        vrp=1
    except:
        vrp=0

    #getting phaseadvance
    #inx1=phasex.indx[bpmleft]
    #inx2=phasex.indx[bpmright]
    try:
        inx1=phasex.indx[bpmleft]
        inx2=phasex.indx[bpmright]
        if (inx1==(inx2-1)) or (inx2==0 and inx1==len(phasex.NAME)-1):
            pxp=1
        else:
            pxp=0
    except:
        pxp=0

    try:
        iny1=phasey.indx[bpmleft]
        iny2=phasey.indx[bpmright]
        if iny1==(iny2-1) or (inx2==0 and inx1==len(phasex.NAME)-1):
            pyp=1
        else:
            pyp=0
    except:
        pyp=0


    #horizontal
    if(hrp==1) and (hlp==1) and (pxp==1):

        sxl=phasex.S[inx1]
        sxr=phasex.S[inx2]
        if (sxl==phasex.S[len(phasex.NAME)-1]): #if start point is in center of IP ! 
            L=sxr
        else:
            L=(sxr-sxl)/2
        
        phix=phasex.PHASEX[inx1]
        
        betastar_h=(2*sqrt(blx)*sqrt(brx)*sin(phix*2*pi))/(blx+brx-2*sqrt(blx)*sqrt(brx)*cos(2*pi*phix))*L
        location_h=((blx-brx)/(blx+brx-2*sqrt(blx)*sqrt(brx)*cos(2*pi*phix)))*L

        

        errbx=0
        errsx=0

        
        betastar_h_phase=(L)/tan(phix*pi)
        betastar_h_phase_e=(1+cos(2*phix*pi))*L/2
        
        
    else:
        print "WARN: Will not calculate horizontal IP"
        phix="NIM"
        betastar_h="NIM"
        location_h="NIM"
        errbx="NIM"
        errsx="NIM"
        betastar_h_phase="NIM"
        betastar_h_phase_e="NIM"
    #vertical
    if (vrp==1) and (vlp==1)and (pyp==1):

        syl=phasey.S[iny1]
        syr=phasey.S[iny2]
        if (syl==phasey.S[len(phasey.NAME)-1]): #if start point is in center of IP ! 
            L=syr
        else:
            L=(syr-syl)/2
        

        phiy=phasey.PHASEY[iny1]
        
        betastar_v=(2*sqrt(bly)*sqrt(bry)*sin(phiy*2*pi))/(bly+bry-2*sqrt(bly)*sqrt(bry)*cos(2*pi*phiy))*L
        location_v=((bly-bry)/(bly+bry-2*sqrt(bly)*sqrt(bry)*cos(2*pi*phiy)))*L

        betastar_v_phase=(L)/tan(phiy*pi)
        betastar_v_phase_e=(1+cos(2*phiy*pi))*L/2

        errby=0
        errsy=0
    else:
        print "WARN: Will not calculate vertical IP"
        phiy="NIM"
        betastar_v="NIM"
        location_v="NIM"
        errby="NIM"
        errsy="NIM"
        betastar_v_phase="NIM"
        betastar_v_phase_e="NIM"        
        
    return betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_phase,betastar_h_phase_e,betastar_v_phase,betastar_v_phase_e

def getAndWriteData(namename,phases,betah,betav,disph,dispv,couple,chromatic,model,modelcor,modelp,modelb,path,switch,accel):
    '''
    Function that returns the optics function at the given element
    
    :Parameters:
        'namename': string
            name of the element
        'phase': list
            horizontal and vertical phase [phasex,phasey,phasext,phaseyt]
        'betah': list
            horizontal beta [betax,betaxminp,betaxmaxp,betaxminb,betaxmaxb,alfaxminp,alfaxmaxp,alfaxminb,alfaxmaxb,ampbetax]
        'betav': list
            vertical beta [betay,betayminp,betaymaxp,betayminb,betaymaxb,alfayminp,alfaymaxp,alfayminb,alfaymaxb,ampbetay]
        'disph': list
            horizontal dispersion [dispx,ndispx,dispminxp,dispmaxp,dispminxb,dispmaxb]
        'dispv': list
            vertical dispersion [dispy,dispminyp,dispmayp,dispminyb,dispmayb]
        'couple': list 
            containing the coupling [coupleme,couplecme,couplemip,couplemap]
        'chromatic': list 
            containing the chromatic ouput [wx,wxminp,wxmqxp,wy,wyminp,wymaxp]
        'model': twiss
            model in twiss format
        'modelcor': twiss
            model containing twiss with corrections values (former known as model_play)
        'modelp': twiss
            model from propagation
        'modelb': twiss
            model from back propagation
        'path': string
            where to save file
        'switch': int
            to tell if it is either element(1)/segment(0)
        'accel': string
            name of the accelerator
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    ###################################################################
    # Function that returns the optics function at the given element
    # Parameters:
    # - name : name of the element
    # - phase: (list) horizontal and vertical phase [phasex,phasey,phasext,phaseyt]
    # - betah : (list) horizontal beta [betax,betaxminp,betaxmaxp,betaxminb,betaxmaxb,alfaxminp,alfaxmaxp,alfaxminb,alfaxmaxb,ampbetax]
    # - betav : (list) vertical beta [betay,betayminp,betaymaxp,betayminb,betaymaxb,alfayminp,alfaymaxp,alfayminb,alfaymaxb,ampbetay]
    # - disph : (list) horizontal dispersion [dispx,ndispx,dispminxp,dispmaxp,dispminxb,dispmaxb]
    # - dispv : (list) vertical dispersion [dispy,dispminyp,dispmayp,dispminyb,dispmayb]
    # - couple : twiss containing the coupling [coupleme,couplecme,couplemip,couplemap]
    # - chromatic : containing the chromatic ouput [wx,wxminp,wxmqxp,wy,wyminp,wymaxp]
    # - model : model in twiss format
    # - modelcor : model containing twiss with corrections values (former known as model_play)
    # - modelp : model from propagation
    # - modelb : model from back propagation
    # - path : path where to save
    # - switch : to tell if it is either element(1)/segment(0)
    ###################################################################

    print "INFO: Start writing files",switch


    #### gathering file
    if switch==1:
        if os.path.isfile(path+"sbs_summary_bet.out"):
            print "INFO: Updating summary file"
            filesum_b=open(path+"sbs_summary_bet.out","a")
            filesum_c=open(path+"sbs_summary_cou.out","a")
            filesum_d=open(path+"sbs_summary_disp.out","a")

        else:
            print "INFO: Creating summary file"
            filesum_b=open(path+"sbs_summary_bet.out","w")
            filesum_c=open(path+"sbs_summary_cou.out","w")
            filesum_d=open(path+"sbs_summary_disp.out","w")


            # beta
            print >> filesum_b, "* NAME S BETXP ERRBETXP BETXMDL ALFXP ERRALFXP ALFXMDL BETY ERRBETY BETYMDL ALFA ERRALFY ALFYMDL MDL_S"
            print >> filesum_b, "$ %s %le %le %le %le %le %le %le %le %le %le %le %le %le %le"

            #coupling
            print >> filesum_c,"* NAME   S   f1001 f1001re  f1001im    f1010   f1010re   f1010im  f1001_PLAY ef1001_play   f1001re_PLAY  f1001im_PLAY    f1010_PLAY ef1010_play   f1010re_PLAY   f1010im_PLAY C11Mo C12Mo C21Mo C22Mo ANDMo C11_cor eC11_cor C12_cor eC12_cor C21_cor eC21_cor C22_cor eC22_cor ANG_cor eANG_cor S_MODEL"
            print >> filesum_c,"$ %s %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le"

            # disp and chromatic
            if len(chromatic)!=0:
        
                print >> filesum_d, "* NAME S DX_MDL DX_PLAY EDX_PLAY DPX_MDL DPX_PLAY EDPX_PLAY DY_MDL DY_PLAY EDY_PLAY DPY_MDL DPY_PLAY EDPY_PLAY WX_MDL WX_PLAY eWX_play PHIX_MDL PHIX_PLAY ePHIX_PLAY WY_MDL WY_PLAY eWY_play PHIY_MDL PHIY_PLAY ePHIY_PLAY S_MODEL"
                print >> filesum_d, "$ %s %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le"

            else:
                print >> filesum_d, "* NAME S DX_MDL DX_PLAY EDX_PLAY DPX_MDL DPX_PLAY EDPX_PLAY  DY_MDL DY_PLAY EDY_PLAY DPY_MDL DPY_PLAY EDPY_PLAY  S_MODEL"
                print >> filesum_d, "$ %s %le %le %le %le %le %le %le %le %le  %le %le %le %le %le"

    phasex=phases[0]
    phasey=phases[1]        

    ##beta
    #-horizontal
    filex=open(path+"sbsbetax_"+namename+".out","w")
    filexa=open(path+"sbsalfax_"+namename+".out","w")

    if switch==0:    
        print >> filex,"* NAME S BETX ERRBETX BETXAMP ERRBETXAMP BETXP ERRBETXP BETXMDL MODEL_S"
        print >> filex,"$ %s %le %le %le %le %le  %le %le %le %le"
    
        print >> filexa,"* NAME S ALFX ERRALFX ALFXP ERRALFXP MODEL_S"
        print >> filexa,"$ %s %le %le %le %le %le %le"
    else:
        print >> filex,"* NAME S BETXP ERRBETXP BETXMDL MODEL_S"
        print >> filex,"$ %s %le %le %le %le %le"
    
        print >> filexa,"* NAME S ALFXP ERRALFXP ALFXMDL MODEL_S"
        print >> filexa,"$ %s %le %le %le %le %le"        
    
    bme=betah[0] # measurement
    bmea=betah[9] # measurement amp
    bmip=betah[1] # min p
    bmap=betah[2] # max p
    bmib=betah[3] # min b
    bmab=betah[4] # max b
    amip=betah[5] # min p
    amap=betah[6] # max p
    amib=betah[7] # min b
    amab=betah[8] # max b

    
    if switch==0: # only enter when it is a segment
        bpms=intersect([bme,bmip,model,modelcor,modelp,modelb,bmea])
    else:
        bpms=intersect([bmip,model,modelcor,modelp,modelb])

    for bpm in bpms:
        s=bpm[0]
        name=bpm[1]
        smo=model.S[model.indx[name]]

        #beta model
        betam=model.BETX[model.indx[name]]

        #beta pro
        betape=abs(bmap.BETX[bmap.indx[name]]-bmip.BETX[bmip.indx[name]])
        betabe=abs(bmab.BETX[bmab.indx[name]]-bmib.BETX[bmib.indx[name]])
        #ebep=sqrt((1/betape)**2+(1/betabe)**2)/2
        ebep=sqrt(((betape)**2+(betabe)**2)/2)
                
        #alfa
        amo=model.ALFX[model.indx[name]]
                
        alfape=abs(amap.ALFX[amap.indx[name]]-amip.ALFX[amip.indx[name]])
        alfabe=abs(amab.ALFX[amab.indx[name]]-amib.ALFX[amib.indx[name]])
        eaep=sqrt(((alfape)**2+(alfabe)**2)/2)

        betap=modelp.BETX[modelp.indx[name]]
        betab=modelb.BETX[modelb.indx[name]]
        bep=(1/(betape+betabe))*(betape*betap+betabe*betab)
                
        alfap=modelp.ALFX[modelp.indx[name]]                
        alfab=modelb.ALFX[modelb.indx[name]]
        aep=(1/(alfape+alfabe))*(alfape*alfap+alfabe*alfab)

                
        if switch==0:
            
                alfame=bme.ALFX[bme.indx[name]]
                ealfame=sqrt((bme.ERRALFX[bme.indx[name]]**2+bme.STDALFX[bme.indx[name]]**2)/2)
                #beta me
                betaa=bmea.BETX[bmea.indx[name]]
                ebetaa=bmea.BETXSTD[bmea.indx[name]]
                betame=bme.BETX[bme.indx[name]]
                ebetame=sqrt((bme.ERRBETX[bme.indx[name]]**2+bme.STDBETX[bme.indx[name]]**2)/2)
                bep=modelcor.BETX[modelcor.indx[name]]
                alfap=modelp.ALFX[modelp.indx[name]]

                
                print >> filexa,name,s,alfame,ealfame,aep,eaep,amo,smo
                print >> filex,name,s,betame,ebetame,betaa,ebetaa,bep,ebep,betam,smo
        else:
                #print "Adding", name," to the summary ",namename

                print >> filexa,name,s,aep,eaep,amo,smo
                print >> filex,name,s,bep,ebep,betam,smo

                if namename in name:
                    fileb1=name+" "+str(s)+" "+str(round(bep,2))+" "+str(round(ebep,2))+" "+str(round(betam,2))+" "+str(round(aep,4))+" "+str(round(eaep,4))+" "+str(round(amo,4))
        
    filex.close()
    filexa.close()

    #-vertical
    filey=open(path+"sbsbetay_"+namename+".out","w")
    fileya=open(path+"sbsalfay_"+namename+".out","w")
                
    if switch==0:
                print >> filey,"* NAME S BETY ERRBETY  BETYAMP ERRBETYAMP BETYP ERRBETYP BETYMDL MDL_S"
                print >> filey,"$ %s %le %le %le %le %le %le %le %le %le "
                print >> fileya,"* NAME S ALFY ERRALFY ALFYP ERRALFYP ALFMDL MODEL_S"
                print >> fileya,"$ %s %le %le %le %le %le %le"
    else:
                print >> filey,"* NAME S BETY ERRBETY BETYMDL MDL_S"
                print >> filey,"$ %s %le %le %le %le %le"        
                print >> fileya,"* NAME S ALFA ERRALFY ALFYMDL MDL_S"
                print >> fileya,"$ %s %le %le %le %le %le"

    
    bme=betav[0] # measurement
    bmea=betav[9] # measurement amp            
    bmip=betav[1] # min p
    bmap=betav[2] # max p
    bmib=betav[3] # min b
    bmab=betav[4] # max b
    amip=betav[5] # min p
    amap=betav[6] # max p
    amib=betav[7] # min b
    amab=betav[8] # max b
                
    if switch==0: # only enter when it is a segment
                bpms=intersect([bme,bmip,model,modelcor,modelp,modelb,bmea])
    else:
                bpms=intersect([bmip,model,modelcor,modelp,modelb])        

    for bpm in bpms:
        s=bpm[0]
        name=bpm[1]
        smo=model.S[model.indx[name]]

        #beta model
        betam=model.BETY[model.indx[name]]

                      
        #beta

        betape=abs(bmap.BETY[bmap.indx[name]]-bmip.BETY[bmip.indx[name]])
        betabe=abs(bmab.BETY[bmab.indx[name]]-bmib.BETY[bmib.indx[name]])
        #ebep=sqrt((1/betape)**2+(1/betabe)**2)/2
        ebep=sqrt((betape)**2+(betabe)**2)/2
 
        #alfa
        amo=model.ALFY[model.indx[name]]
        
        alfape=abs(amap.ALFY[amap.indx[name]]-amip.ALFY[amip.indx[name]])
        alfabe=abs(amab.ALFY[amab.indx[name]]-amib.ALFY[amib.indx[name]])
        eaep=sqrt((alfape)**2+(alfabe)**2)/2

        

        
        if switch==0:
                #beta me               
            betaa=bmea.BETY[bmea.indx[name]]
            ebetaa=bmea.BETYSTD[bmea.indx[name]]
            betame=bme.BETY[bme.indx[name]]                             
            ebetame=sqrt((bme.ERRBETY[bme.indx[name]]**2+bme.STDBETY[bme.indx[name]]**2)/2)
            print >> filey,name,s,betame,ebetame,betaa,ebetaa,bep,ebep,betam,smo            
            #alfa me
            bep=modelcor.BETY[modelcor.indx[name]]
            aep=modelcor.ALFY[modelcor.indx[name]]
               
            alfame=bme.ALFY[bme.indx[name]]
            ealfame=sqrt((bme.ERRALFY[bme.indx[name]]**2+bme.STDALFY[bme.indx[name]]**2)/2)
                
            print >> fileya,name,s,alfame,ealfame,aep,eaep,amo,smo
        else:
            betap=modelp.BETY[modelp.indx[name]]
            betab=modelb.BETY[modelb.indx[name]]
            bep=(1/(betape+betabe))*(betape*betap+betabe*betab)            
            alfap=modelp.ALFY[modelp.indx[name]]
            alfab=modelb.ALFY[modelb.indx[name]]
            aep=(1/(alfape+alfabe))*(alfape*alfap+alfabe*alfab)
            print >> fileya,name,s,aep,eaep,amo,smo
            print >> filey,name,s,bep,ebep,betam,smo

            if namename in name:
                print >>filesum_b,fileb1,round(bep,2),round(ebep,2),round(betam,2),round(aep,4),round(eaep,4),round(amo,4),round(smo,2)
            
        
    filey.close()
    fileya.close()
    

    ##dispersion
    if len(disph)!=0:
        dme=disph[0]
        dnme=disph[1]
        dminp=disph[2]
        dmaxp=disph[3]
        dminb=disph[4]
        dmaxb=disph[5]

            #dispersion
        filex= open(path+'/sbsDx_'+namename+'.out','w')


        if switch==0:
            print >>filex,"* NAME  S  DX  STDDX  DX_MDL DX_PLAY EDX_PLAY DPX DPX_MDL DPX_PLAY EDPX_PLAY MODEL_S"
            print >>filex,"$ %s %le %le %le %le %le %le %le %le %le %le %le"

            bpms=intersect([dme,dminp,modelp,model,modelb])
        else:
            print >>filex,"* NAME  S  DX_MDL DX_PLAY EDX_PLAY DPX_MDL DPX_PLAY EDPX_PLAY  MODEL_S"
            print >>filex,"$ %s %le %le %le %le %le %le %le %le"
            bpms=intersect([dminp,modelp,model,modelb])        
    


        for bpm in bpms:
            s=bpm[0]
            name=bpm[1]

            #model
            dmo=model.DX[model.indx[name]]
            smo=model.S[model.indx[name]]
            dpmo=model.DPX[model.indx[name]]
            # pro
            dpe=abs(dmaxp.DX[dmaxp.indx[name]]-dminp.DX[dminp.indx[name]])
            dbe=abs(dmaxb.DX[dmaxb.indx[name]]-dminb.DX[dminb.indx[name]])
            dppe=abs(dmaxp.DPX[dmaxp.indx[name]]-dminp.DPX[dminp.indx[name]])
            dpbe=abs(dmaxb.DPX[dmaxb.indx[name]]-dminb.DPX[dminb.indx[name]])               
            edep=sqrt((dpe)**2+(dbe)**2)/2
            edepp=sqrt((dppe)**2+(dpbe)**2)/2
            

            if switch==0:
                #mea
                dmea=dme.DX[dme.indx[name]]
                edmea=dme.STDDX[dme.indx[name]]
                dpmea=dme.DPX[dme.indx[name]]
            
                #pro
                dep=modelcor.DX[modelcor.indx[name]]
                depp=modelcor.DPX[modelcor.indx[name]]

                print >>filex,name,s,dmea,edmea,dmo,dep,edep,dpmea,dpmo,depp,edepp,smo
            else:
                dp=modelp.DX[modelp.indx[name]]
                db=modelb.DX[modelb.indx[name]]
                dep=(1/(dpe+dbe))*(dpe*dp+dbe*db)
                
                dpp=modelp.DPX[modelp.indx[name]]
                dpb=modelb.DPX[modelb.indx[name]]
                try:
                    depp=(1/(dppe+dpbe))*(dppe*dpp+dpbe*dpb)
                except:
                    depp=0
                
                print >>filex,name,s,dmo,dep,edep,dpmo,depp,edepp,smo

                if namename in name:
                    filed1=name+" "+str(s)+" "+str(dmo)+" "+str(dep)+" "+str(edep)+" "+str(dpmo)+" "+str(depp)+" "+str(edepp)


        filex.close()
        
            #normalized dispersion
        filex= open(path+'/sbsNDx_'+namename+'.out','w')


        if switch==0:
            print >> filex,"* NAME  S  NDX   STDNDX NDX_MDL NDX_PLAY ENDX_PLAY MODEL_S "
            print >> filex,"$ %s   %le  %le   %le   %le    %le %le    %le "
            bpms=intersect([dnme,dminp,modelp,model,modelb,bme])
        else:
            print >> filex,"* NAME  S  NDX_MDL NDX_PLAY ENDX_PLAY MODEL_S "
            print >> filex,"$ %s   %le   %le    %le %le    %le "
            bpms=intersect([dminp,modelp,model,modelb,modelcor])

        for bpm in bpms:
            s=bpm[0]
            name=bpm[1]

            #model
            dmo=model.DX[model.indx[name]]
            smo=model.S[model.indx[name]]

                # pro
            dpe=(abs(dmaxp.DX[dmaxp.indx[name]]-dminp.DX[dminp.indx[name]]))/sqrt(modelp.BETX[modelp.indx[name]])
            dbe=(abs(dmaxb.DX[dmaxb.indx[name]]-dminb.DX[dminb.indx[name]]))/sqrt(modelp.BETX[modelp.indx[name]])
            try:
                edep=sqrt((1/dpe)**2+(1/dbe)**2)/2
            except:
                edep=0

            if switch==0:
                #mea
                dmea=dnme.NDX[dnme.indx[name]]
                edmea=dnme.STDNDX[dnme.indx[name]]
                #pro
                dp=modelp.DX[modelp.indx[name]]/sqrt(modelp.BETX[modelp.indx[name]])
                db=modelb.DX[modelb.indx[name]]/sqrt(modelb.BETX[modelp.indx[name]])
                dep=(1/(dpe+dbe))*(dpe*dp+dbe*db)
                
                print >>filex,name,s,dmea,edmea,dmo,dep,edep,smo
            else:
                
                dep=modelcor.DX[modelcor.indx[name]]/sqrt(modelcor.BETX[modelcor.indx[name]])
                print >>filex,name,s,dmo,dep,edep,smo

            #if namename in name:

                #filed1=filed1+" "+str(dmo)+" "+str(dep)+" "+str(edep)



    #vertical dispersion
    if len(dispv)!=0:
        dme=dispv[0]
        dminp=dispv[1]
        dmaxp=dispv[2]
        dminb=dispv[3]
        dmaxb=dispv[4]
        filey= open(path+'/sbsDy_'+namename+'.out','w')



        if switch==0:
            bpms=intersect([dme,dminp,modelp,model,modelb])
            print >>filey,"* NAME  S  DY  STDDY  DY_MDL DY_PLAY EDY_PLAY DPY DPY_MDL DPY_PLAY EDPY_PLAY MODEL_S"
            print >>filey,"$ %s %le %le %le %le %le %le %le %le %le %le %le"
        else:
            bpms=intersect([dminp,modelp,model,modelb])
            print >>filey,"* NAME  S  DY_MDL DY_PLAY EDY_PLAY DPY_MDL DPY_PLAY EDPY_PLAY  MODEL_S"
            print >>filey,"$ %s %le %le %le %le %le %le %le %le"
    
    
        
        for bpm in bpms:
            s=bpm[0]
            name=bpm[1]

            #model
            dmo=model.DY[model.indx[name]]
            smo=model.S[model.indx[name]]
            dpmo=model.DPY[model.indx[name]]
            # pro
            dpe=abs(dmaxp.DY[dmaxp.indx[name]]-dminp.DY[dminp.indx[name]])
            dbe=abs(dmaxb.DY[dmaxb.indx[name]]-dminb.DY[dminb.indx[name]])
            dppe=abs(dmaxp.DPY[dmaxp.indx[name]]-dminp.DPY[dminp.indx[name]])
            dpbe=abs(dmaxb.DPY[dmaxb.indx[name]]-dminb.DPY[dminb.indx[name]])               
            edep=sqrt((dpe)**2+(dbe)**2)/2
            edepp=sqrt((dppe)**2+(dpbe)**2)/2
            

            if switch==0:
                #mea
                dmea=dme.DY[dme.indx[name]]
                edmea=dme.STDDY[dme.indx[name]]
                dpmea=dme.DPY[dme.indx[name]]
            
                #pro
                dep=modelcor.DY[modelcor.indx[name]]
                depp=modelcor.DPY[modelcor.indx[name]]

                print >>filey,name,s,dmea,edmea,dmo,dep,edep,dpmea,dpmo,depp,edepp,smo
            else:
                dp=modelp.DY[modelp.indx[name]]
                db=modelb.DY[modelb.indx[name]]
                dep=(1/(dpe+dbe))*(dpe*dp+dbe*db)
                
                dpp=modelp.DPY[modelp.indx[name]]
                dpb=modelb.DPY[modelb.indx[name]]
                depp=(1/(dppe+dpbe))*(dppe*dpp+dpbe*dpb)
                
                print >>filey,name,s,dmo,dep,edep,dpmo,depp,edepp,smo


                if namename in name:

                    if len(chromatic)!=0:
                        filed1=filed1+" "+str(dmo)+" "+str(dep)+" "+str(edep)
                    else:
                        print >>filesum_d,filed1,dmo,dep,edep,dpmo,depp,edepp,smo


 

        filey.close()    
    ##coupling
    couplemi=couple[2]
    couplema=couple[3]

    couplemi.Cmatrix()
    couplema.Cmatrix()
    model.Cmatrix()
    modelp.Cmatrix()
    modelcor.Cmatrix()
    #modelb.Cmatrix()
    


    #if (switch==0) and (len(couple)!=0):
        #bpms=intersect([coupleme,couplemi,model])

        #fterms
        #print >> filex,'* NAME   S   f1001 f1001re  f1001im    f1010   f1010re   f1010im  f1001_EXP f1001err_EXP  f1001re_EXP  f1001im_EXP    f1010_EXP  f1010err_EXP   f1010re_EXP   f1010im_EXP     f1001_PLAY ef1001_play   f1001re_PLAY  f1001im_PLAY    f1010_PLAY ef1010_play   f1010re_PLAY   f1010im_PLAY  S_MODEL'
        #print >> filex,'$  %s    %le   %le    %le  %le   %le    %le    %le  %le   %le    %le  %le   %le    %le    %le   %le  %le   %le    %le    %le   %le  %le %le %le %le'

        #cterms
        #print >>filex2,"NAME S C11Me eC11Me C12Me eC12Me C21Me eC21Me C22Me eC22 ANGMe eANGMe C11Mo C12Mo C21Mo C22Mo ANDMo C11_cor eC11_cor C12_cor eC12_cor C21_cor eC21_cor C22_cor eC22_cor ANG_cor eANG_cor"
        #print >>filex2,"%s %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le"
        
    if switch!=0:
        filex=open(path+'/sbscouple_'+namename+'.out','w')
        filex2=open(path+'/sbscouple2_'+namename+'.out','w')
        bpms=intersect([couplemi,model])

        #fterms
        print >> filex,'* NAME   S   f1001 f1001re  f1001im    f1010   f1010re   f1010im  f1001_PLAY ef1001_play   f1001re_PLAY  f1001im_PLAY    f1010_PLAY ef1010_play   f1010re_PLAY   f1010im_PLAY  S_MODEL'
        print >> filex,'$  %s    %le    %le  %le   %le    %le  %le    %le  %le   %le  %le   %le    %le    %le   %le  %le %le'

        #cterms
        print >>filex2,"NAME S C11Mo C12Mo C21Mo C22Mo ANDMo C11_cor eC11_cor C12_cor eC12_cor C21_cor eC21_cor C22_cor eC22_cor ANG_cor eANG_cor"
        print >>filex2,"%s %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le"        
    else:
        bpms=[]

    for bpm in bpms:
        s=bpm[0]
        name=bpm[1]

        #fterms
        f1001=model.f1001[model.indx[name]]
        f1010=model.f1010[model.indx[name]]
        f1001e=couplema.f1001[couplema.indx[name]]-couplemi.f1001[couplemi.indx[name]]
        f1010e=couplema.f1010[couplema.indx[name]]-couplemi.f1010[couplemi.indx[name]]

        #cterms
        C11=model.C[model.indx[name]][0]
        C12=model.C[model.indx[name]][1]
        C21=model.C[model.indx[name]][2]
        C22=model.C[model.indx[name]][3]
        angmo="x"

        eC11c=couplema.C[couplema.indx[name]][0]-couplemi.C[couplemi.indx[name]][0]
        eC12c=couplema.C[couplema.indx[name]][1]-couplemi.C[couplemi.indx[name]][1]
        eC21c=couplema.C[couplema.indx[name]][2]-couplemi.C[couplemi.indx[name]][3]
        eC22c=couplema.C[couplema.indx[name]][3]-couplemi.C[couplemi.indx[name]][3]
        
        smo=model.S[model.indx[name]]

        C11c=modelp.C[modelp.indx[name]][0]
        C12c=modelp.C[modelp.indx[name]][1]
        C21c=modelp.C[modelp.indx[name]][2] 
        C22c=modelp.C[modelp.indx[name]][3]


        angle="x"

        f1001c=modelp.f1001[modelp.indx[name]]
        f1010c=modelp.f1010[modelp.indx[name]]
        f1001e=couplema.f1001[couplema.indx[name]]-couplemi.f1001[couplemi.indx[name]]
        f1010e=couplema.f1010[couplema.indx[name]]-couplemi.f1010[couplemi.indx[name]]    
            #fterms
        print >>filex,name,s,abs(f1001),f1001.real,f1001.imag,abs(f1010),f1010.real,f1010.imag,abs(f1001c),f1001e,f1001c.real,f1001c.imag,abs(f1010c),f1010e,f1010c.real,f1010c.imag,smo

            #cterms
        print >>filex2,name,s,C11,C12,C21,C22,angmo,C11c,eC11c,C12c,eC12c,C21c,eC21c,C22c,eC22c,angle

        if namename in name:
            print >>filesum_c,name,s,abs(f1001),f1001.real,f1001.imag,abs(f1010),f1010.real,f1010.imag,abs(f1001c),f1001e,f1001c.real,f1001c.imag,abs(f1010c),f1010e,f1010c.real,f1010c.imag,C11,C12,C21,C22,angmo,C11c,eC11c,C12c,eC12c,C21c,eC21c,C22c,eC22c,"x","X",smo
        

    if switch!=0:
        filex.close()
        filex2.close()
    
    ##chromatic beta
    if len(chromatic)!=0:

        # horizontal
        chromame=chromatic[0]
        chromamip=chromatic[1]
        chromamap=chromatic[2]

        filex= open(path+'/sbsWx_'+namename+'.out','w')

        if (switch==0):
            bpms=intersect([chromame,chromamip,model,modelp])

            print >> filex,"* NAME  S  WX   WXERR WX_MDL WX_PLAY eWX_play  PHIX PHIXERR PHIX_MDL PHIX_PLAY ePHIX_PLAY    MODEL_S "
            print >> filex,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "
        else:
            bpms=intersect([chromamip,model,modelp])

            print >> filex,"* NAME  S  WX_MDL WX_PLAY eWX_play PHIX_MDL PHIX_PLAY ePHIX_PLAY    MODEL_S "
            print >> filex,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "


        for bpm in bpms:
            s=bpm[0]
            name=bpm[1] 

            wmo=model.WX[model.indx[name]]
            phmo=model.PHIX[model.indx[name]]            
            smo=model.S[model.indx[name]]

            we=chromamap.WX[chromamap.indx[name]]-chromamip.WX[chromamip.indx[name]]
            phe=chromamap.PHIX[chromamap.indx[name]]-chromamip.PHIX[chromamip.indx[name]]            

            if switch==0:
                wme=chromame.WX[chromame.indx[name]]
                ewme=chromame.WXERR[chromame.indx[name]]
                phme=chromame.PHIX[chromame.indx[name]]
                ephme=chromame.PHIXERR[chromame.indx[name]]
                
                w_cor=modelcor.WX[modelcor.indx[name]]
                ph_cor=modelcor.PHIX[modelcor.indx[name]]

                print >> filex,name,s,wme,ewme,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo
                
            else:
                
                w_cor=modelp.WX[modelp.indx[name]]
                ph_cor=modelp.PHIX[modelp.indx[name]]

                print >> filex,name,s,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

                if namename in name:

                    filed1=filed1+" "+str(wmo)+" "+str(w_cor)+" "+str(we)+" "+str(phme)+" "+str(ephme)+" "+str(phmo)+" "+str(ph_cor)+" "+str(phe)

        # vertical
        chromame=chromatic[3]
        chromamip=chromatic[4]
        chromamap=chromatic[5]

        filey= open(path+'/sbsWy_'+namename+'.out','w')

        if switch==0:
            bpms=intersect([chromame,chromamip,model,modelp])

            print >> filey,"* NAME  S  WY   WYERR WY_MDL WY_PLAY eWY_play  PHIY PHIYERR PHIYMDL PHIY_PLAY ePHIY_PLAY    MODEL_S "
            print >> filey,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "
        else:
            bpms=intersect([chromamip,model,modelp])

            print >> filey,"* NAME  S  WY_MDL WY_PLAY eWY_play PHIY_MDL PHIY_PLAY ePHIY_PLAY    MODEL_S "
            print >> filey,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "

        for bpm in bpms:
            s=bpm[0]
            name=bpm[1]

            wmo=model.WY[model.indx[name]]
            phmo=model.PHIY[model.indx[name]]            
            smo=model.S[model.indx[name]]

            we=chromamap.WY[chromamap.indx[name]]-chromamip.WY[chromamip.indx[name]]
            phe=chromamap.PHIY[chromamap.indx[name]]-chromamip.PHIY[chromamip.indx[name]]            

            if switch==0:
                wme=chromame.WY[chromame.indx[name]]
                ewme=chromame.WYERR[chromame.indx[name]]
                phme=chromame.PHIY[chromame.indx[name]]
                ephme=chromame.PHIYERR[chromame.indx[name]]
                
                w_cor=modelcor.WY[modelcor.indx[name]]
                ph_cor=modelcor.PHIY[modelcor.indx[name]]

                print >> filey,name,s,wme,ewme,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo
                
            else:
                
                w_cor=modelp.WY[modelp.indx[name]]
                ph_cor=modelp.PHIY[modelp.indx[name]]

                print >> filey,name,s,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

                if namename in name:

                    print >>filesum_d,filed1,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

        filesum_b.close()
        filesum_c.close()
        filesum_d.close()
                
        filey.close()

    ## to find IP
    if ("IP" in namename) and (switch==1):
        getIP([betah[9],betav[9]],[modelp,modelb],[betah[1],betah[3],betah[2],betah[4]],[betah[5],betah[7],betah[6],betah[8]],model,phasex,phasey,namename,accel,path)
    elif ("ADT" in namename) and (switch==1):
        errors=[betah[1],betah[3],betah[2],betah[4]]
        TransverseDampers(modelp,modelb,namename,model,path,phases[0],phases[1],errors)
        

def run4mad(path,hor,ver,hore,vere,dp,dpe,startbpm,endbpm,name, fs, exppath,twisspath, method):


    ##
    #  Copy getfterms.py locally to be used by MADx
    cpath=options.bb
    os.system('cp '+cpath+'/SegmentBySegment/getfterms_0.3.py'+' '+path+'/')        
    
    # Copy the modifiers.madx file locally to be used by MADx
    
    if os.path.isfile(twisspath+'/modifiers.madx'):
        os.system('cp '+twisspath+'/modifiers.madx'+' '+path+'/')
    else :   #If the modifiers file does not exist create empty file 
        os.system('touch '+path+'/modifiers.madx')
    
    
    if options.accel=="LHCB2":
        dire=-1
        start="MKI.A5R8.B2"
        beam="B2"
        
    elif options.accel=="LHCB1":
        dire=1
        start="MSIA.EXIT.B1"   #  compatible with repository
        beam="B1"

    
    ### Check chromatic functions in wpath
    if options.wpath!="0":
        print "Chrom path,", options.wpath
        try:
            wx=twiss(options.wpath+"/wx.out")   
            wy=twiss(options.wpath+"/wy.out")
        except:
            print "No Chromatic functions (wx,wy) available at ",options.wpath
            options.wpath="0" # This means set wx,wy=0,0
        
    if options.wpath!="0":
        try:
            wxs=wx.WX[wx.indx[startbpm]]
            phixs=wx.PHIX[wx.indx[startbpm]]
        except:
            print "Start BPM, ",startbpm," not in WX file"
            wxs=0
            phixs=0

        try:
            wys=wy.WY[wy.indx[startbpm]]
            phiys=wy.PHIY[wy.indx[startbpm]]
        except:
            print "Start BPM, ",startbpm," not in WY file"
            wys=0
            phiys=0
    else:
        wxs=0
        phixs=0
        wys=0
        phiys=0

    ### check on error propogation
    errbetx=hor[1]
    betx=hor[0]
    errbety=ver[1]
    bety=ver[0]
    errbetxb=hore[1]
    betxb=hore[0]
    errbetyb=vere[1]
    betyb=vere[0]    
    f1001r=fs[0]
    f1001i=fs[1]
    f1010r=fs[2]
    f1010i=fs[3]
    f1001std=fs[4]

    madfilename=path+'/t_'+str(name)+'.madx'
    
    filename=path+'/var4mad.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%BETX/\''+str(betx)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETY/\''+str(bety)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETX/\''+str(errbetx)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETY/\''+str(errbety)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFX/\''+str(hor[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFY/\''+str(ver[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFX/\''+str(hor[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFY/\''+str(ver[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DX/\''+str(dp[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDX/\''+str(dp[4])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DY/\''+str(dp[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDY/\''+str(dp[5])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPX/\''+str(dp[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPY/\''+str(dp[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDBX/\''+str(betxb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDBY/\''+str(betyb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRENDBX/\''+str(errbetxb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRENDBY/\''+str(errbetyb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFENDX/\''+str(-hore[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFENDY/\''+str(-vere[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFENDX/\''+str(hore[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFENDY/\''+str(vere[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DENDX/\''+str(dpe[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDENDX/\''+str(dpe[4])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DENDY/\''+str(dpe[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDENDY/\''+str(dpe[5])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPENDX/\''+str(-dpe[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPENDY/\''+str(-dpe[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%STARTFROM/\''+str(startbpm.replace("-","_"))+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDAT/\''+str(endbpm.replace("-","_"))+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(name)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DIRE/\''+str(dire)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BEAM/\''+str(beam)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PATH/\'\"'+str(path.replace('/','\/'))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001R/\''+str(f1001r)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001I/\''+str(f1001i)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010R/\''+str(f1010r)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010I/\''+str(f1010i)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001maR/\''+str(f1001r+f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001maI/\''+str(f1001i+f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001miR/\''+str(f1001r-f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001miI/\''+str(f1001i-f1001std)+'\'/g\' \\\n')    
    file4nad.write('    -e \'s/%F1010maR/\''+str(f1010r+f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010maI/\''+str(f1010i+f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010miR/\''+str(f1010r-f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010miI/\''+str(f1010i-f1001std)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%METHOD/\''+method+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%EXP/\'\"'+exppath.replace("/","\/")+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%WX/\''+str(wxs)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PHIX/\''+str(phixs)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%WY/\''+str(wys)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PHIY/\''+str(phiys)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%WPATH/\'\"'+str(options.wpath.replace('/','\/'))+'\"\'/g\' \\\n')
    
    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/job.InterpolateBetas.0_2_dev.mask > '+madfilename+'\n')

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))
    runmad(path,name)


    #Prepare watchdog file command in directory
    watchfilename=path+"/watch_"+str(name)
    fwatch=open(watchfilename,"w")
    print >> fwatch, "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/SegmentBySegment/watch.py "+madfilename+" "+path+"/gplot_"+str(name)
    fwatch.close()
    os.system("chmod +x "+watchfilename)
    
    
    #sys.exit()

   
    
def runmad(path,name):
    
    os.system(options.mad+' < '+path+'t_'+str(name)+'.madx')
    
   
  
def run4plot(path,spos,epos,beta4plot,cpath,meapath,name,qx,qy,accel,method):
    if method=="driven": method=""   # patch to make it work at inj. rogelio
    filename=path+'/var4plot.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%PATH/\'\"'+str(path.replace("/","\/"))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%EndPoint/\''+str(epos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%StartPoint/\''+str(spos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(name)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETA/\''+str(beta4plot)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%QX/\''+str(qx)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%QY/\''+str(qy)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%METHOD/\''+method+'\'/g\' \\\n')  
    file4nad.write('    -e \'s/%MEA/\'\"'+str(meapath.replace("/","\/"))+'\"\'/g\' \\\n')
    if (name=="IP8" and accel=="LHCB2") or (name=="IP2" and accel=="LHCB1"):
        file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.IP2IP8.0_3.mask > '+path+'/gplot_'+name)
    elif "RHIC" in options.accel:
        file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.0_1_RHIC.mask > '+path+'/gplot_'+name)
    else:
        file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.0_3.mask > '+path+'/gplot_'+name)
        

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))
    os.system("gnuplot "+str(path+'/gplot_'+name))

 
   
def GetPhaseEM(exp, mod):
    phasem=[]
    bpm1=[]
    bpm2=[]
    s1=[]
    s2=[]
    phaseexp=[]
    modName=mod.NAME
    names=[]

    for name in modName:

        if "BPM" in name:

            names.append(name)

    
    for elm in names:
        status=1
        try:
            expind=exp.indx[elm]
        except:
            print elm, "not found in exp"
            status=0

        if status==1:
            el2=exp.NAME2[exp.indx[elm]]
            elmind=mod.indx[elm]
            try:
                elm2ind=mod.indx[el2]
            except:
                print el2, "not found in model"
                status=0

            if status==1:
                if "PHXMDL" in  exp.__dict__.keys():
                    modphaseadv=mod.MUX[elm2ind]-mod.MUX[elmind]
                    phaseexp.append(exp.PHASEX[expind])
                    
                elif "PHYMDL" in  exp.__dict__.keys():
                    modphaseadv=mod.MUY[elm2ind]-mod.MUY[elmind]
                    phaseexp.append(exp.PHASEY[expind])
                    
                    
                bpm1.append(elm)
                bpm2.append(mod.NAME[elm2ind])
                s1.append(mod.S[elmind])
                s2.append(mod.S[elm2ind])
                phasem.append(modphaseadv)
                
                
    return bpm1, bpm2, s1, s2, phaseexp, phasem

def writePhase(filename,bpm1, bpm2, s1, s2, phaseexp, phasem ):

    
    f=open(filename, "w")
    for i in range(len(bpm1)):
        print >>f, bpm1[i], bpm2[i], s1[i], s2[i], phaseexp[i], phasem[i]

    f.close()

#delete
def reversetable(path,name):
    newFile=open(path+"/twiss_"+name+"_back_rev.dat",'w')
    base=twiss(path+"/twiss_"+name+"_back.dat")

    newFile.write("* NAME                                S               BETX               ALFX               BETY               ALFY    DX       DPX     DY     DPY   MUX   MUY\n")
    newFile.write("$ %s                                %le                %le                %le                %le                %le      %le                %le      %le   %le            %le                %le \n")
    
    bpms=base.NAME
    s=base.S
    endpos=base.S[base.indx[bpms[len(bpms)-1]]]
    betx=base.BETX
    bety=base.BETY
    alfx=base.ALFX
    alfy=base.ALFY
    dx=base.DX
    dpx=base.DPX
    dy=base.DY
    dpy=base.DPY
    mux=base.MUX
    muy=base.MUY

    for i in range(len(bpms)):
        index=len(bpms)-1-i
        bpm=bpms[index]
        bex=betx[index]
        bey=bety[index]
        alx=alfx[index]
        aly=alfy[index]
        dex=dx[index]
        depx=dpx[index]
        dey=dy[index]
        depy=dpy[index]
        ss=s[index]
        muxx=mux[index]
        muyy=muy[index]

        #print bpm, dex

        newFile.write(str(bpm)+' '+str(endpos-ss)+' '+str(bex)+' '+str(alx)+' '+str(bey)+' '+str(aly)+'  '+str(dex)+' '+str(depx)+' '+str(dey)+' '+str(depy)+' '+str(muxx)+' '+str(muyy)+'\n')
        
    newFile.close()




#delete
def createTables(outputname,path,columnnames,paranames,data,mainvariable,mainvalue):

    filefile=open(path+"/"+outputname,'w')

    # writing main variables
    if len(mainvariable)==len(mainvalue):
        for count in range(len(mainvariable)):
            
            if isinstance(mainvalue[count],str):
                filefile.write('@ '+mainvariable[count]+' %s '+str(mainvalue[count])+'\n')
            
            else:
                filefile.write('@ '+mainvariable[count]+' %le '+str(mainvalue[count])+'\n')
    else:
        print "cannot write main variables ..."

    # writing columnnames and paranames

    print len(columnnames),len(paranames),outputname
    

    if len(columnnames)==len(paranames):


        filefile.write('* '+'  '.join(columnnames)+'\n')
        filefile.write('$ '+'  '.join(paranames)+'\n')

        #writing data
        #
        # data[x][y]
        #
        #   x = data set
        #   y = data in data set
        #   => has to bed added in same order
        #
        
        for y in range(len(data)):
        
            filefile.write(str(data[y])+'\n')


    
    
    else:

        print "cannot create table names for columns are not equal => system exit "
        sys.exit()

    filefile.close()




#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":
    (options, args) = parse_args()
    
    return_value = main(options)
    
    sys.exit(return_value)