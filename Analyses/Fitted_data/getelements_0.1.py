#############################################################################################
#
# Glenn Vanbavinckhove (gvanbavi@cern.ch)
#
#############################################################################################

# => getelements_0.1.py :  - Creation of the program
#                         
#
#

####### imports
from metaclass import twiss
from os import system
import os,sys
from math import sqrt
from optparse import OptionParser


####### optionparser
parser = OptionParser()
parser.add_option("-m", "--model",
                help="Which model to use",
                metavar="twiss", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB1/twiss_elements.dat",dest="twiss")
parser.add_option("-a", "--accel",
                help="Which accel to use",
                metavar="accel", default="./", dest="accel")
parser.add_option("-o", "--output",
                help="Output files",
                metavar="out", default="./", dest="out")
parser.add_option("-i", "--input",
                help="in file (changeparameters)",
                metavar="inn", default="./", dest="inn")
parser.add_option("-b", "--bsource",
                help="beta-beat path",
                metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src", dest="bb")
(options, args) = parser.parse_args()



def runmad(output,corpath,cpath):

    if options.accel=="LHCB2":

        dire=-1
        start="MKI.A5R8.B2"
        beam="B2"

    elif options.accel=="LHCB1":

        dire=1
        start="MSIA.EXIT.B1"
        beam="B1"


    file4mad=open(output+"/madx.sh","w")
    file4mad.write('sed -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4mad.write('    -e \'s/%BV/\''+str(dire)+'\'/g\' \\\n')
    file4mad.write('    -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
    file4mad.write('    -e \'s/%CORPATH/\'\"'+corpath.replace("/","\/")+'\"\'/g\' \\\n')
    file4mad.write('    -e \'s/%PATH/\'\"'+output.replace("/","\/")+'\"\'/g\' \\\n')
    file4mad.write('<'+cpath+'/Analyses//Fitted_data//'+'/job_0_1.mask > '+output+'/fitted.madx \n')

    file4mad.close()

    print "Starting madx instance"
    
    os.system("chmod 777 "+output+"/madx.sh")
    os.system(output+"/madx.sh")

    os.system("madx <"+output+'/fitted.madx')

    #sys.exit()

    print "madx instance finished"


def modelIntersect(expbpms, model):
    bpmsin=[]
        
    for bpm in expbpms:

        #print bpm[1]            
        try:
            check=model.indx[bpm]
            bpmsin.append(bpm)
        except:
            print bpm, "Not in Model"
    if len(bpmsin)==0:
        print "Zero intersection of Exp and Model"
        print "Please, provide a good Dictionary"
        print "Now we better leave!"
        sys.exit()

    return bpmsin


def gatheringdata(path,model):

    
    fit=twiss(path+'/fitted_twiss.dat')

    outputfile=open(path+'/summary_twiss_bet.dat','w')
    outputfile_phase=open(path+'/summary_twiss_pha.dat','w')

    # write tunes

    
    # main output
    print >> outputfile,"* NAME S BETX BETXMDL  ERRBETX BETY BETYMDL ERRBETY  ALFX ALFXMDL ERRALFX ALFY ALFYMDL ERRALFY"
    print >> outputfile,"$ %s   %le %le %le    %le  %le     %le  %le     %le  %le  %le  %le     %le  %le"

    print >> outputfile_phase,"* NAME NAME1 S S1 PHX PHXMDL  ERRPX  PHY   PHYMDL  ERRPY"
    print >> outputfile_phase,"$ %s   %s  %le %le  %le   %le  %le  %le   %le  %le"
    
    bpms=modelIntersect(fit.NAME, model)

    for i in range(len(bpms)): # for beta


        bps=fit.S[fit.indx[bpms[i]]]
        bpm=bpms[i]

        errbx=(fit.BETX[fit.indx[bpm]]**2+model.BETX[model.indx[bpm]]**2)/2
        errbx=sqrt(errbx-((fit.BETX[fit.indx[bpm]]+model.BETX[model.indx[bpm]])/2)**2)#+2.2e-15
        errby=(fit.BETY[fit.indx[bpm]]**2+model.BETY[model.indx[bpm]]**2)/2
        errby=sqrt(errby-((fit.BETY[fit.indx[bpm]]+model.BETY[model.indx[bpm]])/2)**2)#+2.2e-15
        errax=(fit.ALFX[fit.indx[bpm]]**2+model.ALFX[model.indx[bpm]]**2)/2
        errax=sqrt(errax-((fit.ALFX[fit.indx[bpm]]+model.ALFX[model.indx[bpm]])/2)**2)#+2.2e-15
        erray=(fit.ALFY[fit.indx[bpm]]**2+model.ALFY[model.indx[bpm]]**2)/2
        erray=sqrt(erray-((fit.ALFY[fit.indx[bpm]]+model.ALFY[model.indx[bpm]])/2)**2)#+2.2e-15

        print >> outputfile,bpm,bps,fit.BETX[fit.indx[bpm]],model.BETX[model.indx[bpm]],errbx,fit.BETY[fit.indx[bpm]],model.BETY[model.indx[bpm]],errby,fit.ALFX[fit.indx[bpm]],model.ALFX[model.indx[bpm]],errax,fit.ALFY[model.indx[bpm]],model.ALFY[model.indx[bpm]],erray

    outputfile.close()


    for i in range(len(bpms)-1): # for phase

        bps=fit.S[fit.indx[bpms[i]]]
        bps1=fit.S[fit.indx[bpms[i+1]]]
        bpm=bpms[i]
        bpm1=bpms[i+1]


        phfitx=-fit.MUX[fit.indx[bpm]]+fit.MUX[fit.indx[bpm1]]
        phfitx_mdl=-model.MUX[model.indx[bpm]]+model.MUX[model.indx[bpm1]]

        phfity=-fit.MUY[fit.indx[bpm]]+fit.MUY[fit.indx[bpm1]]
        phfity_mdl=-model.MUY[model.indx[bpm]]+model.MUY[model.indx[bpm1]]
        #print ((phfitx**2+phfitx_mdl**2)/2),((phfitx+phfitx_mdl)**2)/2
        errpx=sqrt(((phfitx**2+phfitx_mdl**2)/2)-(((phfitx+phfitx_mdl)/2)**2))
        errpy=sqrt(((phfity**2+phfity_mdl**2)/2)-(((phfity+phfity_mdl)/2)**2))
    

        print >> outputfile_phase, bpm,bpm1,bps,bps1,phfitx,phfitx_mdl,errpx,phfity,errpy,phfity,phfity_mdl


    outputfile_phase.close()


print "Entering mad instance"
runmad(options.out,options.inn,options.bb)
print "Gathering data"
gatheringdata(options.out,twiss(options.twiss))
print "Finished !!"
