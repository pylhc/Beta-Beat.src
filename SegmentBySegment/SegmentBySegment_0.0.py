################################################################
#                                                              #
#  @ Glenn Vanbavinckhove  => Date: 11/09/09                   #
#                                                              #
################################################################
#
#  !=> SegementBySegment_0.0.py : Construction of the program
#
#
#

###### imports
from optparse import OptionParser
from metaclass import twiss
import os




###### optionparser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                metavar="ACCEL", default="LHCB1",dest="accel")
parser.add_option("-d", "--dictionary",
                help="File with the BPM dictionary",
                metavar="DICT", default="0", dest="dict")
parser.add_option("-p", "--path", # assumes that output is same as input
                help="Path to measurement files",
                metavar="PATH", default="./", dest="path")
parser.add_option("-s", "--start",
                help="At which element to start",
                metavar="START", default="./", dest="start")
parser.add_option("-e", "--end",
                help="Where to stop/ find the values",
                metavar="END", default="./", dest="end")
parser.add_option("-l", "--label",
                help="label ofr output",
                metavar="LABEL", default="./", dest="label")

(options, args) = parser.parse_args()



####### some defs
def getValues(path,element):

    filedatax=twiss(path+"getbetax.out")
    filedatay=twiss(path+"getbetay.out")

    betax=filedatax.BETX[filedatax.indx[element]]
    errbetax=filedatax.ERRBETX[filedatax.indx[element]]
    alfx=filedatax.ALFX[filedatax.indx[element]]
    erralfx=filedatax.ERRALFX[filedatax.indx[element]]
    
    betay=filedatay.BETY[filedatay.indx[element]]
    errbetay=filedatay.ERRBETY[filedatay.indx[element]]
    alfy=filedatay.ALFY[filedatay.indx[element]]
    erralfy=filedatay.ERRALFY[filedatay.indx[element]]

    hor=[betax,errbetax,alfx,erralfx]
    ver=[betay,errbetay,alfy,erralfy]

    return[hor,ver]

def getTwiss(file,element):

    filedatax=twiss(file)
  

    betax=filedatax.BETX[filedatax.indx[element]]
    alfx=filedatax.ALFX[filedatax.indx[element]]
    sx=filedatax.S[filedatax.indx[element]]

    betay=filedatax.BETY[filedatax.indx[element]]
    alfy=filedatax.ALFY[filedatax.indx[element]]
    sy=filedatax.S[filedatax.indx[element]]
 
    hor=[betax,alfx,sx]
    ver=[betay,alfy,sy]

    return[hor,ver]

def run4mad(path,hor,ver):

    if options.accel=="LHCB2":

        dire=-1
        start="MKI.A5R8.B2"
        beam="B2"
        
    elif options.accel=="LHCB1":

        dire=1
        start="MKI.A5L2.B1"
        beam="B1"

    cpath=os.getcwd()
     
    filename=path+'/var4mad.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%BETX/\''+str(hor[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETY/\''+str(ver[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETX/\''+str(hor[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETY/\''+str(ver[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFX/\''+str(hor[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFY/\''+str(ver[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFX/\''+str(hor[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFY/\''+str(ver[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%STARTFROM/\''+str(options.start)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDAT/\''+str(options.end)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(options.label)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DIRE/\''+str(dire)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BEAM/\''+str(beam)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PATH/\''+str(path.replace("/",""))+'\'/g\' \\\n')
    file4nad.write('<'+cpath+'/job.InterpolateBetas.mask > '+path+'/t_'+str(options.label)+'.madx \n')

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system("./"+str(filename))
    print "finished"
    
    os.system('madx < '+path+'t_'+str(options.label)+'.madx')
   
  
def run4plot(path,spos,epos):

    cpath=os.getcwd()

    filename=path+'/var4plot.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%PATH/\''+str(path.replace("/",""))+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%EndPoint/\''+str(epos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%StartPoint/\''+str(spos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(options.label)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('<'+cpath+'/gplot.mask > '+path+'/gplot_'+options.label)

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system("./"+filename)

    
   
def GetPhaseEM(exp, mod):
    phasem=[]
    bpm1=[]
    bpm2=[]
    s1=[]
    s2=[]
    phaseexp=[]
    for elm in mod.NAME:
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
                    
                elif "PHYMDL" in  exp.__dict__.keys():
                    modphaseadv=mod.MUY[elm2ind]-mod.MUY[elmind]
                    
                    
                bpm1.append(elm)
                bpm2.append(mod.NAME[elm2ind])
                s1.append(mod.S[elmind])
                s2.append(mod.S[elm2ind])
                phasem.append(modphaseadv)
                phaseexp.append(exp.PHASE[expind])
                
    return bpm1, bpm2, s1, s2, phaseexp, phasem

def writePhase(filename,bpm1, bpm2, s1, s2, phaseexp, phasem ):

    
    f=open(filename, "w")
    for i in range(len(bpm1)):
        print >>f, bpm1[i], bpm2[i], s1[i], s2[i], phaseexp[i], phasem[i]

    f.close()



####### main part
path=options.path
element=options.start

[hor,ver]=getValues(path,element)
run4mad(path,hor,ver)

##################
[hor,ver]=getTwiss(path+"StartPoint.twiss",element)
startpos=hor[2]
[hor,ver]=getTwiss(path+"twiss_"+str(options.label)+".dat",element)
betx=hor[0]
bety=ver[0]
ax=hor[1]
ay=ver[1]
[hor,ver]=getTwiss(path+"twiss_"+str(options.label)+".dat",options.end)
endpos=hor[2]

[hor,ver]=getTwiss(path+"twiss.b+.dat",element)
bplusx=hor[0]
bplusy=ver[0]
[hor,ver]=getTwiss(path+"twiss.b-.dat",element)
bminx=hor[0]
bminy=ver[0]

[hor,ver]=getTwiss(path+"twiss.a+.dat",element)
aplusx=hor[1]
aplusy=ver[1]
[hor,ver]=getTwiss(path+"twiss.a-.dat",element)
aminx=hor[1]
aminy=ver[1]

berrx=bplusx-bminx
aerrx=aplusx-aminx
berry=bplusy-bminy
aerry=aplusy-aminy



##################


phx=twiss(path+'/getphasex.out')
phy=twiss(path+'/getphasey.out')

m=twiss(path+"twiss_"+options.label+".dat")

bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)
writePhase(path+"/phasexEM.out",bpm1, bpm2, s1, s2, phaseexp, phasem)
bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)
writePhase(path+"/phaseyEM.out",bpm1, bpm2, s1, s2, phaseexp, phasem)

m=twiss(path+"twiss_"+options.label+"_play.dat")

bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)
writePhase(path+"/phasexEM_play.out",bpm1, bpm2, s1, s2, phaseexp, phasem)
bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)
writePhase(path+"/phaseyEM_play.out",bpm1, bpm2, s1, s2, phaseexp, phasem)





######################
run4plot(path,startpos,endpos)


############
fileresul=open(path+'/resul_'+options.label+'.tfs','w')
fileresul.write('* NAME   S     BETX    ERRBETX    ALFX   ERRALFX   BETY   ERRBETY  ALFY   ERRALFX\n')
fileresul.write('* %s     %le    %le   %le         %le    %le       %le    %le      %le    %le\n')
fileresul.write(str(options.end)+' '+str(endpos)+' '+str(betx)+' '+str(berrx)+' '+str(ax)+' '+str(aerrx)+' '+str(bety)+' '+str(berry)+' '+str(ay)+' '+str(aerry))
fileresul.close()





