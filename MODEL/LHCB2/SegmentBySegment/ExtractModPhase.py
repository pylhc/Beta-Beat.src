from metaclass import *

from Numeric import *




##############################
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

############################


from optparse import OptionParser
parser = OptionParser()
parser.add_option("-p", "--path",
                help="Path to GetLLM output",
                metavar="PATH", default="/afs/cern.ch/user/r/rtomas/lintrack/Beta-Beat.src/SVD",dest="PATH")
parser.add_option("-m", "--model",    
               help="Twiss File",
                metavar="TwissFile", default="0", dest="Twiss")


(options, args) = parser.parse_args()


path=options.PATH
tw=options.Twiss

print "Path:",path, "Model:",tw

phx=twiss(path+'/getphasex.out')
phy=twiss(path+'/getphasey.out')
m=twiss(tw)

bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)


f=open("phasexEM.out", "w")
for i in range(len(bpm1)):
    print >>f, bpm1[i], bpm2[i], s1[i], s2[i], phaseexp[i], phasem[i]

f.close()


bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)


f=open("phaseyEM.out", "w")
for i in range(len(bpm1)):
    print >>f,bpm1[i], bpm2[i], s1[i], s2[i], phaseexp[i], phasem[i]

f.close()
