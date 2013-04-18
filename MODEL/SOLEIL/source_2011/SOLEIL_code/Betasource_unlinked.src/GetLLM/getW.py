from metaclass import *
from optparse import OptionParser
from linreg import linreg
from math import sqrt, atan, atan2, pi
import sys

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




parser = OptionParser()
#parser.add_option("-a", "--accel",
#                help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
#                metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-m", "--model",
                help="Twiss File",
                metavar="TwissFile", default="0", dest="Twiss")
parser.add_option("-f", "--files",
                help="paths to the results dirs separated by comma",
                metavar="FILES", default="0", dest="files")
parser.add_option("-p", "--dpps",
                help="dp/p's for the results dirs separated by comma and in the same order",
                metavar="DPPS", default="0", dest="dpps")
parser.add_option("-u", "--dpunit",
                help="dpp unit, defaul 0.001",
                metavar="DPUNIT", default="0.001", dest="dpunit")
parser.add_option("-o", "--output",
                help="output path, defaul ./",
                metavar="OUT", default="./", dest="out")


(options, args) = parser.parse_args()



m=twiss(options.Twiss)
dpps=options.dpps.split(",")
for i in range(len(dpps)):
    dpps[i]=float(dpps[i])*float(options.dpunit)
files=options.files.split(",")


#p0=twiss("4:44:25/getbetax.out")
#p1=twiss("6:27:41/getbetax.out")
#p_1=twiss("6:11:56/getbetax.out")
#m=twiss("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB2/acdipole.opt/twiss.dat")


m.chrombeat()

tfilesx=[]
tfilesy=[]
for el in files:
    tfilesx.append(twiss(el+"/getbetax.out"))
    tfilesy.append(twiss(el+"/getbetay.out"))

# Is there a dpp=zero?
zeroi=-1
for i in range(len(dpps)):
    if float(dpps[i])==0.0:
        zeroi=i
        tzerox=tfilesx[i]
        tzeroy=tfilesy[i]
        print "Found dpp=0 for case:", files[i]
if zeroi<0:
    print "dpp=0 not found, better stop"
    sys.exit()


# From MAD manual
#WX = Wx = sqrt(ax2 + bx2),
#ax = (del betax / del PT) / betax,
#bx = (del alphax / del PT) - (alphax / betax) * (del betax / del PT).
# PHIX: Chromatic phase function Phix, [2pi]:
#PHIX = Phix = atan(ax / bx). 

f=open(options.out+"wx.out","w")
print >>f, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WX", "WXERR", "PHIX", "PHIXERR"
print >>f, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le"

bpms=intersect(tfilesx)
for bpm in bpms:
    el=bpm[1]
    indx=[]
    b=[]
    a=[]
    betax0=tzerox.BETX[tzerox.indx[el]]
    alfax0=tzerox.ALFX[tzerox.indx[el]]
    alfax0err=tzerox.STDALFX[tzerox.indx[el]]
    for file in tfilesx:
        ix=file.indx[el]
        indx.append(ix)
        #b.append((file.BETX[ix]-betax0)/betax0)
        b.append(file.BETX[ix])
        a.append(file.ALFX[ix])
    bfit=linreg(dpps, b)
    afit=linreg(dpps, a)
    dbb=bfit[0]/betax0
    dbberr=bfit[3]/betax0
    da=afit[0]
    daerr=afit[3]
    AX=dbb     # FROM and FOR MAD
    AXerr=dbberr
    BX=da-alfax0*dbb       # FROM and FOR MAD 
    BXerr=sqrt(daerr**2 + (alfax0err*dbb)**2 + (alfax0*dbberr)**2)
    wx=sqrt(AX**2+BX**2)
    wxerr=sqrt( (AXerr*AX/wx)**2 + (BXerr*BX/wx)**2  )
    phix=atan2(BX,AX)/2./pi
    phixerr=1./(1.+(AX/BX)**2)*sqrt( (AXerr/BX)**2 + (AX/BX**2*BXerr)**2)/2./pi
    print >>f, file.NAME[ix], file.S[ix],  dbb, dbberr, da, daerr, wx, wxerr, phix, phixerr
    

#        print >>f, el, p0.S[p0i], p0.BETX[p0i],p0.STDBETX[p0i], p1.BETX[p1i],p1.STDBETX[p1i], m.dbx[mi], p_1.BETX[p_1i],p_1.STDBETX[p_1i]
#        print >>g, el, p0.S[p0i], p0.ALFX[p0i],p0.STDALFX[p0i], p1.ALFX[p1i],p1.STDALFX[p1i], m.dbx[mi], p_1.ALFX[p_1i],p_1.STDALFX[p_1i]
        

f.close()





f=open(options.out+"wy.out","w")
print >>f, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WY", "WYERR",    "PHIY", "PHIYERR"
print >>f, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le"


bpms=intersect(tfilesy)
for bpm in bpms:
    el=bpm[1]
    indx=[]
    b=[]
    a=[]
    betay0=tzeroy.BETY[tzeroy.indx[el]]
    alfay0=tzeroy.ALFY[tzeroy.indx[el]]
    alfay0err=tzeroy.STDALFY[tzeroy.indx[el]]
    for file in tfilesy:
        ix=file.indx[el]
        indx.append(ix)
        #b.append((file.BETX[ix]-betax0)/betax0)
        b.append(file.BETY[ix])
        a.append(file.ALFY[ix])
    bfit=linreg(dpps, b)
    afit=linreg(dpps, a)
    dbb=bfit[0]/betay0
    dbberr=bfit[3]/betay0
    da=afit[0]
    daerr=afit[3]
    AX=dbb     # FROM and FOR MAD
    AXerr=dbberr
    BX=da-alfay0*dbb       # FROM and FOR MAD 
    BXerr=sqrt(daerr**2 + (alfay0err*dbb)**2 + (alfay0*dbberr)**2)
    wy=sqrt(AX**2+BX**2)
    wyerr=sqrt( (AXerr*AX/wy)**2 + (BXerr*BX/wy)**2  )
    phiy=atan2(BX,AX)/2./pi
    phiyerr=1./(1.+(AX/BX)**2)*sqrt( (AXerr/BX)**2 + (AX/BX**2*BXerr)**2)/2./pi
    print >>f, file.NAME[ix], file.S[ix],  dbb, dbberr, da, daerr, wy, wyerr, phiy, phiyerr
    

#        print >>f, el, p0.S[p0i], p0.BETX[p0i],p0.STDBETX[p0i], p1.BETX[p1i],p1.STDBETX[p1i], m.dbx[mi], p_1.BETX[p_1i],p_1.STDBETX[p_1i]
#        print >>g, el, p0.S[p0i], p0.ALFX[p0i],p0.STDALFX[p0i], p1.ALFX[p1i],p1.STDALFX[p1i], m.dbx[mi], p_1.ALFX[p_1i],p_1.STDALFX[p_1i]
        

f.close()






sys.exit()

