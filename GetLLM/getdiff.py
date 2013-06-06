from metaclass import *
from optparse import OptionParser
import os

parser = OptionParser()


parser.add_option("-p", "--path",
    help="Path to measurement",
    metavar="PATH", default="0",dest="path")


parser.add_option("-r", "--reference",
        help="Path to the reference measurement",
        metavar="REF", default="0",dest="ref")



parser.add_option("-o", "--output",
                help="output path",
                metavar="out", default="Difference", dest="out")



(options, args) = parser.parse_args()


if options.ref=="0" or options.path=="0":
    print "-r and -p options are mandatory to compute difference"


os.system("mkdir -p "+options.out)

outx=open(options.out+"/betxdiff.dat","w")
outy=open(options.out+"/betydiff.dat","w")

print >>outx, "@ REF %s "+options.ref
print >>outy, "@ REF %s "+options.ref
print >>outx, "@ PATH %s "+options.path
print >>outy, "@ PATH %s "+options.path
print >>outx, "$ NAME S   DIFFBETX  STDX"
print >>outy, "$ NAME S   DIFFBETY  STDY"
print >>outx, "* %s   %le %le       %le"
print >>outy, "* %s   %le %le       %le"


def intersect(ListOfFile): 
    '''Pure intersection of all bpm names in all files '''
    if len(ListOfFile) == 0:
        print "Nothing to intersect!!!!"
        sys.exit()
    z = ListOfFile[0].NAME
    for b in ListOfFile:
        z=filter(lambda x: x in z   , b.NAME)
    #SORT by S
    result=[]
    x0=ListOfFile[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))

    result.sort()
    return result






ty1=twiss(options.ref+'/getbetay_free.out')
ty2=twiss(options.path+'/getbetay_free.out')

tx1=twiss(options.ref+'/getbetax_free.out')
tx2=twiss(options.path+'/getbetax_free.out')


bpms=intersect([tx1,tx2])


for ii in range(len(bpms)):
    el=bpms[ii][1]
    i=tx1.indx[el]
    j=tx2.indx[el]
    err=sqrt((tx2.ERRBETX[j]/tx2.BETXMDL[j])**2+(tx1.ERRBETX[i]/tx1.BETXMDL[i])**2)
    std=sqrt((tx2.STDBETX[j]/tx2.BETXMDL[j])**2+(tx1.STDBETX[i]/tx1.BETXMDL[i])**2)
    diff=(tx2.BETX[j]-tx2.BETXMDL[j])/tx2.BETXMDL[j]-(tx1.BETX[i]-tx1.BETXMDL[i])/tx1.BETXMDL[i]
    print >> outx, el, tx1.S[i], diff, err, std



bpms=intersect([ty1,ty2])

for ii in range(len(bpms)):
    el=bpms[ii][1]
    i=ty1.indx[el]
    j=ty2.indx[el]
    err=sqrt((ty2.ERRBETY[j]/ty2.BETYMDL[j])**2+(ty1.ERRBETY[i]/ty1.BETYMDL[i])**2)
    std=sqrt((ty2.STDBETY[j]/ty2.BETYMDL[j])**2+(ty1.STDBETY[i]/ty1.BETYMDL[i])**2)
    diff=(ty2.BETY[j]-ty2.BETYMDL[j])/ty2.BETYMDL[j]-(ty1.BETY[i]-ty1.BETYMDL[i])/ty1.BETYMDL[i]
    print >> outy, el, ty1.S[i], diff, err, std


outx.close()
outy.close()
