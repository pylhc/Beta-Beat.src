import sys, os
from optparse import OptionParser
try:
	from metaclass25 import *
except:
	from metaclass import *


parser = OptionParser()
parser.add_option("-I", "--IP",
                help="IP to run: 1,2,5,8. This code will only use Q1s",
                metavar="IP", default="1", dest="IP")
parser.add_option("-k", "--deltak",
                help="Value to change quadrupoles",
                metavar="DELTAK", default="1e-5", dest="dk")
parser.add_option("-m", "--modifiers",
                help="Path to modifiers file with beta* definition, default is recommended",
                metavar="MPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/90m/", dest="mpath")
parser.add_option("-c", "--cpath",
                help="Path to K-mod base source files, default is recommended",
                metavar="CPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/90m/", dest="cpath")

parser.add_option("-R", "--RunMAD",
                help="Run MAD yes [1] or No [0]",
                metavar="RUNMAD", default="1", dest="RunMAD")


(options, args) = parser.parse_args()



dicIPsQuads={"1":["MQXA.1R1","MQXA.1L1"],
	     "2":["MQXA.1R2","MQXA.1L2"],
	     "5":["MQXA.1R5","MQXA.1L5"],
	     "8":["MQXA.1R8","MQXA.1L8"]
	     }


try:
    quads=dicIPsQuads[options.IP]
except:
    print "IP in -I option not defined in dictionary: 1,2,5,8"
    sys.exit()

print "IP quads:",  quads 


cpath=options.cpath
mpath=options.mpath
dk=options.dk

path='./'


if options.RunMAD=="1":
	filename=path+'/var4mad.sh'
	file4nad=open(filename,'w')

	file4nad.write('sed -e \'s/%IP/\''+options.IP+'\'/g\' \\\n')
	file4nad.write('    -e \'s/%MPATH/'+mpath.replace("/","\/")+'/g\' \\\n')
	file4nad.write('    -e \'s/%DELTAK/'+dk+'/g\' \\\n')
	file4nad.write('<'+cpath+'/job.statistics.mask > '+path+'/t.madx \n')
	file4nad.close()
	
	print "Runing MAD, this will take a while"
	os.system("chmod +x "+path+'/var4mad.sh')
	os.system(path+'/var4mad.sh')
	os.system("madx < "+path+'/t.madx')





# Proceeding with fit of MAD output as done in Fitit.py
from MultiVariateRegression import *
from metaclass import *


def writefunc(c,name):
    global path
    
    print "writing function: ",name
    g=open(path+"gnuplot.funcs", "a")
    #Function
    print >>g, name,"(x,y)=", c[0], "+(", c[1],")*x","+(", c[2],")*y","+(", c[3],")*x*y", "+(", c[4],")*x*x", "+(", c[5],")*y*y"
    #Error
    print >>g, name+"err","(x,y,xe,ye)=", "abs(", c[1],"*xe","+(", c[3],")*xe*y","+2*(", c[4],")*x*xe)", "+abs(", c[2],"*ye","+(", c[3],")*x*ye",  "+2*(", c[5],")*y*ye)"
    g.close()

    p=open(path+"pythonfuncs.py", "a")
    print >>p, "def ", name, "(x,y):"
    print >>p, "       return ", c[0], "+(", c[1],")*x","+(", c[2],")*y","+(", c[3],")*x*y", "+(", c[4],")*x*x", "+(", c[5],")*y*y"
    #Error
    print >>p, "def ", name+"err", "(x,y):"
    print >>p, "       return ", "abs(", c[1],"*xe","+(", c[3],")*xe*y","+2*(", c[4],")*x*xe)", "+abs(", c[2],"*ye","+(", c[3],")*x*ye",  "+2*(", c[5],")*y*ye)"
    p.close()



def fitit(ql,qr,b,name):
    global label
    name=label+name
    myData=[]
    for i in range(len(ql)):
       myData.append([b[i],ql[i],qr[i]])
    myEquations = []
    myEquations.append(lambda rawItem, coefIndex: 1  ) # independent
    myEquations.append(lambda rawItem, coefIndex: rawItem[1] ) #ql
    myEquations.append(lambda rawItem, coefIndex: rawItem[2] ) #qr
    myEquations.append(lambda rawItem, coefIndex: rawItem[2]*rawItem[1] )
    myEquations.append(lambda rawItem, coefIndex: rawItem[1]*rawItem[1] )
    myEquations.append(lambda rawItem, coefIndex: rawItem[2]*rawItem[2] )

    res=regression(myData, myEquations)

    print "coeffs:",res
    numCoefficients = len(myEquations)

    term = [0] * numCoefficients

    chi2=0

    for i in range(len(myData)):
        tot=0
        for j in range(numCoefficients):
            tot +=  res[j]*float(myEquations[j](myData[i],j))


        chi2+=(myData[i][0]-tot)*(myData[i][0]-tot)

    print "chi2:", chi2
    writefunc(res,name)
    return res, chi2



b1=twiss('TunesVsBetas.IP'+options.IP+'.b1')
b2=twiss('TunesVsBetas.IP'+options.IP+'.b2')

label='ip'+options.IP
cbxb1=fitit(b1.DQXLB1,b1.DQXRB1,b1.BETXB1 ,"bxb1")
cbyb1=fitit(b1.DQYLB1,b1.DQYRB1,b1.BETYB1 ,"byb1")
cwxb1=fitit(b1.DQXLB1,b1.DQXRB1,b1.WXB1   ,"wxb1")
cwyb1=fitit(b1.DQYLB1,b1.DQYRB1,b1.WYB1   ,"wyb1")
     
     
cbxb2=fitit(b2.DQXLB2,b2.DQXRB2,b2.BETXB2 ,"bxb2")
cbyb2=fitit(b2.DQYLB2,b2.DQYRB2,b2.BETYB2 ,"byb2")
cwxb2=fitit(b2.DQXLB2,b2.DQXRB2,b2.WXB2   ,"wxb2")
cwyb2=fitit(b2.DQYLB2,b2.DQYRB2,b2.WYB2   ,"wyb2")




