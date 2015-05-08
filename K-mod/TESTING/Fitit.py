from MultiVariateRegression import *
from metaclass import *
















def writefunc(c,name):
    
    print "writing function: ",name
    g=open("gnuplot.funcs", "a")
    print >>g, name,"(x,y)=", c[0], "+(", c[1],")*x","+(", c[2],")*y",  "+(", c[3],")*x*y", "+(", c[4],")*x*x", "+(", c[5],")*y*y"
    g.close()

    p=open("pythonfuncs.py", "a")
    print >>p, "def ", name, "(x,y):"
    print >>p, "       return ", c[0], "+(", c[1],")*x","+(", c[2],")*y",  "+(", c[3],")*x*y", "+(", c[4],")*x*x", "+(", c[5],")*y*y"    
    p.close()




def fitit(ql,qr,b,name):
    

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






b1=twiss('TunesVsBetas.b1')
b2=twiss('TunesVsBetas.b2')


cbxb1=fitit(b1.DQXLB1,b1.DQXRB1,b1.BETXB1 ,"bxb1")
cbyb1=fitit(b1.DQYLB1,b1.DQYRB1,b1.BETYB1 ,"byb1")
cwxb1=fitit(b1.DQXLB1,b1.DQXRB1,b1.WXB1   ,"wxb1")
cwyb1=fitit(b1.DQYLB1,b1.DQYRB1,b1.WYB1   ,"wyb1")
     
     
cbxb2=fitit(b2.DQXLB2,b2.DQXRB2,b2.BETXB2 ,"bxb2")
cbyb2=fitit(b2.DQYLB2,b2.DQYRB2,b2.BETYB2 ,"byb2")
cwxb2=fitit(b2.DQXLB2,b2.DQXRB2,b2.WXB2   ,"wxb2")
cwyb2=fitit(b2.DQYLB2,b2.DQYRB2,b2.WYB2   ,"wyb2")



