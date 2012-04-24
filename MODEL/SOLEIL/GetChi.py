from metaclass import twiss

x=twiss("twiss.dat")

listofbpms=filter(lambda y: 'BP' in y.upper(), x.NAME)

x.chiterms(listofbpms)

f = open('chi3000.dat','w')
for i in range(len(x.chiBPMs)):
    print >>f, x.chiS[i][0], x.chiS[i][1], x.chiS[i][2], abs(x.chi[i]), x.chi[i].real, x.chi[i].imag
    


f.close()
