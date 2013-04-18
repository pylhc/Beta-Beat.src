from metaclass import twiss

twissno=twiss('twiss_elements.dat')
twissac=twiss('twiss_elements_ac.dat')

filefile=open("acdipole.tfs","w")


betaxIP1=twissno.BETX[twissno.indx['IP1']]-twissac.BETX[twissac.indx ['IP1']]
betayIP1=twissno.BETY[twissno.indx['IP1']]-twissac.BETY[twissac.indx['IP1']]

betaxIP2=twissno.BETX[twissno.indx['IP2']]-twissac.BETX[twissac.indx['IP2']]
betayIP2=twissno.BETY[twissno.indx['IP2']]-twissac.BETY[twissac.indx['IP2']]

betaxIP5=twissno.BETX[twissno.indx['IP5']]-twissac.BETX[twissac.indx['IP5']]
betayIP5=twissno.BETY[twissno.indx['IP5']]-twissac.BETY[twissac.indx['IP5']]

betaxIP8=twissno.BETX[twissno.indx['IP8']]-twissac.BETX[twissac.indx['IP8']]
betayIP8=twissno.BETY[twissno.indx['IP8']]-twissac.BETY[twissac.indx['IP8']]

print >> filefile,"* NAME DBETX DBETY"
print >> filefile,"$ %s %le %le"

print >> filefile,"IP1",betaxIP1,betayIP1
print >> filefile,"IP2",betaxIP2,betayIP2
print >> filefile,"IP5",betaxIP5,betayIP5
print >> filefile,"IP8",betaxIP8,betayIP8

filefile.close()
