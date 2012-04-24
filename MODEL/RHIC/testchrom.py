
from metaclass import *

a=twiss("twiss.c.dat")
a.chrombeat()
for i in range(len(a.S)):
  print a.NAME[i], a.S[i], a.dbx[i]
