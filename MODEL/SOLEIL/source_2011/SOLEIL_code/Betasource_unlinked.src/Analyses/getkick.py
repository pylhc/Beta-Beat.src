
from  rhicdata import rhicdata
from Numeric import *
import sys

datafile=sys.argv[1]

d=rhicdata(datafile)

#kickamp=400
# modified 2007  wed sep 19, Glenn Vanbavinckhove
kickamp=float(sys.argv[2])
# end modification



#(begin)  added 02/10/2007 Glenn Vanbavinckhove

f=open(datafile+'.kick','w')

print >>f, "* NAME", "TURN", "KICK", "RMS", "PLANE"
print >>f, "$ %s      %le      %le      %le     %le"


# (end)

hturn={}
vturn={}

for j in range(len(d.H)):
	found=0
	a=d.H[j].data
	for i in range(len(a)-1):
	       
		if abs(a[i+1]-a[i])>kickamp:
			#print '"'+d.H[j].name+'"', i, abs(a[i+1]-a[i]), sqrt(average(a**2)-average(a)**2), "0"
			print >>f, '"'+d.H[j].name+'"', i, abs(a[i+1]-a[i]), sqrt(average(a**2)-average(a)**2), "0"


	 		found=1
			try:
				tv=hturn[i]
			except:
				tv=0
				hturn[i]=0
			hturn[i]=tv+1
			break
	if found == 0:
		
		#print '"'+d.H[j].name+'"', -1, -1, sqrt(average(a**2)-average(a)**2) , "0"
		print >>f, '"'+d.H[j].name+'"', -1, -1, sqrt(average(a**2)-average(a)**2) , "0" 

for j in range(len(d.V)):
	found=0
	a=d.V[j].data
	for i in range(len(a)-1):
		if abs(a[i+1]-a[i])>kickamp:
			print >>f, '"'+d.V[j].name+'"', i, abs(a[i+1]-a[i]), sqrt(average(a**2)-average(a)**2), "1"
			found=1
			try:
				tv=vturn[i]
			except:
				tv=0
				vturn[i]=0
			vturn[i]=tv+1	
			break
	if found == 0:
		print >>f,'"'+d.V[j].name+'"', -1, -1, sqrt(average(a**2)-average(a)**2), "1"
			



if len(hturn) > 0:
  mv=max(hturn.values())
  for ke in sort(hturn.keys()):
	if mv==hturn[ke]: 
		keh=ke
else:
  keh=10000

if len(hturn) > 0:
  mv=max(vturn.values())
  for ke in sort(vturn.keys()):
        if mv==vturn[ke]:                 
		kev=ke 
else:
  kev=10000


print >>f,"@ HTURNMOST ", keh
print >>f,"@ VTURNMOST ", kev

print >>f,"@ KICKAMP ", kickamp

f.close()



