from metaclass import twiss

beam1=twiss('./twiss_elements_b1.dat')

beam2=twiss('./twiss_elements_b2.dat')

#beam1

filefile=open('instrumentb1','w')


####### beam1
names=beam1.NAME
stringb1="empty"

datab1=open("listb1","w")

for name in names:

	key=beam1.KEYWORD[beam1.indx[name]]


	if ("INSTRUMENT"  in key) or ("COLLIMATOR" in key)  or ("RCOLLIMATOR" in key) or ("TKICKER" in key) or ("RBEND" in key):

		stringb1=stringb1+","+name
		

	if ("IP" in name):
		stringb1=stringb1+","+name

stringb1=stringb1.replace("empty,","")

print >> datab1,stringb1

datab1.close()


####### beam2
names=beam2.NAME
stringb2="empty"

datab2=open("listb2","w")

for name in names:

	key=beam2.KEYWORD[beam2.indx[name]]


	if ("INSTRUMENT"  in key) or ("COLLIMATOR" in key)  or ("RCOLLIMATOR" in key) or ("TKICKER" in key) or ("RBEND" in key):

		stringb2=stringb2+","+name
		

	if ("IP" in name):
		stringb2=stringb2+","+name

stringb2=stringb2.replace("empty,","")

print >> datab2,stringb2

datab2.close()
