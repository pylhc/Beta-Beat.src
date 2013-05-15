#
##
### ! Glenn Vanbavinckhove  !
##
#

filefile=open('mydictionary.py','w')

for i in range(1,96):

    if i<10:
        print >> filefile,"\"M.0"+str(i)+"\":[\"bp0"+str(i)+"\",\"bp0"+str(i)+"\"],"
    else:
        print >> filefile,"\"M."+str(i)+"\":[\"bp"+str(i)+"\",\"bp"+str(i)+"\"],"

filefile.close()
