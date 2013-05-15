import  sddsdata as sdds

filename="LHCBPM2-Mon_Jun_29_16-26-00_CEST_2009.sdds"

a=sdds.sddsdata(filename, 'big')

print a.data[0].keys()
