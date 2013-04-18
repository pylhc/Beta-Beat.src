#  Date : 28/02/2009
#  Created by Glenn Vanbavinckhove
#  Purpose : creates bump for Fullresponse matrix

# importing libraries



# initializing some stuff
  #=> variables to the list quads
variablesQ = ['qp1','qp2','qp3','qp4','qp5','qp6','qp7','qp8','qp44','qp55']

  #=> variables to the list sex
variablesS = ["sx1","sx2","sx3","sx4","sx5","sx6","sx7","sx8","sx9","sx10"]
  
  #=> to write to bump information
quadFile=open('BumpsQuads','w')
sexFile=open('BumpsSex','w')
  #=> to open the changeparameters files
sexChange=open('changeparametersS','r')
quadChange=open('changeparametersQ','r')



# main method
  # => for quads
quadLines=quadChange.readlines()

for i in variablesQ:

    line=quadLines[1].split()

    quadFile.write(i+'->KICK:= '+i+'\n') 

  # => for sextupoles
sexLines=sexChange.readlines()

for i in variablesS:

    line=sexLines[1].split()

    sexFile.write(i+'->KICK:= '+i+'\n') 



#closing the files
quadFile.close()
sexFile.close()
