#  Date : 29/02/2009
#  Created by Glenn Vanbavinckhove
#  Purpose : creates zero changeparameters

# importing libraries



# initializing some stuff
  #=> variables to the list quads
variablesQ = ['qp1','qp2','qp3','qp4','qp5','qp6','qp7','qp8','qp44','qp55']

  #=> variables to the list sex
variablesS = ["sx1","sx2","sx3","sx4","sx5","sx6","sx7","sx8","sx9","sx10"]


sexChange=open('changeparametersS','w')
quadChange=open('changeparametersQ','w')

for i in variablesQ:

    quadChange.write(str(i) +' ='+str(i)+'+ (0.0);\n')
    


for i in variablesS:

    sexChange.write(str(i) +' ='+str(i)+'+ (0.0);\n')



quadChange.close()
sexChange.close()
