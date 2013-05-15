import sys
if "/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/" not in sys.path: # add internal path for python scripts to current environment (tbach, 2012/05)
    sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")


from numpy import *
from os import system, getpid
import math
import pickle
from metaclass import twiss
from optparse import OptionParser
import multiprocessing
from math import ceil





def shell_command(cmd):
    print 'process id:', getpid()
    ret=system(cmd)
    if ret:
        raise ValueError("COMMAND: %s finished with exit value %i" % (cmd,ret))




def parallel_command(n,period,numberOfCases):
    global path
    iterfile=open(path+'/iter.madx', 'r')
    lines=iterfile.readlines()
    iterfile.close()
    casesperprocess=int(ceil(numberOfCases*1.0/n))
    linesperprocess=casesperprocess*period
    
    for i in range(len(lines)):      #split the iter.madx using in final number of processes 
        if (i%linesperprocess ==0):
            proid=i/linesperprocess+1
            f=open(path+'/iter.'+str(proid)+'.madx', 'w')
        print  >>f, lines[i][:-1]
        if i==len(lines)-1 or (i%linesperprocess == linesperprocess-1):
            #print "closing ", i
            f.close()
            
    newn=proid  # Redefine number of process in case it ended up being smaller
    
    # Prepare copies of the job.iterate.madx and all the shell commands
    commandlist=[]
    for i in range(1,newn+1):
        cmd='sed \'s/iter.madx/iter.'+str(i)+'.madx/g\' '+path+'/job.iterate.madx > '+path+'/job.iterate.'+str(i)+'.madx'
        shell_command(cmd)
        cmd='/afs/cern.ch/group/si/slap/bin/madx < '+path+'/job.iterate.'+str(i)+'.madx > ttt'+str(i)
        commandlist.append(cmd)


    # Submit the shell commands in parallel
    pool.map(shell_command, commandlist)
    print commandlist



def loadtwiss_beta(varandpath):
        var, path=varandpath
        print "Reading twiss."+var
        x=twiss(path+"/twiss."+var)
        system('rm '+path+'/twiss.'+var)
        return var, x

def loadtwiss_coup(varandpath):
        var, path = varandpath
        print "Reading twiss."+var
        x=twiss(path+"/twiss."+var)
        x.Cmatrix()
        system('rm '+path+'/twiss.'+var)
        return var, x

def loadtwiss_chrom_coup(varandpathanddpp):

    var, path, dpp = varandpathanddpp
    print  "Reading twiss.dp+."+var
    xp=twiss(path+"/twiss.dp+."+var)
    xp.Cmatrix()
    print  "Reading twiss.dp-."+var
    xm=twiss(path+"/twiss.dp-."+var)
    xm.Cmatrix()
    # Initializing and Calculating chromatic coupling for every BPM
    xp.Cf1001r=[]
    xp.Cf1001i=[]
    xp.Cf1010r=[]
    xp.Cf1010i=[]
    for j in range(len(xp.NAME)):
        
        vvv=(xp.F1001R[j]-xm.F1001R[j])/(2*dpp)
        xp.Cf1001r.append(vvv)
        
        vvv=(xp.F1001I[j]-xm.F1001I[j])/(2*dpp)
        xp.Cf1001i.append(vvv)

        vvv=(xp.F1001R[j]-xm.F1001R[j])/(2*dpp)
        xp.Cf1010r.append(vvv)
        
        vvv=(xp.F1010I[j]-xm.F1010I[j])/(2*dpp)
        xp.Cf1010i.append(vvv)

    #FullResponse[var]=xp
    system('rm '+path+'/twiss.dp+.'+var)
    system('rm '+path+'/twiss.dp-.'+var)
    return var, xp






    
if __name__ == '__main__':
    
    numberofCPUs=multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=numberofCPUs)


    ##### optionparser
    parser = OptionParser()
    parser.add_option("-a", "--accel",
                      help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                      metavar="ACCEL", default="LHCB1",dest="accel")
    parser.add_option("-p", "--path",
                      help="path to save",
                metavar="path", default="./",dest="path")
    parser.add_option("-c", "--core",
                      help="core files",
                      metavar="core", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/",dest="core")
    parser.add_option("-k", "--deltak",
                      help="delta k to be applied to quads for sensitivity matrix",
                      metavar="core", default="0.00002",dest="k")


    (options, args) = parser.parse_args()

    # paths
    corepath=options.core
    accel=options.accel
    path=options.path


    #
    # Chromatic coupling
    #

    FullResponse={}   #Initialize FullResponse
    execfile(corepath+"/"+accel+'/AllLists_chromcouple.py')
    exec('variables=kss()')           #Define variables
    delta1=zeros(len(variables))*1.0   #Zero^th of the variables
    incr=ones(len(variables))*0.05    #increment of variables
    dpp=0.0001
    FullResponse['incr']=incr           #Store this info for future use
    FullResponse['delta1']=delta1


######## loop over normal variables
    f=open(path+'/iter.madx','w')
    for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        var=variables[i]
        print >>f, var,"=", var, "+(",delta[i],");"
        print >>f, "twiss, deltap= "+str(dpp)+",file=\""+path+"/twiss.dp+."+var+"\";"
        print >>f, "twiss, deltap=-"+str(dpp)+",file=\""+path+"/twiss.dp-."+var+"\";"
        print >>f, var,"=", var, "-(",delta[i],");"


    print >>f, "twiss, deltap= "+str(dpp)+",file=\""+path+"/twiss.dp+.0\";"
    print >>f, "twiss, deltap=-"+str(dpp)+",file=\""+path+"/twiss.dp-.0\";"
    f.close()
    print "Runing MADX"
    #shell_command('/afs/cern.ch/group/si/slap/bin/madx < '+path+'/job.iterate.madx')
    parallel_command(n=5,period=4,numberOfCases=len(delta1)+1) # period=4 since there are 4 lines in iter.madx per case, numberofcases has +1 since there is the 0 case

    
    
    varsforloop=variables+['0']
    newvarsforloop=[]
    for xxx in varsforloop:
        newvarsforloop.append([xxx,path, dpp])
    a=pool.map(loadtwiss_chrom_coup, newvarsforloop)
    for var, xxx in a:
        FullResponse[var]=xxx
    
    pickle.dump(FullResponse,open(path+'/FullResponse_chromcouple','w'),-1)

    

    #
    # Coupling
    #

    FullResponse={}   #Initialize FullResponse
    execfile(corepath+"/"+accel+'/AllLists_couple.py')
    exec('variables=Qs()')           #Define variables
    delta1=zeros(len(variables))*1.0   #Zero^th of the variables
    incr=ones(len(variables))*0.0001    #increment of variables


    FullResponse['incr']=incr           #Store this info for future use
    FullResponse['delta1']=delta1       #"     "     "

    ######## loop over normal variables
    f=open(path+'/iter.madx','w')
    for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        var=variables[i]
        print >>f, var,"=", var, "+(",delta[i],");"
        print >>f, "twiss, file=\""+path+"/twiss."+var+"\";"
        print >>f, var,"=", var, "-(",delta[i],");"

    print >>f, "twiss, file=\""+path+"/twiss.0\";"
    f.close()
    #Sending the mad jobs in parallel
    parallel_command(n=7,period=3,numberOfCases=len(delta1)+1)
    #Loading the twiss files into fullresp in parallel
    varsforloop=variables+['0']
    newvarsforloop=[]
    for xxx in varsforloop:
        newvarsforloop.append([xxx,path])
    a=pool.map(loadtwiss_coup, newvarsforloop)
    for var, xxx in a:
        FullResponse[var]=xxx
    pickle.dump(FullResponse,open(path+'/FullResponse_couple','w'),-1)

    

    #
    #
    # Beta
    #
    #
    FullResponse={}   #Initialize FullResponse
    execfile(corepath+"/"+accel+'/AllLists.py')
    exec('variables=Q()')           #Define variables
    delta1=zeros(len(variables))*1.0   #Zero^th of the variables
    #incr=ones(len(variables))*0.00005    #increment of variables    #### when squeeze low twiss fails because of to big delta
    incr=ones(len(variables))*float(options.k)


    FullResponse['incr']=incr           #Store this info for future use
    FullResponse['delta1']=delta1       #"     "     "

    ######## loop over normal variables
    f=open(path+'/iter.madx','w')
    for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        var=variables[i]
        print >>f, var,"=", var, "+(",delta[i],");"
        print >>f, "twiss, file=\""+path+"/twiss."+var+"\";"
        print >>f, var,"=", var, "-(",delta[i],");"

    print >>f, "twiss,file=\""+path+"/twiss.0\";"
    f.close()

    parallel_command(n=numberofCPUs,period=3,numberOfCases=len(delta1)+1)
    #
    varsforloop=variables+['0']
    newvarsforloop=[]
    for xxx in varsforloop:
        newvarsforloop.append([xxx,path])

    a=pool.map(loadtwiss_beta, newvarsforloop)
    for var, xxx in a:
        FullResponse[var]=xxx

    pickle.dump(FullResponse,open(path+'/FullResponse','w'),-1)
