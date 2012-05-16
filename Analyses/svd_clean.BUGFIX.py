#!/usr/bin/env python

#---------------------------------
# Use RhicBpmAnalysis to get asc version of sdds file
# Reads asc sdds files and generates a matrix (turns x bpms)
# performs hardware (status 1), peak to peak, and svd localized mode cut
# to remove malfunctioning bpms
# Also removes the noise floor of the data to enhance physical signal
# Ram Calaga, Dec 30, 2003

# usage: ./svd_clean filename
# outputs filename.new

# @ Glenn Vanbavinckhove => made it possible to filter on space
#---------------------------------

#------------------------No Changes below this----------------------------------------------
import sys,os
#try:
#	
#	from Numeric import *
#	from LinearAlgebra import singular_value_decomposition
#	from MLab import mean
#except:
from numpy import *
from numpy.linalg import svd as singular_value_decomposition
from numpy import dot as matrixmultiply

from string import split



from random import gauss, seed
from optparse import OptionParser

#-------------Base Class for parsing sdds file and make nturnsxnbpms matrix with tbt--------



# option parser
parser = OptionParser()
parser.add_option("-t", "--turn",
                help="Turn number to start",
                metavar="startturn", default="40",dest="startturn")
parser.add_option("-p", "--pk-2-pk", 
                help="Peak to peak amplitude cut",
                metavar="peak", default="0.00000001", dest="peak")
parser.add_option("-s", "--sumsquare",
                help="sqroot(sumsquare of the 4 bpm values) should be > 0.9",
                metavar="sun", default="0.925", dest="sum")
parser.add_option("-v", "--sing_val", 
                help="keep svdcut number of singular values in decreasing order",
                metavar="sing", default="100000", dest="sing")
parser.add_option("-d", "--std_dev",
                help="standard deviation for random ditribution",
                metavar="std", default="0.0", dest="std")
parser.add_option("-f", "--file",
                help="file to clean",
                metavar="file", default="./", dest="file")


(options, args) = parser.parse_args()


#----INPUT VARIABLES
turn = int(options.startturn)                # turn number to start
pk_pk_cut = float(options.peak)          # peak to peak amplitude 
sumsquare = float(options.sum)         # sqroot(sumsquare of the 4 bpm values) should be > 0.9
sing_val = int(float(options.sing))         # keep svdcut number of singular values in decreasing order
stdev = float(options.std)             # standard deviation for random ditribution
fle= options.file
#----------------------------------------




class rhicBpmTbt:

    def __init__(self):

        fleb=open(fle.replace("new","bad"),'w');
        print >> fleb, "@  FILE %s ",fle
        print >> fleb, "*  NAME    S    PLANE"
        print >> fleb, "$   %s     %le   %s "
        fleb.close()
        
        self.extractTbt()
        self.makeBpmMatrix()
        self.getValidData()
        #self.randomGauss()
        
        
    def openSddsFile(self):
        
        if len(options.file) > 1:
            f = open(options.file, 'r')
            data = f.readlines()
            f.close()
        else:
            sys.exit('Exit, please enter Filename')
        return data

    
    def extractTbt(self):
        global tx, ty, coordx, coordy, nturns, header,nturns,timestamp,datestamp
        
        data = self.openSddsFile()
        
        tx = []
        ty = []
        coordx = []
        coordy = []
        nturns=timestamp=datestamp=0.0
        
        len_x = len(data)
        #print "length of dataaaaaaaaaa "+str(len_x)
        header = data[0]
        print 'extracting tbt from file\n'
        for i in range(1,len_x):
            s = split(data[i])

            d=s[0]
            if d.startswith('#'):
		if 'NTURNS' in d:
			nturns=s[1]
		elif 'ACQ-TIME' in d:
			timestamp=s[1]
		elif 'ACQ-DATE' in d:
			datestamp=s[1]
                print 'header line '+ s[0]
                 
            else:
                 #print s[1]
                 pl = int(s[0])
                 bpm_name = s[1]
                 if 'null' not in s[2]:
                     s_coord = float(s[2])
                 else:
                     pl=3  # Quick bug fix to avoid appending this BPM
            
                 # first row is scoordinate of each bpm
                 if pl == 0:
                     coordx.append((float(s[2]),s[1]))
                     #coordx.append((str(round(float(s[2]),3)),s[1]))
                     for oo in range(2,len(s)):
                         tx.append(float(s[oo])) 

                 elif pl == 1:
                     coordy.append((float(s[2]),s[1]))
                     #coordy.append((str(round(float(s[2]),3)),s[1]))
                     #coordy.append((s[2],s[1]))
                     for oo in range(2,len(s)):
                         ty.append(float(s[oo]))

                 else:
                     print 'transverse plane unknown for -- ', bpm_name
                
        #num of turns = length(each line) - 2 (plane and bpm name)
        nturns = len(s)-2
        
        #create dictionary keys - s coordinates & values - bpm names
        coordx = dict(coordx)
        coordy = dict(coordy)

        

    def makeBpmMatrix(self):
        global tx, ty, nbpmx, nbpmy, coordx, coordy, nturns

        nbpmx = len(coordx)
        nbpmy = len(coordy)
        
        tx = transpose(reshape(array(tx), (nbpmx,nturns)))
        ty = transpose(reshape(array(ty), (nbpmy,nturns)))
        
        print 'horizontal bpms --',nturns-1, nbpmx
        print 'vertical bpms   --',nturns-1, nbpmy, '\n'
        
        

    def getValidData(self):
        global tx, ty, nbpmx, nbpmy, nturns

        g_indx = []; g_indy = []

        print 'executing hardware status cut\n'
        print turn
        for k in range(nbpmx):
            
            if 99.000 not in tx[turn:,k] or 0.0 not in tx[turn:,k]:
            #if 99.000 not in tx[turn:,k] and 0.0 not in tx[turn:,k]:
                g_indx.append(k)
            #else:
             #   print coordx[str(tx[0,k])]
    
        for l in range(nbpmy):
            #if 99.000 not in tx[turn:,k] and 0.0 not in tx[turn:,k]:
            if 99.000 not in ty[turn:,l] or 0.0 not in ty[turn:,l]:
                g_indy.append(l)
            #else:
             #   print coordy[str(ty[0,l])]
        
        #---remove bpms which are not status bit 1
        tx = take(tx, (g_indx), 1)
        ty = take(ty, (g_indy), 1)
        nbpmx = shape(tx)[1] 
        nbpmy = shape(ty)[1]

        print 'status 1 horizontal bpms --',nturns-1, nbpmx
        print 'status 1 vertical bpms   --',nturns-1, nbpmy, '\n'


    def randomGauss(self):
        global tx, ty, nturns, nbpmx, nbpmy
            
        seed(101010)
        print 'adding gaussian noise to tbt, stdev = ', stdev, '\n'
        for kk in range(nbpmx):
            for mm in range(1,nturns):
                tx[mm,kk] = tx[mm,kk] + gauss(0,stdev)
                
        for nn in range(nbpmy): 
            for pp in range(1,nturns):
                ty[pp,nn] = ty[pp, nn] + gauss(0,stdev)



#-------Derived class to compute SVD of tbt data and clean-----------

class svdBadBpm(rhicBpmTbt):

    def __init__(self):

        rhicBpmTbt.__init__(self)
        self.removeBadBpms('x')
        self.removeBadBpms('y')
        print 'removing noise floor'
        self.svdClean('x')
        self.svdClean('y')
        self.writeFile()
        
        
    def normalizeTbt(self, plane):
        global nturns, tx, ty

        if plane == 'x':
            print 'HORIZONTAL PLANE:'
            b = tx[turn:,:]  #truncate by the first 5 turns
            n_turns = shape(b)[0]

        elif plane == 'y':
            print 'VERTICAL PLANE'
            b = ty[turn:,:]  #truncate by the first 5 turns
            n_turns = shape(b)[0]
        else:
            print "no tbt data acquired"
        print 'Start Turn:',turn
        if shape(b)[0]==0 or shape(b)[1]==0 or len(b)<turn:
            print "No turn by turn data, check start turn \n exiting"; sys.exit() 
        b = (b-mean(b))/sqrt(n_turns)
        return b

            
    def peak_peak(self, plane):
        global nturns, tx, ty
        
        bpm_matrix = self.normalizeTbt(plane)
        n_bpms = shape(bpm_matrix)[1]
        n_turns = shape(bpm_matrix)[0]
        
        #--------------find bpms with signal < 0.4 -- flat signal
        print 'executing peak to peak cut'
        
        pk_pk_ind = []
        for z in range(n_bpms):
            pk_pk = max(bpm_matrix[:,z]) - min(bpm_matrix[:,z])
            if pk_pk > pk_pk_cut/sqrt(n_turns):
                pk_pk_ind.append(z)
        bpm_matrix = take(bpm_matrix, (pk_pk_ind), 1)

        n_bpms = shape(bpm_matrix)[1]
        
        if plane == 'x':
            tx = take(tx, (pk_pk_ind), 1)
            nbpmx = n_bpms
            #print 'horizontal bpms above pk_pk --',nbpmx
            
        if plane == 'y':
            ty = take(ty, (pk_pk_ind), 1)
            nbpmy = n_bpms
            #print '# of vertical bpms above pk_pk --',nbpmx
    
        return bpm_matrix

            
    def svdBpmMatrix(self,plane): #Computes SVD of the bpm_matrix
        global nturns

        bpm_matrix = self.peak_peak(plane)
        n_bpms = shape(bpm_matrix)[1]
        print 'performing SVD'
        #----svd for matrix with bpms >10
        if n_bpms > 10:
            A = singular_value_decomposition(bpm_matrix, full_matrices=0)
        else:
            sys.exit('Exit, # of bpms < 10')
            
        return A
    
        
        
    def removeBadBpms(self, plane):
        global nturns, tx, ty, coordx, coordy, nbpmx, nbpmy

        A = self.svdBpmMatrix(plane)

        #------record spatial vectors to analyze
        V = transpose(array(A[2]))
        V_sort = sort(abs(V), 0)
        V_argsort = argsort(abs(V), 0)
        
        n_bpms = shape(V)[0]

        print 'analyzing spatial vectors', n_bpms, len(V_sort), len(V_sort[0])
        num_rows = []
        for k in range(len(V_sort[0])): # JUST modified, before for k in range(n_bpms)
            for j in range(n_bpms):
                sumsq = sqrt(sum(V_sort[n_bpms-j:n_bpms,k]*V_sort[n_bpms-j:n_bpms,k]))
                if sumsq > sumsquare:
                    num_rows.append(j)
                    break

        bad_indx = []
    
        for r in range(len(num_rows)):
            if num_rows[r] < 5:
                bad_indx.append(r)

        #---record each spatial vector with bad bpms
        V_bad_ind = ravel(take(V_argsort[n_bpms-1:n_bpms,:], (bad_indx), 1))


        #---remove redundancy
        bad_bpms = []
        for qq in V_bad_ind:
            if qq not in bad_bpms:
                bad_bpms.append(qq)


        #---record spatial vectors with physical modes
        good_ind = []
        for ee in range(n_bpms):
            if ee not in bad_bpms:
                good_ind.append(ee)

        if plane == 'x':
            #for bb in range(len(bad_bpms)):
             #   print coordx[str(tx[0,bad_bpms[bb]])]
            fleb=open(fle.replace("new","bad"),'a');#--  append bad bpms to file
#            sys.exit()
            for bb in range(len(bad_bpms)):
                print >> fleb, coordx[tx[0,bad_bpms[bb]]], tx[0,bad_bpms[bb]], "H"
            fleb.close()
            tx = take(tx, (good_ind), 1)
            nbpmx = shape(tx)[1]
            print 'good horizontal bpms --', nbpmx, '\n'
            
        elif plane == 'y':
            #for cc in range(len(bad_bpms)):
             #   print coordy[str(ty[0,bad_bpms[bb]])]
            fleb=open(fle.replace("new","bad"),'a'); #--  append bad bpms to file
            for cc in range(len(bad_bpms)):
                print >> fleb, coordy[ty[0,bad_bpms[cc]]], ty[0,bad_bpms[cc]], "V"
            fleb.close()
            ty = take(ty, (good_ind), 1)
            nbpmy = shape(ty)[1]
            print 'good vertical bpms   --', nbpmy, '\n'

            
    def svdClean(self, plane):
        global nturns, tx, ty

        
        if plane == 'x':
            b = tx[turn:,:]  #truncate by the first 5 turns
            n_turns = shape(b)[0]

        elif plane == 'y':
            b = ty[turn:,:]  #truncate by the first 5 turns
            n_turns = shape(b)[0]
        else:
            print "no tbt data acquired"
        
        b_mean = mean(b)
        b = (b-b_mean)/sqrt(n_turns)
        n_bpms = shape(b)[1]
        #----svd for matrix with bpms >10
        if n_bpms > 10:
            A = singular_value_decomposition(b,full_matrices=0)
            #print "Singular values:",A[1]
        else:
            sys.exit('Exit, # of bpms < 10')
        
        #----SVD cut for noise floor
        if sing_val > n_bpms:
            svdcut = n_bpms
            print 'requested more singular values than available'
            print '# of sing_val used for', plane, '=', n_bpms
        else:
            svdcut = int(sing_val)
            print '# of sing_val used for', plane, '=', svdcut
        #print A[1][0]
	A[1][svdcut:] = 0.
       	#temp=matrixmultiply(identity(len(A[1]))*A[1], A[2])
	temp=matrixmultiply(diag(A[1]), A[2])
	b = matrixmultiply(A[0],temp) ### check
        b = (b *sqrt(n_turns))+b_mean
        #b = b*sqrt(n_turns)
        if plane == 'x':
            tx[turn:,:] = b
        elif plane == 'y':
            ty[turn:,:] = b
        else:
            print "no tbt data to analyze"
        nturns = shape(tx)[0]


    def writeFile(self):
        global tx, ty, coordx, coordy, nbpmx, nbpmy, nturns, header

        f = open(options.file+'.tmp','w') # will change
        f.write(header)
	print >>f,"#NTURNS",nturns
	print >>f,"#ACQ-TIME",timestamp
	print >>f,"#ACQ-DATE",datestamp
        for ff in range(nbpmx):
            f.write('0'+' '+coordx[tx[0,ff]]+'     ')
            #f.write('0'+' '+coordx[str(round(tx[0,ff],3))]+'     ')
            for gg in range(nturns):
                f.write(str(round(tx[gg,ff],5))+'  ')
            f.write('\n')
    
        for hh in range(nbpmy):
            f.write('1'+' '+coordy[ty[0,hh]]+'     ')
            #f.write('1'+' '+coordy[str(round(ty[0,hh],3))]+'     ')
            for jj in range(nturns):
                f.write(str(round(ty[jj,hh],5))+'  ')
            f.write('\n')           
        f.close()
        
        os.system('mv '+options.file+'.tmp'+' '+options.file) # will change
        
svdBadBpm()


