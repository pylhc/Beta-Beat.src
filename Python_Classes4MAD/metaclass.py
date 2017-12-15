# version 1 unknown/unknown
# version 2 20120905 (tbach):
# - major refactoring
# - cleared indentation
# - moved help section in comments to top and docstring
# - removed all exec statements (replaced with setattr) (still not good)
# - added class attributes for all unknown members
# - declared constants for PI, E, I
# - cleaned up imports
# - removed all trailing ";"
# - tested it for getllm/lhc, produces exactly same results
# (getllm called only Cmatrix and chiterms functions))
# 

"""
Read the twiss class from the twiss file
x=twiss('twiss')
use it as:
print x.Q1, x.Q2, x.BETX[0]


run beaMatrix for example:
x.beatMatrix()
print x.RM[0]


BETA-BEAT CORRECTION
first compute the response matrix by:
x.beatMatrix()
Define targetbeat as an array containing the desired changed in Dbeta/beta (x,y)
targetbeat=Dbeta/beta
dkl gives the required integrated strengths by:
dkl=matrixmultiply(generalized_inverse(x.RM,0.003),targetbeat)


Want to explore the singular values?:
svd=singular_value_decomposition(x.RM)


Computing SEXTUPOLAR RESONANCE TERMS:
x.fterms()
The fterms are arrays evaluated at all the elements:
print x.f3000 , x.f2100  , x.f1020, x.f1002


COUPLING
Compute the Cmatrix, gamma, f1001 and f1010 from the Twiss file at all elements
x.Cmatrix()
print x.C[0]   (four components of C at the first elements)
print x.f1001
...
"""

import numpy
from numpy import dot as matrixmultiply
from numpy.linalg import inv as inverse
from numpy.linalg import det as determinant
from math import factorial


import sys

I = complex(0, 1)
E = numpy.e
PI = numpy.pi

class twiss:
    """Twiss parameters from madx output (with free choice of select items)"""

    def forknames(self, dictionary):
        NAME = getattr(self, "NAME")
        for n in dictionary:
            if n in NAME:
                for m in dictionary[n]:
                    self.indx[m] = self.indx[n]
                    self.indx[m.upper()] = self.indx[n]
            else:
                print "skipped value from dictionary because not in NAME. value: ", n

    def __init__(self, filename, dictionary=None):
        if dictionary is None:
            dictionary = {}
            
        self.filename = filename # Added to see which file it is during debugging (vimaier)
        self.__has_parsed_a_table_row = False 
        self.indx = {}
        self.keys = []
        alllabels = []

        if filename.endswith(".gz"):
            import gzip
            f = gzip.open(filename, 'rb')
        else:
            f = open(filename, 'r')
        for line in f:
            is_line_parsed = False # Check if line was parsed otherwise print info (vimaier)
            
            if line.startswith("#"): # comment line
                continue
            
            if ("@ " not in line and "@" in line):
                line = line.replace("@" , "@ ")
            split_line = line.split()
            
            if ("@ " in line and "%" in line and "s" not in split_line[2]):
            # Float-Descriptor-line
                label = split_line[1]
                try:
                    setattr(self, label, float(split_line[3].replace("\"", "")))
                except:
                    print "Problem parsing:", line, 
                    print "Going to be parsed as string"
                    try:
                        setattr(self, label, split_line[3].replace("\"", ""))
                    except:
                        print "Problem persists, let's ignore it!"
                is_line_parsed = True
            elif ("@ " in line and "s" in split_line[2]):
            # String-Descriptor-line
                label = split_line[1].replace(":", "")
                setattr(self, label, " ".join(split_line[3:]).replace("\"", ""))
                is_line_parsed = True
                
            if ("* " in line or "*\t" in line):
            # Columns-names-line
                alllabels = split_line
                for alllabels_item in alllabels[1:]:
                    setattr(self, alllabels_item, [])
                    self.keys.append(alllabels_item)
                is_line_parsed = True

            if ("$ " in line or "$\t" in line):
            # Columns-datatypes-line
                alltypes = split_line
                is_line_parsed = True

            if ("@" not in line and "*" not in line and "$" not in line and "#" not in line):
            # Table-entry-line
                values = split_line
                for j in range(0,len(values)):
                    if ("%hd" in alltypes[j + 1] or "%d" in alltypes[j + 1] ):
                        getattr(self, alllabels[j + 1]).append(int(values[j]))
                    if ("%le" in alltypes[j + 1]):
                        getattr(self, alllabels[j + 1]).append(float(values[j]))
                    if ("%im" in alltypes[j + 1]):
                        getattr(self, alllabels[j + 1]).append(complex(values[j]))
                    if ("s" in alltypes[j+1]):
                        getattr(self, alllabels[j + 1]).append(values[j].replace("\"", ""))
                        if "NAME" == alllabels[j + 1]:
                            NAME = getattr(self, "NAME")
                            self.indx[values[j].replace("\"", "")] = len(NAME) - 1
                            self.indx[values[j].replace("\"", "").upper()] = len(NAME) - 1
                            self.indx[values[j].replace("\"", "").lower()] = len(NAME) - 1
                self.__has_parsed_a_table_row = True
                is_line_parsed = True
            
            if not is_line_parsed:
                print >> sys.stderr,"Did not parse line ("," ".join(split_line),") in ",filename

        f.close()
        try:
            alltypes
        except:
            print >> sys.stderr, "From Metaclass: Bad format or empty file ", filename
            raise ValueError


        for j in range(1, len(alllabels)):
            if (("%le" in alltypes[j]) | ("%hd" in alltypes[j])):
                setattr(self, alllabels[j], numpy.array(getattr(self, alllabels[j])))

        if len(dictionary) > 0:
            self.forknames(dictionary)
            
    def has_bpm_data(self):
        return self.__has_parsed_a_table_row
    
    def has_no_bpm_data(self):
        return not self.has_bpm_data()

    def chrombeat(self):
        '''
         Add dbx/dby to the twiss table
        '''
        self.dbx = []
        self.dby = []
        S = getattr(self, "S")
        WX = getattr(self, "WX")
        WY = getattr(self, "WY")
        PHIX = getattr(self, "PHIX")
        PHIY = getattr(self, "PHIY")
        for i in range(0, len(S)):
            ax = WX[i] * numpy.cos(PHIX[i] * 2 * PI)
            ay = WY[i] * numpy.cos(PHIY[i] * 2 * PI)
            self.dbx.append(ax)
            self.dby.append(ay)

    def fterms(self):
        '''
         Add f terms to the twiss table
        '''

        self.f3000 = []
        self.f2100 = []
        self.f1020 = []
        self.f1002 = []
        self.f20001 = []
        self.f1011 = []
        self.f2000 = []
        self.f4000 = []
        self.f3100 = []
        self.f2020 = []
        self.f1102 = []
        self.f2011 = []
        self.f2002 = []
        self.f0031 = []
        self.f0040 = []

        S = getattr(self, "S")
        MUX = getattr(self, "MUX")
        MUY = getattr(self, "MUY")
        Q1 = getattr(self, "Q1")
        Q2 = getattr(self, "Q2")
        K1L = getattr(self, "K1L")
        K2L = getattr(self, "K2L")
        K3L = getattr(self, "K3L")
        BETX = getattr(self, "BETX")
        BETY = getattr(self, "BETY")
        DX = getattr(self, "DX")
        for i in range(0, len(S)):
            phix = MUX - MUX[i]
            phiy = MUY - MUY[i]
            for j in range(0, i):
                phix[j] += Q1
                phiy[j] += Q2
            dumm = -sum(K2L * BETX ** 1.5 * E ** (3 * I * 2 * PI * phix)) / 48.
            self.f3000.append(dumm / (1. - E ** (3 * I * 2 * PI * Q1)))
            dumm = -sum(K2L * BETX ** 1.5 * E ** (I * 2 * PI * phix)) / 16.
            self.f2100.append(dumm / (1. - E ** (I * 2 * PI * Q1)))
            dumm = sum(K2L * BETX ** 0.5 * BETY * E ** (I * 2 * PI * (phix + 2 * phiy))) / 8.
            self.f1020.append(dumm / (1. - E ** (I * 2 * PI * (Q1 + 2 * Q2))))
            dumm = sum(K2L * BETX ** 0.5 * BETY * E ** (I * 2 * PI * (phix - 2 * phiy))) / 8.
            self.f1002.append(dumm / (1. - E ** (I * 2 * PI * (Q1 - 2 * Q2))))
            dumm = sum((K1L - 2 * K2L * DX) * BETX * E ** (2 * I * 2 * PI * phix)) / 8.
            self.f20001.append(dumm / (1. - E ** (2 * I * 2 * PI * Q1)))
            dumm = sum(K2L * BETX ** 0.5 * BETY * E ** (I * 2 * PI * (phix))) / 4.
            self.f1011.append(dumm / (1. - E ** (I * 2 * PI * Q1)))
            dumm = -sum(K1L * BETX ** 1 * E ** (2 * I * 2 * PI * phix)) / 32.
            self.f2000.append(dumm / (1. - E ** (2 * I * 2 * PI * Q1)))
            dumm = -sum(K3L * BETX ** 2 * E ** (4 * I * 2 * PI * (phix))) / 384.
            self.f4000.append(dumm / (1. - E ** (4 * I * 2 * PI * Q1)))
            dumm = -sum(K3L * BETX ** 2 * E ** (2 * I * 2 * PI * phix)) / 96
            self.f3100.append(dumm / (1. - E ** (2 * I * 2 * PI * Q1)))
            dumm = sum(K3L * BETX * BETY * E ** (2 * I * 2 * PI * (phix + phiy))) / 64
            self.f2020.append(dumm / (1. - E ** (2 * I * 2 * PI * (Q1 + Q2))))
            dumm = sum(K3L * BETX * BETY * E ** (-2 * I * 2 * PI * phiy)) / 32
            self.f1102.append(dumm / (1. - E ** (-2 * I * 2 * PI * Q2)))
            dumm = sum(K3L * BETX * BETY * E ** (2 * I * 2 * PI * phix)) / 32
            self.f2011.append(dumm / (1. - E ** (2 * I * 2 * PI * Q1)))            
            dumm = sum(K3L * BETX * BETY * E ** (2 * I * 2 * PI * (phix - phiy))) / 64
            self.f2002.append(dumm / (1. - E ** (2 * I * 2 * PI * (Q1 - Q2))))
            dumm = -sum(K3L * BETY ** 2 * E ** (2 * I * 2 * PI * phiy)) / 96
            self.f0031.append(dumm / (1. - E ** (2 * I * 2 * PI * Q2)))
            dumm = -sum(K3L * BETY ** 2 * E ** (4 * I * 2 * PI * (phiy))) / 384.
            self.f0040.append(dumm / (1. - E ** (4 * I * 2 * PI * Q2)))

        self.f3000 = numpy.array(self.f3000)
        self.f2100 = numpy.array(self.f2100)
        self.f1020 = numpy.array(self.f1020)
        self.f1002 = numpy.array(self.f1002)
        self.f1011 = numpy.array(self.f1011)
        self.f4000 = numpy.array(self.f4000)
        self.f3100 = numpy.array(self.f3100)
        self.f2020 = numpy.array(self.f2020)
        self.f1102 = numpy.array(self.f1102)
        self.f2011 = numpy.array(self.f2011)
        self.f2002 = numpy.array(self.f2002)
        self.f0031 = numpy.array(self.f0031)
        self.f0040 = numpy.array(self.f0040)

        self.f0120 = numpy.conjugate(self.f1002)
        self.f0111 = numpy.conjugate(self.f1011)
        self.f1200 = numpy.conjugate(self.f2100)

        self.fRS3 = 3 * self.f3000 - self.f2100
        self.fRS2 = self.f1020 - self.f0120
        self.fRS1 = 2 * self.f1020 - self.f1011
        self.fRS1 = 2 * self.f0120 - self.f0111


    def fterms_generic(self,fterms_list=['3000','2100','1020','1002','1011','2000','4000','0120','0111','1200'], conjugate_list=['0120','0111','1200']):
   
        for fterm in fterms_list:
            setattr(self, 'f'+fterm, [])

        self.f20001 = []

        S = getattr(self, "S")
        MUX = getattr(self, "MUX")
        MUY = getattr(self, "MUY")
        Q1 = getattr(self, "Q1")
        Q2 = getattr(self, "Q2")
        K1L = getattr(self, "K1L") 
        K2L = getattr(self, "K2L")
        K3L = getattr(self, "K3L")
        K1SL = getattr(self, "K1SL") 
        K2SL = getattr(self, "K2SL")
        K3SL = getattr(self, "K3SL")
        BETX = getattr(self, "BETX")
        BETY = getattr(self, "BETY")
        DX = getattr(self, "DX")
        KL = [K1L, K2L, K3L]
        KSL = [K1SL, K2SL, K3SL]
        
        print 'sum of K2L:  ' , sum(abs(K2L))
        print 'sum of K3L:  ' , sum(abs(K3L))
        
        for i in range(0, len(S)):
            phix = MUX - MUX[i]
            phiy = MUY - MUY[i]
            for j in range(0, i):
                phix[j] += Q1
                phiy[j] += Q2
            
            for fterm in fterms_list:
                [j,k,l,m] = [int(i) for i in fterm]
                n = j+k+l+m
                factor = -(I**(l+m)) / ( factorial(j)*factorial(k)*factorial(l)*factorial(m) * 2**n )
                if (l+m)%2==0:
                    factor = factor.real
                    dumm = sum(KL[n-2]  * BETX**((j+k)/2.) * BETY**((l+m)/2.) * E**(I*2*PI*((j-k)*phix + (l-m)*phiy)) ) * factor
                    getattr(self, 'f'+fterm).append(dumm / (1. - E**(I*2*PI*((j-k)*Q1 + (l-m)*Q2))) )
                elif (l+m)%2==1:
                    factor = factor.imag
                    dumm = sum(KSL[n-2] * BETX**((j+k)/2.) * BETY**((l+m)/2.) * E**(I*2*PI*((j-k)*phix + (l-m)*phiy)) ) * factor     
                    getattr(self, 'f'+fterm).append(dumm / (1. - E**(I*2*PI*((j-k)*Q1 + (l-m)*Q2))) )
        

        #Legacy from old non-Generic version for backwards compatibility
        dumm = sum((K1L - 2 * K2L * DX) * BETX * E ** (2 * I * 2 * PI * phix)) / 8.
        self.f20001.append(dumm / (1. - E ** (2 * I * 2 * PI * Q1)))
        self.f20001 = numpy.array(self.f20001)

        for fterm in fterms_list:
            setattr(self, 'f'+fterm, numpy.array(getattr(self, 'f'+fterm)))
       
        self.fRS3 = 3 * self.f3000 - self.f2100
        self.fRS2 = self.f1020 - self.f0120
        self.fRS1 = 2 * self.f1020 - self.f1011
        self.fRS1 = 2 * self.f0120 - self.f0111
        

    def fterms_generic_test(self,fterms_list=['3000','2100','1020','1002','1011','2000','4000','0120','0111','1200'], conjugate_list=['0120','0111','1200']):
   
        for fterm in fterms_list:
            setattr(self, 'f'+fterm, [])

        self.f20001 = []

        S = getattr(self, "S")
        MUX = getattr(self, "MUX")
        MUY = getattr(self, "MUY")
        Q1 = getattr(self, "Q1")
        Q2 = getattr(self, "Q2")
        K1L = getattr(self, "K1L") 
        K2L = getattr(self, "K2L")
        K3L = getattr(self, "K3L")
        BETX = getattr(self, "BETX")
        BETY = getattr(self, "BETY")
        DX = getattr(self, "DX")
        KL = [K1L, K2L, K3L]
        
        print 'sum of K2L:  ' , sum(abs(K2L))
        
        for i in range(len(S)):
            phix = MUX - MUX[i]
            phiy = MUY - MUY[i]
            
            tune_mask = numpy.concatenate([numpy.ones(i), numpy.zeros(len(S)-i)])
            phix = phix + tune_mask*Q1 
            phiy = phiy + tune_mask*Q2 
            
            for fterm in fterms_list:
                [j,k,l,m] = [int(i) for i in fterm]
                n = j+k+l+m
                factor = -(I ** (l+m) ) / (factorial(j) * factorial(k) * factorial(l) * factorial(m) * 2 ** n )
                if (l+m)%2==0:
                    factor = factor.real
                    dumm = sum(KL[n-2] * BETX ** ((j+k)/2.) * BETY ** ((l+m)/2.) * E ** (I * 2 * PI * ( (j-k) * phix + (l-m) * phiy ) ) ) * factor
                    getattr(self, 'f'+fterm).append(dumm / (1. - E ** (I * 2 * PI * ((j-k) * Q1 + (l-m) * Q2)) ) )
                elif (l+m)%2==1:
                    factor = factor.imag
                    dumm = sum(KL[n-2] * BETX ** ((j+k)/2.) * BETY ** ((l+m)/2.) * E ** (I * 2 * PI * ((j-k) * phix + (l-m) * phiy)) ) * factor     
                    getattr(self, 'f'+fterm).append(dumm / (1. - E ** (I * 2 * PI * ((j-k) * Q1 + (l-m) * Q2)) ) )
                    
        dumm = sum((K1L - 2 * K2L * DX) * BETX * E ** (2 * I * 2 * PI * phix)) / 8.
        self.f20001.append(dumm / (1. - E ** (2 * I * 2 * PI * Q1)))
        self.f20001 = numpy.array(self.f20001)

        for fterm in fterms_list:
            setattr(self, 'f'+fterm, numpy.array(getattr(self, 'f'+fterm)))
       
        self.fRS3 = 3 * self.f3000 - self.f2100
        self.fRS2 = self.f1020 - self.f0120
        self.fRS1 = 2 * self.f1020 - self.f1011
        self.fRS1 = 2 * self.f0120 - self.f0111
        

    def chiterms(self, ListOfBPMS=None):
        '''
         Add chi terms to the twiss table
        '''
        if ListOfBPMS is None:
            ListOfBPMS = []
        
        factMADtoSix = 0.0005
        self.chi3000 = []
        self.chi4000 = []
        self.chi2000 = []
        NAME = getattr(self, "NAME")
        S = getattr(self, "S")
        MUX = getattr(self, "MUX")
        K1L = getattr(self, "K1L")
        K2L = getattr(self, "K2L")
        K3L = getattr(self, "K3L")
        BETX = getattr(self, "BETX")
        
        if len(ListOfBPMS) == 0:
            print "Assuming that BPM elements are named as BP and H"
            for el in NAME:
                if "BP" in el and "H" in el:
                    ListOfBPMS.append(el)

        print "Found ", len(ListOfBPMS), "BPMs for chiterms computation"
        if len(ListOfBPMS) < 3:
            print "Error, not enough H BPMs in ListOfBPMs"
            sys.exit(1)

        self.chi = []
        self.chiBPMs = []
        self.chiS = []
        for i in range(len(ListOfBPMS) - 2):
            name = ListOfBPMS[i]
            name1 = ListOfBPMS[i + 1]
            name2 = ListOfBPMS[i + 2]
            self.chiBPMs.append([name, name1, name2])
            indx = self.indx[name]
            indx1 = self.indx[name1]
            indx2 = self.indx[name2]
            bphmii = MUX[indx]
            bphmii1 = MUX[indx1]
            bphmii2 = MUX[indx2]
            bphs = S[indx]
            bphs1 = S[indx1]
            bphs2 = S[indx2]
            self.chiS.append([bphs, bphs1, bphs2])
            d1 = (bphmii1 - bphmii) * 2 * PI - PI / 2
            d2 = (bphmii2 - bphmii1) * 2 * PI - PI / 2
            f1 = numpy.sqrt(1 + (numpy.sin(d1) / numpy.cos(d1)) ** 2)
            f2 = numpy.sqrt(1 + (numpy.sin(d2) / numpy.cos(d2)) ** 2)
            quadr = 0
            quadi = 0
            sexr = 0
            sexi = 0
            octr = 0
            octi = 0
            for j in range(len(NAME)):
                k1l = K1L[j]
                k2l = K2L[j]
                k3l = K3L[j]
                bx = BETX[j]
                m = MUX[j]
                if S[j] > bphs and S[j] < bphs1 and k2l ** 2 > 0:
                    quadr += numpy.cos(-1 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI) * k1l * bx ** 1 * f1
                    quadi += numpy.sin(-1 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI) * k1l * bx ** 1 * f1

                    sexr += numpy.cos(-2 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI) * k2l * bx ** 1.5 * f1
                    sexi += numpy.sin(-2 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI) * k2l * bx ** 1.5 * f1

                    octr += numpy.cos(-3 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI) * k3l * bx ** 2 * f1
                    octi += numpy.sin(-3 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI) * k3l * bx ** 2 * f1

                if S[j] > bphs1 and S[j] < bphs2 and k2l ** 2 > 0:
                    quadr += numpy.cos(-1 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI - d1 - d2) * k1l * bx ** 1 * f2
                    quadi += numpy.sin(-1 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI - d1 - d2) * k1l * bx ** 1 * f2

                    sexr += numpy.cos(-2 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI - d1 - d2) * k2l * bx ** 1.5 * f2
                    sexi += numpy.sin(-2 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI - d1 - d2) * k2l * bx ** 1.5 * f2

                    octr += numpy.cos(-3 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI - d1 - d2) * k3l * bx ** 2 * f2
                    octi += numpy.sin(-3 * (m - bphmii) * 2 * PI) * numpy.sin((m - bphmii) * 2 * PI - d1 - d2) * k3l * bx ** 2 * f2

                if S[j] > bphs2:
                    break
            self.chi.append(complex(sexr, sexi) / 4 * factMADtoSix)
            self.chi4000.append(complex(octr, octi) / 4 * factMADtoSix)
            self.chi2000.append(complex(quadr, quadi) / 4 * factMADtoSix)

    def Cmatrix(self):
        '''
         Calculate the C matrix
        '''
        self.C = []
        self.gamma = []
        self.f1001 = []
        self.f1010 = []
        S = getattr(self, "S")
        R11 = getattr(self, "R11")
        R12 = getattr(self, "R12")
        R21 = getattr(self, "R21")
        R22 = getattr(self, "R22")
        BETX = getattr(self, "BETX")
        BETY = getattr(self, "BETY")
        ALFX = getattr(self, "ALFX")
        ALFY = getattr(self, "ALFY")

        J = numpy.reshape(numpy.array([0, 1, -1, 0]), (2, 2))
        for j in range(0, len(S)):
            R = numpy.array([[R11[j], R12[j]], [R21[j], R22[j]]])
            
            C = matrixmultiply(-J, matrixmultiply(numpy.transpose(R), J))
            C = (1 / numpy.sqrt(1 + determinant(R))) * C

            g11 = 1 / numpy.sqrt(BETX[j])
            g12 = 0
            g21 = ALFX[j] / numpy.sqrt(BETX[j])
            g22 = numpy.sqrt(BETX[j])
            Ga = numpy.reshape(numpy.array([g11, g12, g21, g22]), (2, 2))

            g11 = 1 / numpy.sqrt(BETY[j])
            g12 = 0
            g21 = ALFY[j] / numpy.sqrt(BETY[j])
            g22 = numpy.sqrt(BETY[j])
            Gb = numpy.reshape(numpy.array([g11, g12, g21, g22]), (2, 2))
            C = matrixmultiply(Ga, matrixmultiply(C, inverse(Gb)))
            gamma = 1 - determinant(C)
            self.gamma.append(gamma)
            C = numpy.ravel(C)
            self.C.append(C)
            self.f1001.append(((C[0] + C[3]) * 1j + (C[1] - C[2])) / 4 / gamma)
            self.f1010.append(((C[0] - C[3]) * 1j + (-C[1] - C[2])) / 4 / gamma)

        self.F1001R = numpy.array(self.f1001).real
        self.F1001I = numpy.array(self.f1001).imag
        self.F1010R = numpy.array(self.f1010).real
        self.F1010I = numpy.array(self.f1010).imag
        self.F1001W = numpy.sqrt(self.F1001R ** 2 + self.F1001I ** 2)
        self.F1010W = numpy.sqrt(self.F1010R ** 2 + self.F1010I ** 2)

    def beatMatrix(self):
        '''
         Add RM to the twiss table
        '''
        self.RM = []
        S = getattr(self, "S")
        MUX = getattr(self, "MUX")
        MUY = getattr(self, "MUY")
        Q1 = getattr(self, "Q1")
        Q2 = getattr(self, "Q2")
        BETX = getattr(self, "BETX")
        BETY = getattr(self, "BETY")
        
        for j in range(0, len(S)):
            self.RM.append(-BETX * numpy.cos(2 * PI * (Q1 - 2 * abs(MUX[j] - MUX))) / numpy.sin(2 * PI * Q1))
        for j in range(0, len(S)):
            self.RM.append(-BETY * numpy.cos(2 * PI * (Q2 - 2 * abs(MUY[j] - MUY))) / numpy.sin(2 * PI * Q2))
        self.RM = numpy. array(self.RM)

    def abh(self, bet1, alf1, KL, K):
        """ bet1 and alf1 at the end of the element """
        gamma1 = (1. + alf1 ** 2) / bet1
        KL2 = 2.*KL
        sinhc = numpy.sinh(KL2) / KL2
        res = 0.5 * bet1 * (1. + sinhc) + alf1 * numpy.sinh(KL) ** 2. / KL / K + (sinhc - 1.) / (2.*K ** 2.) * gamma1
        return res

    def ab(self, bet1, alf1, KL, K):
        """ bet1 and alf1 at the end of the element """
        gamma1 = (1. + alf1 ** 2) / bet1
        KL2 = 2.*KL
        sinc = numpy.sin(KL2) / KL2
        res = 0.5 * bet1 * (1. + sinc) + alf1 * numpy.sin(KL) ** 2. / KL / K + (1. - sinc) / (2.*K ** 2.) * gamma1
        return res

    def AveBetas(self):
        totx = 0
        toty = 0
        totl = 0
        S = getattr(self, "S")
        L = getattr(self, "L")
        K1L = getattr(self, "K1L")
        BETX = getattr(self, "BETX")
        BETY = getattr(self, "BETY")
        ALFX = getattr(self, "ALFX")
        ALFY = getattr(self, "ALFY")
        NAME = getattr(self, "NAME")
        
        for i in range(len(S)):
            if L[i] > 0:
                k = numpy.sqrt(abs(K1L[i]) / L[i])
                kL = k * L[i]
                if K1L[i] == 0:
                    bxs = BETX[i] * L[i] + ALFX[i] * L[i] ** 2 + (1 + ALFX[i] ** 2) / BETX[i] * L[i] ** 3 / 3.
                    bys = BETY[i] * L[i] + ALFY[i] * L[i] ** 2 + (1 + ALFY[i] ** 2) / BETY[i] * L[i] ** 3 / 3.
                    totx = totx + bxs
                    toty = toty + bys
                    totl = totl + L[i]
                    print NAME[i], S[i], L[i], bxs, bys
                if K1L[i] > 0.0:
                    bxs = self.ab(BETX[i], ALFX[i], kL, k) * L[i]
                    totx = totx + bxs
                    bys = self.abh(BETY[i], ALFY[i], kL, k) * L[i]
                    toty = toty + bys
                    totl = totl + L[i]
                    abx=self.ab(BETX[i], ALFX[i], kL, k)
                    aby= self.abh(BETY[i], ALFY[i], kL, k)
                    print NAME[i], S[i], L[i], abx, aby, abx*K1L[i],aby*K1L[i] 
                if K1L[i] < 0:
                    bxs = self.abh(BETX[i], ALFX[i], kL, k) * L[i]
                    totx = totx + bxs
                    bys = self.ab(BETY[i], ALFY[i], kL, k) * L[i]
                    toty = toty + bys
                    totl = totl + L[i]
                    abx=self.abh(BETX[i], ALFX[i], kL, k)
                    aby=self.ab(BETY[i], ALFY[i], kL, k)
                    print NAME[i], S[i], L[i], abx, aby, abx*K1L[i],aby*K1L[i] 
            else:
                print NAME[i], S[i], L[i], BETX[i], BETY[i]
            print "TOTAL", S[i], totl, totx, toty

    def I5(self):
        H = 0
        NAME = getattr(self, "NAME")
        DX = getattr(self, "DX")
        DPX = getattr(self, "DPX")
        BETX = getattr(self, "BETX")
        ALFX = getattr(self, "ALFX")
        ANGLE = getattr(self, "ANGLE")
        L = getattr(self, "L")
        
        for i in range(len(NAME)):
            H = H + (DX[i] ** 2 + (DPX[i] * BETX[i] + DX[i] * ALFX[i]) ** 2) / BETX[i] * (abs(ANGLE[i])) ** 3 / L[i] ** 2

        return H
