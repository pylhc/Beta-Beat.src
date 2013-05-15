from math import *
from optparse import OptionParser
from string import split
import sys, copy


parser = OptionParser()
#parser.add_option("-f", "--file",
#    help="input file",
#    metavar="FILE", default="all.csv", dest="file")

parser.add_option("-Q", "--Tune",
    help="Tune in usual units, default=0.31",
    metavar="TUNE", default="0.31", dest="Q")

parser.add_option("-L", "--Lstar",
    help="L* in meters, default=22.965000",
    metavar="LSTAR", default="22.965000", dest="Lstar")

parser.add_option("-K", "--K1",
    help="K of MQXA1,absolute value in m^-2, default=0.008730196766",
    metavar="K", default="0.008730196766", dest="K")

parser.add_option("-D", "--DeltaK",
    help="DeltaK used during K-mod m^-2",
    metavar="DK", default="", dest="DK")

parser.add_option("-l", "--quadlength",
    help="Length of MQXA1,default=6.37" ,
    metavar="LENGTH", default="6.37", dest="l")

parser.add_option("-d", "--tuneshifts",
    help="tuneshifts from left and right quads separated by commas",
    metavar="DQ", default="0,0", dest="dqs")

parser.add_option("-e", "--errorintunes",
    help="Error bar of tuneshifts from left and right quads separated by commas",
    metavar="EDQ", default="0,0", dest="edqs")

parser.add_option("-t", "--label",
        help="Measurement lable, B1VIP1, B2HIP5, etc",
        metavar="LABEL", default="B", dest="label")


parser.add_option("-b", "--beta",
        help="Design beta* in order to compute ratio",
        metavar="BETA", default="0.6", dest="beta")



def abh(b,L,w,KL,K):
    bet1=b+(L-w)**2/b
    alf1=(L-w)/b
    KL2=2*KL
    sinhc=sinh(KL2)/KL2
    res=0.5*bet1*(1+sinhc)+alf1*sinh(KL)**2/KL/K+(sinhc-1)/(2.*b*K**2)
    return res


def ab(b,L,w,KL,K):
    bet1=b+(L-w)**2/b
    alf1=(L-w)/b
    KL2=2*KL
    sinc=sin(KL2)/KL2
    res=0.5*bet1*(1+sinc)+alf1*sin(KL)**2/KL/K+(1-sinc)/(2.*b*K**2)
    return res



class Simplex:
    def __init__(self, testfunc, guess, increments, kR = -1, kE = 2, kC = 0.5):
        """Initializes the simplex.
        INPUTS
        ------
        testfunc      the function to minimize
        guess[]       an list containing initial guesses
        increments[]  an list containing increments, perturbation size
        kR            reflection constant  (alpha =-1.0)
        kE            expansion constant   (gamma = 2.0)
        kC            contraction constant (beta  = 0.5)
        """
        self.testfunc = testfunc
        self.guess = guess
        self.increments = increments
        self.kR = kR
        self.kE = kE
        self.kC = kC
        self.numvars = len(self.guess)
        self.simplex = []

        self.lowest = -1
        self.highest = -1
        self.secondhighest = -1

        self.errors = []
        self.currenterror = 0

        # Initialize vertices
        # MV: the first vertex is just the initial guess
        #     the other N vertices are the initial guess plus the individual increments
        #     the last two vertices will store the centroid and the reflected point
        #     the compute errors at the ... vertices
        
        for vertex in range(0, self.numvars + 3):
            self.simplex.append(copy.copy(self.guess))

        for vertex in range(0, self.numvars + 1):
            for x in range(0, self.numvars):
                if x == (vertex - 1):
                    self.simplex[vertex][x] = self.guess[x] + self.increments[x]
            self.errors.append(0)
        self.calculate_errors_at_vertices()

    def minimize(self, epsilon = 0.0001, maxiters = 400, monitor = 0):
        """Walks to the simplex down to a local minima.
        INPUTS
        ------
        epsilon       convergence requirement
        maxiters      maximum number of iterations
        monitor       if non-zero, progress info is output to stdout  

        OUTPUTS
        -------
        an array containing the final values
        lowest value of the error function
        number of iterations taken to get here
        """
        
        iter = 0
        
        for iter in range(0, maxiters):
            # Identify highest and lowest vertices
            
            self.highest = 0
            self.lowest = 0
            for vertex in range(1, self.numvars + 1):
                if self.errors[vertex] > self.errors[self.highest]:
                    self.highest = vertex
                if self.errors[vertex] < self.errors[self.lowest]:
                    self.lowest = vertex

            # Identify second-highest vertex

            self.secondhighest = self.lowest
            for vertex in range(0, self.numvars + 1):
                if vertex == self.highest:
                    continue
                elif vertex == self.secondhighest:
                    continue
                elif self.errors[vertex] > self.errors[self.secondhighest]:
                    self.secondhighest = vertex

            # Test for convergence:
            #   compute the average merit figure (ugly)

            S = 0.0
            for vertex in range(0, self.numvars + 1):
                S = S + self.errors[vertex]
            F2 = S / (self.numvars + 1)

            #   compute the std deviation of the merit figures (ugly)

            S1 = 0.0
            for vertex in range(0, self.numvars + 1):
                S1 = S1 + (self.errors[vertex] - F2)**2
            T = sqrt(S1 / self.numvars)
            
            # Optionally, print progress information

            if monitor:
                print '\r' + 72 * ' ',
                print '\rIteration = %d   Best = %f   Worst = %f' % \
                      (iter,self.errors[self.lowest],self.errors[self.highest]),
                sys.stdout.flush()
                
            if T <= epsilon:
                # We converged!  Break out of loop!
                
                break;
            else:
                # Didn't converge.  Keep crunching.
                
                # Calculate centroid of simplex, excluding highest vertex
                # store centroid in element N+1

                # loop over coordinates
                for x in range(0, self.numvars):
                    S = 0.0
                    for vertex in range(0, self.numvars + 1):
                        if vertex == self.highest:
                            continue
                        S = S + self.simplex[vertex][x]
                    self.simplex[self.numvars + 1][x] = S / self.numvars

                # reflect the simplex across the centroid
                # store reflected point in elem. N + 2 (and self.guess)
                
                self.reflect_simplex()
                self.currenterror = self.testfunc(self.guess)

                if self.currenterror < self.errors[self.highest]:
                    self.accept_reflected_point()

                if self.currenterror <= self.errors[self.lowest]:
                    self.expand_simplex()
                    self.currenterror = self.testfunc(self.guess)

                    # at this point we can assume that the highest
                    # value has already been replaced once
                    if self.currenterror < self.errors[self.highest]:
                        self.accept_expanded_point()
                elif self.currenterror >= self.errors[self.secondhighest]:
                    # worse than the second-highest, so look for
                    # intermediate lower point

                    self.contract_simplex()
                    self.currenterror = self.testfunc(self.guess)

                    if self.currenterror < self.errors[self.highest]:
                        self.accept_contracted_point()
                    else:
                        self.multiple_contract_simplex()
                        
        # Either converged or reached the maximum number of iterations.
        # Return the lowest vertex and the currenterror.

        for x in range(0, self.numvars):
            self.guess[x] = self.simplex[self.lowest][x]
        self.currenterror = self.errors[self.lowest]
        return self.guess, self.currenterror, iter

    # same as expand, but with alpha < 1; kC = 0.5 fine with NR

    def contract_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kC * self.simplex[self.highest][x] + (1 - self.kC) * self.simplex[self.numvars + 1][x]
        return

    # expand: if P is vertex and Q is centroid, alpha-expansion is Q + alpha*(P-Q),
    #         or (1 - alpha)*Q + alpha*P; default alpha is 2.0; agrees with NR
    def expand_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kE * self.guess[x]                 + (1 - self.kE) * self.simplex[self.numvars + 1][x]
        return

    # reflect: if P is vertex and Q is centroid, reflection is Q + (Q-P) = 2Q - P,
    #          which is achieved for kR = -1 (default value); agrees with NR
    def reflect_simplex(self):
        # loop over variables
        for x in range(0, self.numvars):
            self.guess[x] = self.kR * self.simplex[self.highest][x] + (1 - self.kR) * self.simplex[self.numvars + 1][x]
            # store reflected point in elem. N + 2
            self.simplex[self.numvars + 2][x] = self.guess[x]
        return

    # multiple contraction: around the lowest point; agrees with NR

    def multiple_contract_simplex(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                continue
            for x in range(0, self.numvars):
                self.simplex[vertex][x] = 0.5 * (self.simplex[vertex][x] + self.simplex[self.lowest][x])
        self.calculate_errors_at_vertices()
        return

    def accept_contracted_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        return

    def accept_expanded_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        return

    def accept_reflected_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.simplex[self.numvars + 2][x]
        return

    def calculate_errors_at_vertices(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                # compute the error unless we're the lowest vertex
                continue
            for x in range(0, self.numvars):
                self.guess[x] = self.simplex[vertex][x]
            self.currenterror = self.testfunc(self.guess)
            self.errors[vertex] = self.currenterror
        return

global mindb, maxdb
def chi2(d):
    global L,KL,K,mindb,maxdb
    #print mindb, maxdb
    b=d[0]
    w=d[1]    
    c2= (ab(b,L,w,KL,K)-mindb)**2+(abh(b,L,-w,KL,K)-maxdb)**2
    
    
    #print c2
    return c2






(options, args) = parser.parse_args()


L=float(options.Lstar)
K=sqrt(float(options.K))
l=float(options.l)
KL=K*l
dqs=options.dqs.split(",")
dql=float(dqs[0])
dqr=float(dqs[1])
edqs=options.edqs.split(",")
edql=float(edqs[0])
edqr=float(edqs[1])
DK=float(options.DK)
Q=float(options.Q)
guess=[1000,-0.]
incr=[0.01,0.001]




def betasfromtunes(Tdql,Tdqr):
    global Q,dql,dqr,l,DK,maxdb,mindb
    dbl= 2*(1/tan(2*pi*Q)*(1-cos(2*pi*Tdql))+sin(2*pi*Tdql))/(l*DK)
    dbr= 2*(1/tan(2*pi*Q)*(1-cos(2*pi*Tdqr))+sin(2*pi*Tdqr))/(l*DK)
    mindb=min(abs(dbl),abs(dbr))
    maxdb=max(abs(dbl),abs(dbr))
    return mindb, maxdb


idealab=ab(float(options.beta),L,0.0,KL,K)
idealabh=abh(float(options.beta),L,0.0,KL,K)




resb=[]
resw=[]
resabbeat=[]
resabhbeat=[]

mindb, maxdb=betasfromtunes(dql,dqr)
s = Simplex(chi2, guess, incr)
values, err, iter = s.minimize(epsilon=0.00000001,maxiters=10000)

resb.append(values[0])
resw.append(values[1])
abbeat= (mindb-idealab)/idealab
abhbeat= (maxdb-idealabh)/idealabh
resabbeat.append(abbeat)
resabhbeat.append(abhbeat)


mindb, maxdb=betasfromtunes(dql+edql,dqr)
s = Simplex(chi2, guess, incr)
values, err, iter = s.minimize(epsilon=0.00000001,maxiters=10000)

resb.append(values[0])
resw.append(values[1])
abbeat= (mindb-idealab)/idealab
abhbeat= (maxdb-idealabh)/idealabh
resabbeat.append(abbeat)
resabhbeat.append(abhbeat)

mindb, maxdb=betasfromtunes(dql-edql,dqr)
s = Simplex(chi2, guess, incr)
values, err, iter = s.minimize(epsilon=0.00000001,maxiters=10000)

resb.append(values[0])
resw.append(values[1])
abbeat= (mindb-idealab)/idealab
abhbeat= (maxdb-idealabh)/idealabh
resabbeat.append(abbeat)
resabhbeat.append(abhbeat)

mindb, maxdb=betasfromtunes(dql,dqr+edqr)
s = Simplex(chi2, guess, incr)
values, err, iter = s.minimize(epsilon=0.00000001,maxiters=10000)

resb.append(values[0])
resw.append(values[1])
abbeat= (mindb-idealab)/idealab
abhbeat= (maxdb-idealabh)/idealabh
resabbeat.append(abbeat)
resabhbeat.append(abhbeat)

mindb, maxdb=betasfromtunes(dql,dqr-edqr)
s = Simplex(chi2, guess, incr)
values, err, iter = s.minimize(epsilon=0.00000001,maxiters=10000)
resb.append(values[0])
resw.append(values[1])
abbeat= (mindb-idealab)/idealab
abhbeat= (maxdb-idealabh)/idealabh
resabbeat.append(abbeat)
resabhbeat.append(abhbeat)




print

print 'betastar = ', resb
print 'ws = ', resw

IPb=[]
for i in range(len(resb)):
    IPb.append(resb[i]+resw[i]**2/resb[i])

print 'IP betas=', IPb
#ave=sum(resb)/len(resb)
#std=sqrt(sum(map(lambda x: x*x,resb))/len(resb)-ave*ave)
#resb=IPb


stdIPb=sqrt(max(abs(IPb[1]-IPb[0]), abs(IPb[2]-IPb[0]))**2+
         max(abs(IPb[3]-IPb[0]), abs(IPb[4]-IPb[0]))**2)

std=sqrt(max(abs(resb[1]-resb[0]), abs(resb[2]-resb[0]))**2+
         max(abs(resb[3]-resb[0]), abs(resb[4]-resb[0]))**2)
stdw=sqrt(max(abs(resw[1]-resw[0]), abs(resw[2]-resw[0]))**2+
             max(abs(resw[3]-resw[0]), abs(resw[4]-resw[0]))**2)


stdab=sqrt(max(abs(resabbeat[1]-resabbeat[0]), abs(resabbeat[2]-resabbeat[0]))**2+
             max(abs(resabbeat[3]-resabbeat[0]), abs(resabbeat[4]-resabbeat[0]))**2)

stdabh=sqrt(max(abs(resabhbeat[1]-resabhbeat[0]), abs(resabhbeat[2]-resabhbeat[0]))**2+
             max(abs(resabhbeat[3]-resabhbeat[0]), abs(resabhbeat[4]-resabhbeat[0]))**2)


print 'IP beta = ', resb[0], "+/-", std
#print 'iterations = ', iter

f=open("BetaStarResults.dat","a")
print >>f, options.label, IPb[0], stdIPb, resw[0], stdw, resb[0], std, resabbeat[0], stdab,resabhbeat[0], stdabh
f.close()





