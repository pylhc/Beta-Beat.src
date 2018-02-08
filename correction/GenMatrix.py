
import datetime
import time
import os

import numpy as np
from numpy.linalg import pinv as  generalized_inverse
from numpy import dot as matrixmultiply
import Python_Classes4MAD.metaclass as metaclass

def MakeDispList(x, m, modelcut=0.1, errorcut=0.027):
    intersected_names = []
    num_of_removed_bpms = 0
    if x == []:
        return []
    keys = x.__dict__.keys()
    ndmdl = "NDXMDL"
    STD = x.STDNDX
    NDX = x.NDX
    print "Number of x BPMs", len(x.NAME)
    for i in range(len(x.NAME)):
        ndm = x.__dict__[ndmdl][i]
        if (STD[i] < errorcut and abs(NDX[i] - ndm) < modelcut):
            try:
                m.indx[x.NAME[i].upper()]
            except KeyError:
                print "Not in Response:", x.NAME[i].upper()
                num_of_removed_bpms += 1
            else:
                intersected_names.append(x.NAME[i])
        else:
            num_of_removed_bpms += 1
    if num_of_removed_bpms > 0:
        print "Warning: ", num_of_removed_bpms, "BPMs removed from data for not being in the model or having too large error"
    return intersected_names

def MakeBetaList(x, m, modelcut=0.1, errorcut=0.027):      #Beta-beating and its error RELATIVE as shown in GUI
    intersected_names = []
    num_of_removed_bpms = 0
    keys = x.__dict__.keys()
    if "BETY" in keys:
        bmdl = "BETYMDL"
        BET = x.BETY
        if hasattr(x, 'ERRBETY'):
            ERR = x.ERRBETY
            if hasattr(x, "STDBETY"):
                STD = x.STDBETY
            else:
                STD = np.zeros(len(BET))
        else:
            ERR = np.zeros(len(BET))
            STD = x.BETYSTD
    else:
        bmdl = "BETXMDL"
        BET = x.BETX
        if hasattr(x, 'ERRBETX'):
            ERR = x.ERRBETX
            if hasattr(x, "STDBETX"):
                STD = x.STDBETX
            else:
                STD = np.zeros(len(BET))
        else:
            STD = x.BETXSTD
            ERR = np.zeros(len(BET))

    print "Number of x BPMs", len(x.NAME)
    for i in range(len(x.NAME)):
        bm = x.__dict__[bmdl][i]
        relerr = np.sqrt(STD[i]**2 + ERR[i]**2) / bm
        if (relerr < errorcut and (abs(BET[i] - bm)/bm) < modelcut):
            try:
                m.indx[x.NAME[i].upper()]
            except KeyError:
                print "Not in Response:", x.NAME[i].upper()
                num_of_removed_bpms += 1
            else:
                intersected_names.append(x.NAME[i])
        else:
            num_of_removed_bpms += 1
    if num_of_removed_bpms > 0:
        print "Warning: ", num_of_removed_bpms, "BPMs removed from data for not being in the model or having too large error"
    return intersected_names


def MakePairs(x, m, modelcut=0.1, errorcut=0.027):
    t = []
    num_removed_bpm_pairs = 0
    keys = x.__dict__.keys()
    if "PHYMDL" in keys:
        phmdl = "PHYMDL"
        STD = x.STDPHY
        PHASE = x.PHASEY

    else:
        phmdl = "PHXMDL"
        STD = x.STDPHX
        PHASE = x.PHASEX

    for i in range(len(x.NAME)-1):    # The -1 comes from the fact that the last BPM is again the first BPM, model not ready for this
        phm = x.__dict__[phmdl][i]

        if (STD[i] < errorcut and abs(PHASE[i]-phm) < modelcut):
            try:
                m.indx[x.NAME[i].upper()]
                m.indx[x.NAME2[i].upper()]
            except KeyError:
                print "Not in Response:", x.NAME[i].upper(), x.NAME2[i].upper()
                num_removed_bpm_pairs += 1
            else:
                t.append(x.NAME[i]+' '+x.NAME2[i])
        else:
            num_removed_bpm_pairs += 1
    if num_removed_bpm_pairs > 0:
        print "Warning: ", num_removed_bpm_pairs, "BPM pairs removed from data for not beeing in the model or having too large error deviations: ", phmdl, modelcut, "STDPH", errorcut, "LEN", len(t)
    return t





def writeparams(deltafamilie, variables, accel_path, app=0, path="./"):
    if (app == 0):
        mode = 'w'
    if (app == 1):
        mode = 'a'
    a = datetime.datetime.fromtimestamp(time.time())
    if accel_path[-5:] == "LHCB1":
        path_B1tfs = os.path.join(accel_path, "varsKLvsK_B1.tfs")
        twiss_kl = metaclass.twiss(path_B1tfs)
    if accel_path[-5:] == "LHCB2":
        path_B2tfs = os.path.join(accel_path, "varsKLvsK_B2.tfs")
        twiss_kl = metaclass.twiss(path_B2tfs)
    g = open (path+'/changeparameters', mode)
    f = open (path+'/changeparameters.tfs', mode)
    print >> f, "@", "APP", "%le", app
    print >> f, "@", "PATH", "%s", path
    print >> f, "@", "DATE", "%s", a.ctime()
    print >> f, "*", "NAME", "DELTA"
    print >> f, "$", "%s", "%le"
    for i, var in enumerate(variables):
        if var[0]=="l":
            g.write(var[1:]+' = '+ var[1:] +' + ( '+str((deltafamilie[i]/twiss_kl.LENGTH[twiss_kl.indx[var]]))+' );\n')
            f.write(var[1:]+'   '+str((deltafamilie[i]/twiss_kl.LENGTH[twiss_kl.indx[var]]))+'\n')
        else:    
            g.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
            f.write(var+'   '+str(deltafamilie[i])+'\n')
    g.close()
    f.close()

def betabeat(a, b):
    rmsx = np.sqrt(  sum( ((a.BETX-b.BETX)/b.BETX)**2 )  /  len(b.BETX)  )
    rmsy = np.sqrt(  sum( ((a.BETY-b.BETY)/b.BETY)**2 )  /  len(b.BETY)  )
    peakx = np.max( abs((a.BETX-b.BETX)/b.BETX) )
    peaky = np.max( abs((a.BETY-b.BETY)/b.BETY) )
    rmsdx = np.sqrt(sum(((a.DX/np.sqrt(a.BETX)-b.DX/np.sqrt(b.BETX)))**2)/len(b.DX))
    peakdx = np.max(abs((a.DX/np.sqrt(a.BETX)-b.DX/np.sqrt(b.BETX))))
    dphixa = []
    dphiya = []
    dphixb = []
    dphiyb = []
    for i in range(1, len(a.MUX)):
        dphixa.append(a.MUX[i]-a.MUX[i-1])
        dphiya.append(a.MUY[i]-a.MUY[i-1])
        dphixb.append(b.MUX[i]-b.MUX[i-1])
        dphiyb.append(b.MUY[i]-b.MUY[i-1])
    dphixa = np.array(dphixa)
    dphiya = np.array(dphiya)
    dphixb = np.array(dphixb)
    dphiyb = np.array(dphiyb)
    rmsphix = np.sqrt(np.sum((dphixa-dphixb)**2)/len(dphixa))
    rmsphiy = np.sqrt(np.sum((dphiya-dphiyb)**2)/len(dphiya))
    peakphix = np.max(abs(dphixa-dphixb))
    peakphiy = np.max(abs(dphiya-dphiyb))

    return np.array([rmsx, rmsy, peakx, peaky, rmsphix, rmsphiy,
                  peakphix, peakphiy, rmsdx, peakdx])


def betabeatEXP(x, y, b):
    rmsx = 0
    rmsy = 0
    peakx = 0
    peaky = 0
    rmsdx = 0
    peakdx = 0
    phasexlist = MakePairs(x, b)
    phaseylist = MakePairs(y, b)

    dphixa = []
    dphiya = []
    dphixb = []
    dphiyb = []
    for ph in phasexlist:
        [bpm1, bpm2] = ph.upper().split()
        dphixa.append(x.PHASE[x.indx[bpm1]])
        dphixb.append(b.MUX[b.indx[bpm2]]-b.MUX[b.indx[bpm1]])

    for ph in phaseylist:
        [bpm1, bpm2] = ph.upper().split()
        dphiya.append(y.PHASE[y.indx[bpm1]])
        dphiyb.append(b.MUY[b.indx[bpm2]]-b.MUY[b.indx[bpm1]])
    dphixa = np.array(dphixa)
    dphiya = np.array(dphiya)
    dphixb = np.array(dphixb)
    dphiyb = np.array(dphiyb)
    if len(dphixa)>0:
        rmsphix = np.sqrt( sum((dphixa-dphixb)**2) / len(dphixa) )
        peakphix = np.max( abs(dphixa-dphixb) )
    else:
        rmsphix = 0
        peakphix = 0
        print "Warning: No horizontal BPMs matching when computing beta-beating"
    if len(dphiya) > 0:
        rmsphiy = np.sqrt(np.sum((dphiya-dphiyb)**2)/len(dphiya))
        peakphiy = np.max(abs(dphiya-dphiyb))
    else:
        rmsphiy = 0
        peakphiy = 0
        print "Warning: No vertical BPMs matching when computing beta-beating"

    return np.array([rmsx, rmsy, peakx, peaky, rmsphix, rmsphiy,
                  peakphix, peakphiy, rmsdx, peakdx])


def correctbeat(a, beat_input, cut=0.01, app=0, path="./"):
    R = np.transpose(beat_input.sensitivity_matrix)
    vector = beat_input.computevector(a)
    wg = beat_input.wg
    weisvec = np.array(np.concatenate([np.sqrt(wg[0])*np.ones(len(beat_input.phasexlist)),
							            np.sqrt(wg[1])*np.ones(len(beat_input.phaseylist)),
							            np.sqrt(wg[2])*np.ones(len(beat_input.betaxlist)),
							            np.sqrt(wg[3])*np.ones(len(beat_input.betaylist)),
							            np.sqrt(wg[4])*np.ones(len(beat_input.displist)),
							            np.sqrt(wg[5])*np.ones(2)]))
    Rnew = np.transpose(np.transpose(R)*weisvec)
    delta = -matrixmultiply(  generalized_inverse(Rnew, cut), (vector-beat_input.zerovector)/beat_input.normvector  )
    writeparams(delta, beat_input.varslist, app,  path=path)


def correctbeatEXP(x, y, dx, beat_input, cut=0.01, app=0, path="./", xbet=[], ybet=[]):
    R =   np.transpose(beat_input.sensitivity_matrix)
    vector = beat_input.computevectorEXP(x, y, dx, xbet, ybet)
    errwg = beat_input.errwg
    if errwg == 1:
        weisvec = beat_input.computeweightvectorEXP(x, y, dx, xbet, ybet)
    else:
        wg = beat_input.wg
        weisvec = np.array(np.concatenate([wg[0]*np.ones(len(beat_input.phasexlist)),
							            wg[1]*np.ones(len(beat_input.phaseylist)),
							            wg[2]*np.ones(len(beat_input.betaxlist)),
							            wg[3]*np.ones(len(beat_input.betaylist)),
							            wg[4]*np.ones(len(beat_input.displist)),
							            wg[5]*np.ones(2)]))

    Rnew = np.transpose(np.transpose(R)*weisvec)
    delta = -matrixmultiply(generalized_inverse(Rnew, cut), (vector-beat_input.zerovector)*weisvec/beat_input.normvector)
    
    writeparams(delta, beat_input.varslist, beat_input.accel_path, app, path)
    return [delta, beat_input.varslist]


class beat_input:
    def __init__(self, varslist, accel_path, phasexlist=[], phaseylist=[], betaxlist=[], betaylist=[], displist=[], wg=[1,1,1,1,1,10], errwg = 1):
        self.varslist = varslist
        self.phasexlist = phasexlist
        self.phaseylist = phaseylist
        self.betaxlist = betaxlist
        self.betaylist = betaylist
        print "Length at entry level of beta ", len(self.betaxlist), len(self.phasexlist)
        self.displist = displist
        self.wg = wg
        self.errwg = errwg
        self.accel_path = accel_path
        self.zerovector = []
        self.normvector = []
        self.sensitivity_matrix = []
    
    def computeweightvectorEXP(self, x, y, dx, xbet=None, ybet=None):
        phix = []
        phiy = []
        betx = []
        bety = []
        disp = []
        tune = []
        wg = self.wg
        print "Error weighed vector used"
        for ph in self.phasexlist:
            bpm1 = ph.split()[0]
            phix.append(wg[0] / x.STDPHX[x.indx[bpm1]])

        for ph in self.phaseylist:
            bpm1 = ph.split()[0]
            phiy.append(wg[1] / y.STDPHY[y.indx[bpm1]])
        for b in self.displist:
            disp.append(wg[4] / dx.STDNDX[dx.indx[b]])
        if xbet is not None:
            for b in self.betaxlist:
                try:
                    betx.append(wg[2] * xbet.BETXMDL[xbet.indx[b]] / np.sqrt(xbet.STDBETX[xbet.indx[b]]**2 + xbet.ERRBETX[xbet.indx[b]]**2))
                except AttributeError:
                    betx.append(wg[2] * xbet.BETXMDL[xbet.indx[b]] / xbet.ERRBETX[xbet.indx[b]])
        if ybet is not None:
            for b in self.betaylist:
                try:
                    bety.append(wg[3] * ybet.BETYMDL[ybet.indx[b]] / np.sqrt(ybet.STDBETY[ybet.indx[b]]**2 + ybet.ERRBETY[ybet.indx[b]]**2))
                except AttributeError:
                    bety.append(wg[3] * ybet.BETYMDL[ybet.indx[b]] / ybet.ERRBETY[ybet.indx[b]])   
        tune.append(wg[5]/0.001)  #Hardcoded precision of tunes, do we have the proper error somewhere?
        tune.append(wg[5]/0.001)

        return np.array(np.concatenate([phix, phiy, betx, bety, disp, tune]))


    def computevectorEXP(self, x, y, dx, xbet=None, ybet=None):
        phix = []
        phiy = []
        betx = []
        bety = []
        disp = []
        print "Exp tunes in beta_input class: ", x.Q1, y.Q2
        for ph in self.phasexlist:
            bpm1 = ph.split()[0]
            phix.append(x.PHASEX[x.indx[bpm1]])

        for ph in self.phaseylist:
            bpm1 = ph.split()[0]
            phiy.append(y.PHASEY[y.indx[bpm1]])
        for b in self.displist:
            disp.append(dx.NDX[dx.indx[b]])
        if xbet is not None:
            for b in self.betaxlist:
                betx.append(xbet.BETX[xbet.indx[b]])
        if ybet is not None:
            for b in self.betaylist:
                bety.append(ybet.BETY[ybet.indx[b]])

        return np.array(np.concatenate([phix, phiy, betx, bety, disp, [x.Q1], [y.Q2]]))


    def computevector(self,x):
        phix = []
        phiy = []
        betx = []
        bety = []
        disp = []
        for ph in self.phasexlist:
            [bpm1, bpm2] = ph.upper().split()
            phix.append(x.MUX[x.indx[bpm2]]-x.MUX[x.indx[bpm1]])
        for ph in self.phaseylist:
            [bpm1, bpm2] = ph.upper().split()
            phiy.append(x.MUY[x.indx[bpm2]]-x.MUY[x.indx[bpm1]])
        for b in self.betaxlist:
            betx.append(x.BETX[x.indx[b]])
        for b in self.betaylist:
            bety.append(x.BETY[x.indx[b]])
        for b in self.displist:
            disp.append(x.DX[x.indx[b]]/np.sqrt(x.BETX[x.indx[b]]))
        return np.array(np.concatenate([phix, phiy, betx, bety, disp, [x.Q1], [x.Q2]]))


    def computeSensitivityMatrix(self, x):
        self.zerovector = self.computevector(x['0'])
        if "incr_dict" in x:
            incr_dict = x["incr_dict"]
        else:
            incr = x['incr'][0]

        nph = len(self.phasexlist) + len(self.phaseylist)
        nbet = len(self.betaxlist) + len(self.betaylist)
        print "Length of beta vector ", nph
        ndisp = len(self.displist)
        self.normvector = np.array(np.concatenate([np.ones(len(self.phasexlist)), np.ones(len(self.phaseylist)), self.zerovector[nph:nph+nbet], np.ones(ndisp+2)]))*1.0  #To take dbeta/beta later
        for var in self.varslist:
            if "incr_dict" in x:
                incr = incr_dict[var]
            vector = self.computevector(x[var])
            self.sensitivity_matrix.append( (vector-self.zerovector)/self.normvector/incr )
        self.sensitivity_matrix = np.array(self.sensitivity_matrix)
        return self.sensitivity_matrix


