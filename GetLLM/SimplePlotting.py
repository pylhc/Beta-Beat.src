from __future__ import unicode_literals
import __init__  # @UnusedImport
from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from  Python_Classes4MAD import metaclass

Butter1 = '#FCE94F'
Butter2 = '#EDD400'
Butter3 = '#C4A000'
Orange1 = '#FCAF3E'
Orange2 = '#F57900'
Orange3 = '#CE5C00'
Chocolate1 = '#E9B96E'
Chocolate2 = '#C17D11'
Chocolate3 = '#8F5902'
Chameleon1 = '#8AE234'
Chameleon2 = '#73D216'
Chameleon3 = '#4E9A06'
SkyBlue1 = '#729FCF'
SkyBlue2 = '#3465A4'
SkyBlue3 = '#204A87'
Plum1 = '#AD7FA8'
Plum2 = '#75507B'
Plum3 = '#5C3566'
ScarletRed1 = '#EF2929'
ScarletRed2 = '#CC0000'
ScarletRed3 = '#A40000'
Aluminium1 = '#EEEEEC'
Aluminium2 = '#D3D7CF'
Aluminium3 = '#BABDB6'
Aluminium4 = '#888A85'
Aluminium5 = '#555753'
Aluminium6 = '#2E3436'

#colors = ['red', 'blue', 'green', 'yellow']
colors = [SkyBlue1,ScarletRed1, Butter1 , Aluminium1]
markeredgecolors = [SkyBlue3,ScarletRed3, Butter3,  Aluminium3]


def setParams(labelfontsize=15, legendfontsize=15):
    rc('font', **{'family': 'sans-serif', 'serif': ['Computer Modern']})
    params = {'backend': 'pdf',
          'axes.labelsize': labelfontsize,
          'font.size': labelfontsize,
          'legend.fontsize': 18,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
#          'axes.title': 12,
          'text.usetex': False,
          'axes.unicode_minus': True,
          'xtick.major.pad': 15,
          'ytick.major.pad': 15,
          'xtick.minor.pad': 15,
          'ytick.minor.pad': 15,
    }
    plt.rcParams.update(params)


def splitFiles(files):
    return files.split(',')


def setLimits(accel, vmin, vmax, hmin, hmax, plot):
    plot.set_ylim(float(vmin), float(vmax))
    plot.set_xlim(float(hmin), float(hmax))


def getTwiss(path):
    tw = metaclass.twiss(path)
    return tw


#relative beating and relative errors
def getWRms(betabeat, betabeaterr):
    meanerr = np.mean(betabeaterr)
    wrms = np.sqrt(np.average(np.square(betabeat), weights=1 / np.square(betabeaterr)))
    print "Weighted rms: " + str(wrms) + "\tMean error: " + str(meanerr)
    return wrms, meanerr


def getweightrms(betabeat, betabeaterr):  # relative beating and relative errors
    return np.sqrt(np.average(np.square(betabeat), weights=1 / np.square(betabeaterr)))


def getbeta(path, plane, mainnode):
    if (mainnode == 'Beta'):
        if (not os.path.exists(os.path.join(path, "getbeta" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getbeta" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getbeta" + plane + "_free.out"))
        try:
            betabeaterr = np.sqrt(getattr(t, "STDBET" + plane.upper()) ** 2 + getattr(t, "ERRBET" + plane.upper()) ** 2) / getattr(t, "BET" + plane.upper() + "MDL")
        except AttributeError:
            betabeaterr = getattr(t, "ERRBET" + plane.upper()) / getattr(t, "BET" + plane.upper() + "MDL")
    elif(mainnode == 'Beta_Amp'):
        if (not os.path.exists(os.path.join(path, "getampbeta" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getampbeta" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getampbeta" + plane + "_free.out"))
        betabeaterr = getattr(t, "BET" + plane.upper() + "STD") / getattr(t, "BET" + plane.upper() + "MDL")
    s = t.S
    betabeat = (getattr(t, "BET" + plane.upper()) - getattr(t, "BET" + plane.upper() + "MDL")) / getattr(t, "BET" + plane.upper() + "MDL")
    print 'Beta-beating from phase plane ' + plane + ':'
    getWRms(betabeat, betabeaterr)
    return s, betabeat, betabeaterr


def getbetamdl(path, plane, mainnode):
    if (mainnode == 'Beta'):
        if (not os.path.exists(os.path.join(path, "getbeta" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getbeta" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getbeta" + plane + "_free.out"))
        try:
            err = np.sqrt(getattr(t, "STDBET" + plane.upper()) ** 2 + getattr(t, "ERRBET" + plane.upper()) ** 2) / getattr(t, "BET" + plane.upper() + "MDL")
        except AttributeError:
            #betabeaterr = errbet or errbet/betmdl?
            err = getattr(t, "ERRBET" + plane.upper())
    elif(mainnode == 'Beta_Amp'):
        if (not os.path.exists(os.path.join(path, "getampbeta" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getampbeta" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getampbeta" + plane + "_free.out"))
        err = getattr(t, "BET" + plane.upper() + "STD")
    s = t.S
    meas = getattr(t, "BET" + plane.upper())
    mdl = getattr(t, "BET" + plane.upper() + "MDL")
    return s, meas, mdl, err


def getcouple(path, plane, subnode):
    if (not os.path.exists(os.path.join(path, "getcouple_free.out"))):
        t = getTwiss(os.path.join(path, "getcouple.out"))
    else:
        t = getTwiss(os.path.join(path, "getcouple_free.out"))
    s = t.S
    if (subnode == "amp"):
        if(plane == "x"):
            f = t.F1001W
            fstd = t.FWSTD1
        if(plane == "y"):
            f = t.F1010W
            fstd = t.FWSTD2
    elif (subnode == "real"):
        if(plane == "x"):
            f = t.F1001R
            fstd = t.FWSTD1
        if(plane == "y"):
            f = t.F1010R
            fstd = t.FWSTD2
    elif (subnode == "imaginary"):
        if(plane == "x"):
            f = t.F1001I
            fstd = t.FWSTD1
        if(plane == "y"):
            f = t.F1010I
            fstd = t.FWSTD2
    return s, f, fstd


def getphasediff(path, plane, mainnode):
    if(mainnode == 'Phase_Total'):
        if (not os.path.exists(os.path.join(path, "getphasetot" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getphasetot" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getphasetot" + plane + "_free.out"))
    elif(mainnode == 'Phase'):
        if (not os.path.exists(os.path.join(path, "getphase" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getphase" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getphase" + plane + "_free.out"))
    s = t.S
    phdiff = getattr(t, "PHASE" + plane.upper()) - getattr(t, "PH" + plane.upper() + "MDL")
    pherr = getattr(t, "STDPH" + plane.upper())
    return s, phdiff, pherr


def getphasemdl(path, plane, mainnode):
    if(mainnode == 'Phase_Total'):
        if (not os.path.exists(os.path.join(path, "getphasetot" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getphasetot" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getphasetot" + plane + "_free.out"))
    elif(mainnode == 'Phase'):
        if (not os.path.exists(os.path.join(path, "getphase" + plane + "_free.out"))):
            t = getTwiss(os.path.join(path, "getphase" + plane + ".out"))
        else:
            t = getTwiss(os.path.join(path, "getphase" + plane + "_free.out"))
    s = t.S
    ph = getattr(t, "PHASE" + plane.upper())
    phstd = getattr(t, "STDPH" + plane.upper())
    phmdl = getattr(t, "PH" + plane.upper() + "MDL")
    return s, ph, phmdl, phstd


def getclosedorbit(path, plane):
    t = getTwiss(os.path.join(path, "getCO" + plane + ".out"))
    s = t.S
    if(plane == "x"):
        planedata = t.X
        std = t.STDX
    elif(plane == "y"):
        planedata = t.Y
        std = t.STDY
    else:
        print 'Unknown plane'
    return s, planedata, std


def getdisperssionmdl(path, plane):
    t = getTwiss(os.path.join(path, "getD" + plane + ".out"))
    s = t.S
    yvaluex1 = getattr(t, "D" + plane.upper())
    yvaluex2 = getattr(t, "D" + plane.upper() + "MDL")
    error = getattr(t, "STDD" + plane.upper())
    return s, yvaluex1, yvaluex2, error


def getdisperssiondiff(path, plane):
    t = getTwiss(os.path.join(path, "getD" + plane + ".out"))
    s = t.S
    diff = getattr(t, "D" + plane.upper()) - getattr(t, "D" + plane.upper() + "MDL")
    error = getattr(t, "STDD" + plane.upper())
    return s, diff, error


def getNormDispMdl(path):
    t = getTwiss(os.path.join(path, "getNDx.out"))
    s = t.S
    yvaluex1 = getattr(t, "NDX")
    yvaluex2 = getattr(t, "NDXMDL")
    error = getattr(t, "STDNDX")
    return s, yvaluex1, yvaluex2, error


def getNormDispDiff(path):
    t = getTwiss(os.path.join(path, "getNDx.out"))
    s = t.S
    diff = getattr(t, "NDX") - getattr(t, "NDXMDL")
    error = getattr(t, "STDNDX")
    return s, diff, error


def getChromaticAmp(path, plane):
    t = getTwiss(os.path.join(path, "chrombeta" + plane + ".out"))
    s = t.S
    w = getattr(t, "W" + plane.upper())
    wm = getattr(t, "W" + plane.upper() + "M")
    werr = getattr(t, "W" + plane.upper() + "ERR")
    return s, w, wm, werr


#TODO: Test! could not test because there is no chromatic coupling plot in the GUI
def getChromaticCoup(path, plane, subnode):
    if (not os.path.exists(os.path.join(path, "chromcoupling_free.out"))):
        t = getTwiss(os.path.join(path, "chromcoupling.out"))
    else:
        t = getTwiss(os.path.join(path, "chromcoupling_free.out"))
    s = t.S
    if (subnode == "ChromaticCouplingReal"):
        if(plane == "x"):
            f = t.Cf1001r
            fstd = t.Cf1001rERR
        if(plane == "y"):
            f = t.Cf1010r
            fstd = t.Cf1010rERR
    elif (subnode == "ChromaticCouplingImaginary"):
        if(plane == "x"):
            f = t.Cf1001i
            fstd = t.Cf1001iERR
        if(plane == "y"):
            f = t.Cf1010i
            fstd = t.Cf1010iERR
    elif (subnode =='ChromaticCouplingAmp'):
        amp = []
        #for i in range(0,len(t.Cf1001i))
            
    return s, f, fstd


def plotAnalysis(path, label, accel, subnode, mainnode, minx, maxx, miny, maxy, hminx, hmaxx, hminy, hmaxy, legendx, legendy, legendh):
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    plx = plt.subplot(gs[0])
    
    paths = splitFiles(path)
    sx, valuex, xerr = [], [], []
    sy, valuey, yerr = [], [], []
    for i in range(len(paths)):
        if(subnode == "Beta_beat"):
            getsx, getvaluex, getxerr = getbeta(paths[i], "x", mainnode)
            getsy, getvaluey, getyerr = getbeta(paths[i], "y", mainnode)
        if(subnode == "Diff_Phase_PhMdl"):
            getsx, getvaluex, getxerr = getphasediff(paths[i], "x", mainnode)
            getsy, getvaluey, getyerr = getphasediff(paths[i], "y", mainnode)
        if(mainnode == "Coupling"):
            getsx, getvaluex, getxerr = getcouple(paths[i], "x", subnode)
            getsy, getvaluey, getyerr = getcouple(paths[i], "y", subnode)
        if(subnode == "CO"):
            getsx, getvaluex, getxerr = getclosedorbit(paths[i], "x")
            getsy, getvaluey, getyerr = getclosedorbit(paths[i], "y")
        if(subnode == "diff_Disp_DMdl"):
            getsx, getvaluex, getxerr = getdisperssiondiff(paths[i], "x")
            getsy, getvaluey, getyerr = getdisperssiondiff(paths[i], "y")
        if(mainnode == "ChromaticCoupling"):
            getsx, getvaluex, getxerr = getChromaticCoup(paths[i], "x", subnode)
            getsy, getvaluey, getyerr = getChromaticCoup(paths[i], "y", subnode)
        if(subnode == "diff_NDisp_NDMdl"):
            getsx, getvaluex, getxerr = getNormDispDiff(paths[i])
        sx.append(getsx)
        valuex.append(getvaluex)
        xerr.append(getxerr)
        setLimits(accel, minx, maxx, hminx, hmaxx, plx)
        if label == 'None':
            labels = paths
        else:
            labels = label.split(',')
        plx.errorbar(sx[i], valuex[i], yerr=xerr[i], fmt='o', color=colors[i],  markersize=4, markeredgecolor=markeredgecolors[i], label=labels[i].rsplit('/', 1)[-1])
        if(mainnode != "Normalized_Dispersion"):
            plx.axes.get_xaxis().set_visible(False)
            ply = plt.subplot(gs[1])
            sy.append(getsy)
            valuey.append(getvaluey)
            yerr.append(getyerr)
            setLimits(accel, miny, maxy, hminy, hmaxy, ply)
            ply.errorbar(sy[i], valuey[i], yerr=yerr[i], fmt='o', color=colors[i],  markersize=4, markeredgecolor=markeredgecolors[i], label=labels[i].rsplit('/', 1)[-1])
            setYAxisLabel(subnode, 'y', ply)
            showIRs(accel, float(maxy), [ply])
            setXAxisLabel(ply)
        elif(mainnode == "Normalized_Dispersion"):
            plx.axes.get_xaxis().set_visible(True)
            plx.axes.get_xaxis().set_visible(True)
    plt.grid(False)
    setYAxisLabel(subnode, 'x', plx)
    if(int(float(legendh)) > 12):
        showLegend(plx, int(float(legendx)), int(float(legendy)))
    return gs


def plotMdlAnalysis(path, label, accel, subnode, mainnode, minx, maxx, miny, maxy, hminx, hmaxx, hminy, hmaxy,legendx, legendy, legendh):
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    plx = plt.subplot(gs[0])
    setLimits(accel, minx, maxx, hminx, hmaxx, plx)
    if(subnode == "Beta_BMdl"):
        sx, measx, mdlx, errx = getbetamdl(path, "x", mainnode)
        sy, measy, mdly, erry = getbetamdl(path, "y", mainnode)
    if(subnode == "Phase_PhMdl"):
        sx, measx, mdlx, errx = getphasemdl(path, "x", mainnode)
        sy, measy, mdly, erry = getphasemdl(path, "y", mainnode)
    if(subnode == "Disp_DMdl"):
        sx, measx, mdlx, errx = getdisperssionmdl(path, "x")
        sy, measy, mdly, erry = getdisperssionmdl(path, "y")
    if(subnode == "ChromaticAmplitude"):
        sx, measx, mdlx, errx = getChromaticAmp(path, "x")
        sy, measy, mdly, erry = getChromaticAmp(path, "y")
    if(subnode == "NDisp_NDMdl"):
        sx, measx, mdlx, errx = getNormDispMdl(path)
    if label == 'None':
        labels = ["mo_" + path.rsplit('/', 1)[-1], "me_" + path.rsplit('/', 1)[-1] ]
    else:
        labels = label.split(',')
    plx.errorbar(sx, mdlx, yerr=errx, fmt='o',color=colors[1], markersize=4, markeredgecolor= markeredgecolors[1], label=labels[0])
    plx.errorbar(sx, measx, yerr=errx, fmt='o',color=colors[0], markersize=4, markeredgecolor= markeredgecolors[0], label=labels[1])
    #plx.tick_params(labelsize=6)
    if("NDisp" not in subnode):
        plx.axes.get_xaxis().set_visible(False)
        ply = plt.subplot(gs[1])
        setLimits(accel, miny, maxy, hminy, hmaxy, ply)
        ply.errorbar(sy, mdly, yerr=erry, fmt='o',color=colors[1], markersize=4, markeredgecolor= markeredgecolors[1], label=labels[0])
        ply.errorbar(sy, measy, yerr=erry, fmt='o',color=colors[0], markersize=4, markeredgecolor= markeredgecolors[0], label=labels[1])
        #ply.tick_params(labelsize=6)
        setYAxisLabel(subnode, 'y', ply)
        setXAxisLabel(ply)
        showIRs(accel, float(maxy), [ply])
    else:
        setXAxisLabel(plx)
        plx.axes.get_xaxis().set_visible(True)
    plt.grid(False)
    setYAxisLabel(subnode, 'x', plx)
    if(int(float(legendh)) > 12):
        showLegend(plx, int(float("53.0")), int(float("7.0")))
    return gs


#################################################################################


#################   LABELS, LEGEND AND TEXT   ###################################

def setYAxisLabel(subnode, axis, p1):
    if (subnode == 'Phase_PhMdl'):
        p1.set_ylabel(r'$\phi_' + axis + '[2*' + r'\pi' + ']$')
    if (subnode == 'Diff_Phase_PhMdl'):
        p1.set_ylabel(r'$\Delta \phi_' + axis + '[2**' + r'\pi]$')
    if (subnode == 'Beta_beat'):
        p1.set_ylabel(r'$\Delta \beta_' + axis + r'/ \beta_' + axis +'$')
    if (subnode == 'Beta_BMdl'):
        p1.set_ylabel(r'$\beta_' + axis + '  [m]$')
    if (subnode == 'amp'):
        if(axis == 'x'):
            p1.set_ylabel(r'abs(F1001)')
        if(axis == 'y'):
            p1.set_ylabel(r'F1010W')
    if (subnode == 'real'):
        if(axis == 'x'):
            p1.set_ylabel(r're(F1001R)')
        if(axis == 'y'):
            p1.set_ylabel(r'F1010R')
    if (subnode == 'imaginary'):
        if(axis == 'x'):
            p1.set_ylabel(r'im(F1001)')
        if(axis == 'y'):
            p1.set_ylabel(r'F1010I')
    if (subnode == 'Disp_DMdl'):
        p1.set_ylabel(r'D' + axis + ' [m]')
    if (subnode == 'diff_Disp_DMdl'):
        p1.set_ylabel(r'$\Delta D / \beta [m]$')
    if (subnode == 'NDisp_NDMdl'):
        p1.set_ylabel(r'$ND' + axis + ' [sqrt(m)]$')
    if (subnode == 'diff_NDisp_NDMdl'):
        p1.set_ylabel(r'$\frac{\Delta D_x}{\beta_vx}  [m]$')
    if (subnode == 'CO'):
        p1.set_ylabel(r'$\Delta' + axis + ' [m]$')
    if (subnode == 'ChromaticAmplitude'):
        p1.set_ylabel(r'W' + axis)
    if (subnode == 'ChromaticCouplingReal'):
        if(axis == 'x'):
            p1.set_ylabel(r'$\Delta ' + 'Real(f1001)/' + '\Delta' + '\delta$')
        if(axis == 'y'):
            p1.set_ylabel(r'$\Delta ' + 'Real(f1010)/' + '\Delta' + '\delta$')
    if (subnode == 'ChromaticCouplingImaginary'):
        if(axis == 'x'):
            p1.set_ylabel(r'$\Delta ' + 'Imaginary(f1001)/' + '\Delta' + '\delta$')
        if(axis == 'y'):
            p1.set_ylabel(r'$\Delta ' + 'Imaginary(f1010)/' + '\Delta' + '\delta$')


def setXAxisLabel(p1):
    p1.set_xlabel(r'Longitudinal location [m]')


def showLegend(p, legendx, legendy, frameon=False, numpoints=1, ncol=1):
    handles, labels = p.get_legend_handles_labels()
    #p.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.02, 1.35), frameon=frameon, numpoints=numpoints, ncol=ncol)
    p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),
          fancybox=True, shadow=True, ncol=3)

def showIRs(accel, maxy, plots=[]):
    IRposition = maxy * 1.25   
    if accel[0:4] == "LHCB":
        for i in plots:
            theSizeofFont = 14
            if accel == "LHCB1":
                i.text(22800, IRposition, 'IP1', fontsize=theSizeofFont)
                i.text(0, IRposition, 'IP2', fontsize=theSizeofFont)
                i.text(3000, IRposition, 'IP3', fontsize=theSizeofFont)
                i.text(6000, IRposition, 'IP4', fontsize=theSizeofFont)
                i.text(9100, IRposition, 'IP5', fontsize=theSizeofFont)
                i.text(12500, IRposition, 'IP6', fontsize=theSizeofFont)
                i.text(15800, IRposition, 'IP7', fontsize=theSizeofFont)
                i.text(19600, IRposition, 'IP8', fontsize=theSizeofFont)
            else:
                i.text(2700, IRposition, 'IP1', fontsize=theSizeofFont)
                i.text(6000, IRposition, 'IP2', fontsize=theSizeofFont)
                i.text(9300, IRposition, 'IP3', fontsize=theSizeofFont)
                i.text(12700, IRposition, 'IP4', fontsize=theSizeofFont)
                i.text(16000, IRposition, 'IP5', fontsize=theSizeofFont)
                i.text(19300, IRposition, 'IP6', fontsize=theSizeofFont)
                i.text(22700, IRposition, 'IP7', fontsize=theSizeofFont)
                i.text(26000, IRposition, 'IP8', fontsize=theSizeofFont)


def setTicksNoLabels(plots=[]):
    #plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
    for i in plots:
        i.tick_params(axis='x', which='both', labelbottom='off')