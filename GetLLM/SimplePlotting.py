import __init__  # @UnusedImport
from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from  Python_Classes4MAD import metaclass

colors = ['red', 'blue', 'green', 'yellow']


def setParams(labelfontsize=8, legendfontsize=8):
    rc('font', **{'family': 'sans-serif', 'serif': ['Computer Modern']})
    params = {'backend': 'pdf',
          'axes.labelsize': labelfontsize,
          'font.size': labelfontsize,
          'axes.titlesize': labelfontsize,
          'legend.fontsize': 6,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
#          'axes.title': 12,
          'text.usetex': False,
          'axes.unicode_minus': True,
          'xtick.major.pad': 8,
          'ytick.major.pad': 8,
          'xtick.minor.pad': 8,
          'ytick.minor.pad': 8,
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
    return s, f, fstd


def plotAnalysis(path, accel, subnode, mainnode, minx, maxx, miny, maxy, hminx, hmaxx, hminy, hmaxy, title, legendx, legendy, legendh):
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    plx = plt.subplot(gs[0])
    ply = plt.subplot(gs[1])
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
        plx.plot(sx[i], valuex[i], color=colors[i], linestyle='solid', linewidth=0.01, marker='.', markerfacecolor=colors[i], markersize=4, label=paths[i].rsplit('/', 1)[-1])
        plx.errorbar(sx[i], valuex[i], yerr=xerr[i], linewidth=0.1, color=colors[i], marker='.', markersize=1, markeredgecolor=colors[i], label='_nolegend_')
        if(mainnode != "Normalized_Dispersion"):
            sy.append(getsy)
            valuey.append(getvaluey)
            yerr.append(getyerr)
            setLimits(accel, miny, maxy, hminy, hmaxy, ply)
            ply.plot(sy[i], valuey[i], color=colors[i], linestyle='solid', linewidth=0.01, marker='.', markerfacecolor=colors[i], markersize=4, label=paths[i].rsplit('/', 1)[-1])
            ply.errorbar(sy[i], valuey[i], yerr=yerr[i], linewidth=0.1, color=colors[i], marker='.', markersize=1, markeredgecolor=colors[i], label='_nolegend_')
            setYAxisLabel(subnode, 'y', ply)
    plt.suptitle(title)
    plx.axes.get_xaxis().set_visible(False)
    setXAxisLabel(ply)
    showIRs(accel, float(maxx), [plx])
    plt.grid(False)
    setYAxisLabel(subnode, 'x', plx)
    if(int(float(legendh)) > 12):
        showLegend(plx, int(float(legendx)), int(float(legendy)))
    return gs


def plotMdlAnalysis(path, accel, subnode, mainnode, minx, maxx, miny, maxy, hminx, hmaxx, hminy, hmaxy, title, legendx, legendy, legendh):
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    plx = plt.subplot(gs[0])
    ply = plt.subplot(gs[1])
    setLimits(accel, minx, maxx, hminx, hmaxx, plx)
    setLimits(accel, miny, maxy, hminy, hmaxy, ply)
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
    plx.plot(sx, measx, color='red', linestyle='solid', linewidth=0.1, marker='.', markerfacecolor='red', markersize=4)
    plx.plot(sx, mdlx, color='blue', linestyle='solid', linewidth=0.1, marker='.', markerfacecolor='blue', markersize=4)
    plx.errorbar(sx, measx, yerr=errx, fmt='', color='red', linewidth=0.1, marker='.', markersize=1, markeredgecolor='red', label='_nolegend_')
    #plx.tick_params(labelsize=6)
    if("NDisp" not in subnode):
        ply.plot(sy, measy, color='red', linestyle='solid', linewidth=0.1, marker='.', markerfacecolor='red', markersize=4, label="me_" + path.rsplit('/', 1)[-1])
        ply.plot(sy, mdly, color='blue', linestyle='solid', linewidth=0.1, marker='.', markerfacecolor='blue', markersize=4, label="mo_" + path.rsplit('/', 1)[-1])
        ply.errorbar(sy, measy, yerr=erry, color='red', linewidth=0.1, marker='.', markersize=1, markeredgecolor='red', label='_nolegend_')
        #ply.tick_params(labelsize=6)
        setYAxisLabel(subnode, 'y', ply)
    plt.suptitle(title)
    plx.axes.get_xaxis().set_visible(False)
    setXAxisLabel(ply)
    showIRs(accel, float(maxx), [plx])
    plt.grid(False)
    setYAxisLabel(subnode, 'x', plx)
    if(int(float(legendh)) > 12):
        showLegend(plx, int(float("53.0")), int(float("7.0")))
    return gs


#################################################################################


#################   LABELS, LEGEND AND TEXT   ###################################

def setYAxisLabel(subnode, axis, p1):
    if (subnode == 'Phase_PhMdl'):
        p1.set_ylabel(u'\u03c6' + axis + '[2*' + u'\u03c0' + ']')
    if (subnode == 'Diff_Phase_PhMdl'):
        p1.set_ylabel(u'\u0394' + u'\u03c6' + axis + '[2**' + u'\u03c0]')
    if (subnode == 'Beta_beat'):
        p1.set_ylabel('(' + u'\u0394' + u'\u03b2' + axis + ')\\' + u'\u03b2' + axis)
    if (subnode == 'Beta_BMdl'):
        p1.set_ylabel(u'\u03b2' + axis + '  [m]')
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
        p1.set_ylabel(u'\u0394D' + axis + ' [m]')
    if (subnode == 'NDisp_NDMdl'):
        p1.set_ylabel(r'ND' + axis + ' [sqrt(m)]')
    if (subnode == 'diff_NDisp_NDMdl'):
        p1.set_ylabel(u'\u0394D' + axis + ' [m]')
    if (subnode == 'CO'):
        p1.set_ylabel(u'\u0394' + axis + ' [m]')
    if (subnode == 'ChromaticAmplitude'):
        p1.set_ylabel(r'W' + axis)
    if (subnode == 'ChromaticCouplingReal'):
        if(axis == 'x'):
            p1.set_ylabel(u'\u0394' + 'Re(f1001)/' + u'\u0394' + u'\u03b4')
        if(axis == 'y'):
            p1.set_ylabel(u'\u0394' + 'Re(f1010)/' + u'\u0394' + u'\u03b4')
    if (subnode == 'ChromaticCouplingImaginary'):
        if(axis == 'x'):
            p1.set_ylabel(u'\u0394' + 'Im(f1001)/' + u'\u0394' + u'\u03b4')
        if(axis == 'y'):
            p1.set_ylabel(u'\u0394' + 'Im(f1010)/' + u'\u0394' + u'\u03b4')


def setXAxisLabel(p1):
    p1.set_xlabel(r'Longitudinal location [m]')


def showLegend(p, legendx, legendy, frameon=True, numpoints=1, ncol=3):
    handles, labels = p.get_legend_handles_labels()
    if(legendx < 1000 & legendy >= 200):
        loc = 'lower left'
    elif(legendx < 1000):
        loc = 'upper left'
    elif(legendx >= 1000 & legendy >= 200):
        loc = 'lower right'
    elif (legendx >= 1000 & legendy < 200):
        loc = 'upper right'
    else:
        print "Unknown parameter"
    p.legend(handles, labels, loc=loc, frameon=frameon, numpoints=numpoints, ncol=ncol)


def showIRs(accel, limit, plots=[]):
    if accel[0:4] == "LHCB":
        for i in plots:
            if accel == "LHCB1":
                i.text(22800, limit * 1.05, 'IP1', fontsize=8)
                i.text(0, limit * 1.05, 'IP2', fontsize=8)
                i.text(3000, limit * 1.05, 'IP3', fontsize=8)
                i.text(6000, limit * 1.05, 'IP4', fontsize=8)
                i.text(9100, limit * 1.05, 'IP5', fontsize=8)
                i.text(12500, limit * 1.05, 'IP6', fontsize=8)
                i.text(15800, limit * 1.05, 'IP7', fontsize=8)
                i.text(19600, limit * 1.05, 'IP8', fontsize=8)
            else:
                i.text(2700, limit * 1.05, 'IP1', fontsize=8)
                i.text(6000, limit * 1.05, 'IP2', fontsize=8)
                i.text(9300, limit * 1.05, 'IP3', fontsize=8)
                i.text(12700, limit * 1.05, 'IP4', fontsize=8)
                i.text(16000, limit * 1.05, 'IP5', fontsize=8)
                i.text(19300, limit * 1.05, 'IP6', fontsize=8)
                i.text(22700, limit * 1.05, 'IP7', fontsize=8)
                i.text(26000, limit * 1.05, 'IP8', fontsize=8)


def setTicksNoLabels(plots=[]):
    #plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
    for i in plots:
        i.tick_params(axis='x', which='both', labelbottom='off')
