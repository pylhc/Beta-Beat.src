from __future__ import unicode_literals
import __init__  # @UnusedImport
from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
import os
from Python_Classes4MAD import metaclass
from Utilities import tfs_pandas
from Utilities.plotting import plot_style as ps

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

colors = [SkyBlue1, ScarletRed1, Butter1, Aluminium1]
markeredgecolors = [SkyBlue3, ScarletRed3, Butter3,  Aluminium3]

IR_POS = {
    "LHCB1": {
        'IP1': 22800,
        'IP2': 0,
        'IP3': 3000,
        'IP4': 6000,
        'IP5': 9100,
        'IP6': 12500,
        'IP7': 15800,
        'IP8': 19600,
    },
    "LHCB2": {
        'IP1': 2700,
        'IP2': 6000,
        'IP3': 9300,
        'IP4': 12700,
        'IP5': 16000,
        'IP6': 19300,
        'IP7': 22700,
        'IP8': 26000,
    }
}

MANUAL_STYLE = {
    u'font.size': 15,
    u'legend.fontsize': 16,
    u'font.weight': u'normal',
    u'axes.labelweight': u'normal',
    u"font.family": 'sans-serif',
    u"font.serif": ['Computer Modern'],
    u'axes.grid': False,
    u'lines.marker': u'o',
    u'lines.markersize': 5.0,
    u'lines.linestyle': u'',
}

IP_FONT = {
    'size': 14,
    'weight': 'bold'
}


def splitFiles(files):
    return files.split(',')


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
        if(plane == "x"):
            f =[]
            fstd = t.Cf1001iERR
            for i in range(len(t.Cf1001r)):
                f.append(np.sqrt(t.Cf1001i[i]**2 +t.Cf1001r[i]**2))
        
        if(plane == "y"):
            f =[]
            fstd = t.Cf1010iERR
            for i in range(len(t.Cf1010r)):
                f.append(np.sqrt(t.Cf1010i[i]**2 +t.Cf1010r[i]**2))
        
        #for i in range(0,len(t.Cf1001i))
            
    return s, f, fstd


def plotAnalysis(opt):
    ps.set_style("standard", MANUAL_STYLE)
    xmin = min(opt.xplot_xmin, opt.yplot_xmin)
    xmax = max(opt.xplot_xmax, opt.yplot_xmax)

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    plx = plt.subplot(gs[0])
    
    paths = splitFiles(opt.path)
    sx, valuex, xerr = [], [], []
    sy, valuey, yerr = [], [], []
    for idx, path in enumerate(paths):
        if opt.subnode == "Beta_beat":
            getsx, getvaluex, getxerr = getbeta(path, "x", opt.mainnode)
            getsy, getvaluey, getyerr = getbeta(path, "y", opt.mainnode)
        if opt.subnode == "Diff_Phase_PhMdl":
            getsx, getvaluex, getxerr = getphasediff(path, "x", opt.mainnode)
            getsy, getvaluey, getyerr = getphasediff(path, "y", opt.mainnode)
        if opt.mainnode == "Coupling":
            getsx, getvaluex, getxerr = getcouple(path, "x", opt.subnode)
            getsy, getvaluey, getyerr = getcouple(path, "y", opt.subnode)
        if opt.subnode == "CO":
            getsx, getvaluex, getxerr = getclosedorbit(path, "x")
            getsy, getvaluey, getyerr = getclosedorbit(path, "y")
        if opt.subnode == "diff_Disp_DMdl":
            getsx, getvaluex, getxerr = getdisperssiondiff(path, "x")
            getsy, getvaluey, getyerr = getdisperssiondiff(path, "y")
        if opt.mainnode == "ChromaticCoupling":
            getsx, getvaluex, getxerr = getChromaticCoup(path, "x", opt.subnode)
            getsy, getvaluey, getyerr = getChromaticCoup(path, "y", opt.subnode)
        if opt.subnode == "diff_NDisp_NDMdl":
            getsx, getvaluex, getxerr = getNormDispDiff(path)
        sx.append(getsx)
        valuex.append(getvaluex)
        xerr.append(getxerr)
        plx.set_xlim(xmin, xmax)
        plx.set_ylim(opt.xplot_ymin, opt.xplot_ymax)

        if opt.label == 'None':
            labels = paths
        else:
            labels = opt.label.split(',')

        plx.errorbar(sx[idx], valuex[idx], yerr=xerr[idx], fmt=rcParams['lines.marker'],
                     color=colors[idx], markeredgecolor=markeredgecolors[idx],
                     label=labels[idx].rsplit('/', 1)[-1])

        if opt.mainnode != "Normalized_Dispersion":
            plx.axes.get_xaxis().set_visible(False)
            ply = plt.subplot(gs[1])
            sy.append(getsy)
            valuey.append(getvaluey)
            yerr.append(getyerr)
            ply.set_xlim(xmin, xmax)
            ply.set_ylim(opt.yplot_ymin, opt.yplot_ymax)
            ply.errorbar(sy[idx], valuey[idx], yerr=yerr[idx], fmt=rcParams['lines.marker'],
                         color=colors[idx], markeredgecolor=markeredgecolors[idx],
                         label=labels[idx].rsplit('/', 1)[-1])
            set_yaxis_label(ply, 'y', opt.subnode)
            show_irs(ply, opt.accel, opt.yplot_ymax, xmin, xmax)
            ps.set_xaxis_label(ply)
        elif opt.mainnode == "Normalized_Dispersion":
            plx.axes.get_xaxis().set_visible(True)
            plx.axes.get_xaxis().set_visible(True)

    set_yaxis_label(plx, 'x', opt.subnode)
    if int(opt.legendh) > 12:
        show_legend(plx, int(opt.legendx), int(opt.legendy))
    return gs


def plotMdlAnalysis(opt):
    ps.set_style("standard", MANUAL_STYLE)
    xmin = min(opt.xplot_xmin, opt.yplot_xmin)
    xmax = max(opt.xplot_xmax, opt.yplot_xmax)

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    plx = plt.subplot(gs[0])
    plx.set_xlim(xmin, xmax)
    plx.set_ylim(opt.xplot_ymin, opt.xplot_ymax)
    if opt.subnode == "Beta_BMdl":
        sx, measx, mdlx, errx = getbetamdl(opt.path, "x", opt.mainnode)
        sy, measy, mdly, erry = getbetamdl(opt.path, "y", opt.mainnode)
    if opt.subnode == "Phase_PhMdl":
        sx, measx, mdlx, errx = getphasemdl(opt.path, "x", opt.mainnode)
        sy, measy, mdly, erry = getphasemdl(opt.path, "y", opt.mainnode)
    if opt.subnode == "Disp_DMdl":
        sx, measx, mdlx, errx = getdisperssionmdl(opt.path, "x")
        sy, measy, mdly, erry = getdisperssionmdl(opt.path, "y")
    if opt.subnode == "ChromaticAmplitude":
        sx, measx, mdlx, errx = getChromaticAmp(opt.path, "x")
        sy, measy, mdly, erry = getChromaticAmp(opt.path, "y")
    if opt.subnode == "NDisp_NDMdl":
        sx, measx, mdlx, errx = getNormDispMdl(opt.path)
    if opt.label == 'None':
        labels = ["mo_" + opt.path.rsplit('/', 1)[-1], "me_" + opt.path.rsplit('/', 1)[-1] ]
    else:
        labels = opt.label.split(',')
    if "ChromaticAmplitude" in opt.subnode:
        plx.plot(sx, mdlx,
                 color=colors[1], markeredgecolor=markeredgecolors[1],
                 linewidth=2, label=labels[0])
        plx.errorbar(sx, measx, yerr=errx, fmt=rcParams['lines.marker'],
                     color=colors[0], markeredgecolor=markeredgecolors[0],
                     label=labels[1])
    else:
        plx.errorbar(sx, mdlx, yerr=errx, fmt=rcParams['lines.marker'],
                     color=colors[1], markeredgecolor=markeredgecolors[1],
                     label=labels[0])
        plx.errorbar(sx, measx, yerr=errx, fmt=rcParams['lines.marker'],
                     color=colors[0], markeredgecolor=markeredgecolors[0],
                     label=labels[1])
    #plx.tick_params(labelsize=6)

    if "NDisp" not in opt.subnode:
        plx.axes.get_xaxis().set_visible(False)
        ply = plt.subplot(gs[1])
        ply.set_xlim(xmin, xmax)
        ply.set_ylim(opt.yplot_ymin, opt.yplot_ymax)
        if "ChromaticAmplitude" in opt.subnode:
            ply.plot(sy, mdly,
                     color=colors[1], markeredgecolor=markeredgecolors[1],
                     label=labels[0], linewidth=2,)
            ply.errorbar(sy, measy, yerr=erry, fmt=rcParams['lines.marker'],
                         color=colors[0], markeredgecolor=markeredgecolors[0],
                         label=labels[1])
        else:
            ply.errorbar(sy, mdly, yerr=erry, fmt=rcParams['lines.marker'],
                         color=colors[1], markeredgecolor=markeredgecolors[1],
                         label=labels[0])
            ply.errorbar(sy, measy, yerr=erry, fmt=rcParams['lines.marker'],
                         color=colors[0], markeredgecolor=markeredgecolors[0],
                         label=labels[1])
        set_yaxis_label(ply, 'y', opt.subnode)
        ps.set_xaxis_label(ply)
        show_irs(ply, opt.accel, opt.yplot_ymax, xmin, xmax)
    else:
        ps.set_xaxis_label(plx)
        plx.axes.get_xaxis().set_visible(True)

    set_yaxis_label(plx, 'x', opt.subnode)
    if(int(opt.legendh) > 12):
        show_legend(plx, int(53.0), int(7.0))
    return gs


#################################################################################


#################   LABELS, LEGEND AND TEXT   ###################################

def set_yaxis_label(plot, axis, subnode):
    labels_map = {
        'Phase_PhMdl': 'phase',
        'Diff_Phase_PhMdl': 'phase',
        'Beta_beat': 'betabeat',
        'Beta_BMdl': 'beta',
        'Disp_DMdl': 'dispersion',
        'diff_Disp_DMdl': 'dispersion',
        'NDisp_NDMdl': 'norm_dispersion',
        'diff_NDisp_NDMdl': 'norm_dispersion',
        'CO': 'co',
        'ChromaticAmplitude': 'chromamp',
        'ChromaticCouplingReal': 'real',
        'ChromaticCouplingImaginary': 'imag',
        'ChromaticCouplingAmp': 'absolute',
        'amp': 'absolute',
        'real': 'real',
        'imaginary': 'imag',
    }

    axis_map = {'x': r'f_{{1001}}', 'y': r'f_{{1010}}'}

    delta = False
    chromcoup = False

    if subnode.lower().startswith('diff') or subnode.lower() == "co":
        delta = True

    if subnode.lower().startswith('ChromaticCoupling'):
        delta = True
        chromcoup = True
        axis = axis_map[axis]

    if subnode.lower() in ['amp', 'real', 'imaginary']:
        axis = axis_map[axis]

    ps.set_yaxis_label(param=labels_map[subnode],
                       plane=axis,
                       ax=plot,
                       delta=delta,
                       chromcoup=chromcoup,
                       )


def show_legend(p, legendx, legendy, frameon=False, numpoints=1, ncol=1):
    handles, labels = p.get_legend_handles_labels()
    #p.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.02, 1.35), frameon=frameon, numpoints=numpoints, ncol=ncol)
    p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),
             fancybox=True, shadow=True, ncol=3)


def show_irs(plot, accel, ymax, xmin, xmax):
    ypos = ymax * 1.05
    if accel in IR_POS:
        for ir in IR_POS[accel]:
            if xmin <= IR_POS[accel][ir] <= xmax:
                plot.text(IR_POS[accel][ir], ypos, ir, fontdict=IP_FONT)
