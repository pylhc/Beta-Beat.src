from __future__ import unicode_literals
import __init__

import argparse
import os
import re

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages

from utils import logging_tools
from utils import tfs_pandas as tfs
from utils.plotting import plot_style as ps

LOG = logging_tools.get_logger(__name__)

"""
=============================   Constants, Style and Arguments   =============================
"""

IR_POS_DEFAULT = {
    "LHCB1": {
        'IP1': 23519.36962,
        'IP2': 192.923,
        'IP3': 3525.207216,
        'IP4': 6857.491433,
        'IP5': 10189.77565,
        'IP6': 13522.21223,
        'IP7': 16854.64882,
        'IP8': 20175.8654,
    },
    "LHCB2": {
        'IP1': 3195.252584,
        'IP2': 6527.5368,
        'IP3': 9859.973384,
        'IP4': 13192.40997,
        'IP5': 16524.84655,
        'IP6': 19857.13077,
        'IP7': 23189.41498,
        'IP8': 26510.4792,
    }
}

MANUAL_STYLE = {
    # differences to the standard style
    u'legend.fontsize': 16,
    u'font.weight': u'normal',
    u'axes.labelweight': u'normal',
    u'axes.grid': False,
    u'lines.markersize': 5.0,
    u'lines.linestyle': u'',
}

mdl_subnodes = ['Phase_PhMdl', 'Beta_BMdl', 'Disp_DMdl', 'NDisp_NDMdl', 'ChromaticAmplitude']


def _parse_options(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--accel",
                        help="What accelerator: ESRF LHCB1 LHCB2",
                        metavar="ACCEL", required=True, dest="accel", type=str)
    parser.add_argument("-p", "--path",
                        help="Path to plot data",
                        metavar="PATH", required=True, dest="path", type=str)
    parser.add_argument("-f", "--folder",
                        help="Folder for pdf file",
                        metavar="FOLDER", default="./", dest="folder", type=str)
    parser.add_argument("-x", "--maxx",
                        help="Upper limit x-plane chart",
                        metavar="MAXX", required=True, dest="xplot_ymax", type=float)
    parser.add_argument("-z", "--minx",
                        help="Lower limit x-plane chart",
                        metavar="MINX", required=True, dest="xplot_ymin", type=float)
    parser.add_argument("-y", "--maxy",
                        help="Upper limit y-plane chart",
                        metavar="MAXY", required=True, dest="yplot_ymax", type=float)
    parser.add_argument("-q", "--miny",
                        help="Lower limit y-plane chart",
                        metavar="MINY", required=True, dest="yplot_ymin", type=float)
    parser.add_argument("-u", "--hmaxx",
                        help="Upper limit x-plane chart x-axis",
                        metavar="HMAXX", required=True, dest="xplot_xmax", type=float)
    parser.add_argument("-o", "--hminx",
                        help="Lower limit x-plane chart x-axis",
                        metavar="HMINX", required=True, dest="xplot_xmin", type=float)
    parser.add_argument("-b", "--hmaxy",
                        help="Upper limit y-plane chart x-axis",
                        metavar="HMAXY", required=True, dest="yplot_xmax", type=float)
    parser.add_argument("-c", "--hminy",
                        help="Lower limit y-plane chart x-axis",
                        metavar="HMINY", required=True, dest="yplot_xmin", type=float)
    parser.add_argument("-s", "--plot",
                        help="Selected plot",
                        metavar="PLOT", default="", dest="subnode", type=str)
    parser.add_argument("-n", "--mainnode",
                        help="Main node",
                        metavar="NODE", default="", dest="mainnode", type=str)
    parser.add_argument("-w", "--legendx",
                        help="X-position of legend",
                        metavar="LEGENDX", default=0, dest="legendx", type=float)
    parser.add_argument("-d", "--legendy",
                        help="Y-Posotion of legend",
                        metavar="LEGENDY", default=0, dest="legendy", type=float)
    parser.add_argument("-l", "--legendh",
                        help="Height of legend",
                        metavar="LEGENDH", default=0, dest="legendh", type=float)
    parser.add_argument("-t", "--labels",
                        help="Label",
                        metavar="LABEL", default="None", dest="label", type=str)

    options = parser.parse_args(args)
    return options


"""
=============================   Main Wrapper   =============================
"""


def pdf_export(args=None):
    """ Exports to pdf """
    opt = _parse_options(args)
    plot_analysis(opt)

    # PDF export
    pdf_path = os.path.abspath(opt.folder)
    mpdf = PdfPages(pdf_path)
    plt.tight_layout()

    try:
        mpdf.savefig(bbox_inches='tight')
        LOG.debug("Exported GetLLM results to PDF '{:s}'".format(pdf_path))
    finally:
        mpdf.close()


def main(args=None):
    """ Plots and shows the plot """
    opt = _parse_options(args)
    plot_analysis(opt)
    plt.show()

"""
=============================   Getter Functions   =============================
"""


def get_path(path, name, ext):
    """ Returns the path to _free if available, without otherwise"""
    path_free = os.path.join(path, name + "_free" + ext)
    path = os.path.join(path, name + ext)
    if os.path.exists(path_free):
        return path_free
    elif os.path.exists(path):
        return path
    else:
        raise IOError("'{f:s}' not found in {p:s}".format(f=name+ext, p=path))


def get_tfs(path, name, ext=".out"):
    """ Returns the appropriate tfs dataframe (free or not) """
    return tfs.read_tfs(get_path(path, name, ext), index="NAME")


def get_beta(path, plane, mainnode):
    p = plane.upper()
    if mainnode == 'Beta':
        data = get_tfs(path, "getbeta" + plane)
        data_out = data.loc[:, ["S"]]

        try:
            # beta beat error
            data_out.loc[:, "Error"] = (
                np.sqrt(data.loc[:, "STDBET" + p] ** 2 + data.loc[:, "ERRBET" + p] ** 2) /
                data.loc[:, "BET" + p + "MDL"])
        except KeyError:
            data_out.loc[:, "Error"] = data.loc[:, "ERRBET" + p] / data.loc[:, "BET" + p + "MDL"]

    elif mainnode == 'Beta_Amp':
        data = get_tfs(path, "getampbeta" + plane)
        data_out = data.loc[:, ["S"]]
        data_out.loc[:, "Error"] = data.loc[:, "BET" + p + "STD"] / data.loc[:, "BET" + p + "MDL"]
        # TODO: Plase fix this STDBETX vs BETXSTD confusion! (in GetLLM)

    data_out.loc[:, "Value"] = ((data.loc[:, "BET" + p] - data.loc[:, "BET" + p + "MDL"]) /
                                data.loc[:, "BET" + p + "MDL"])
    return data_out


def get_beta_mdl(path, plane, mainnode):
    p = plane.upper()
    if mainnode == 'Beta':
        data = get_tfs(path, "getbeta" + plane)
        data_out = data.loc[:, ["S"]]

        try:
            data_out.loc[:, "Error"] = np.sqrt(
                data.loc[:, "STDBET" + p] ** 2 + data.loc[:, "ERRBET" + p] ** 2)
        except KeyError:
            data_out.loc[:, "Error"] = data.loc[:, "ERRBET" + p]

    elif mainnode == 'Beta_Amp':
        data = get_tfs(path, "getampbeta" + plane)
        data_out = data.loc[:, ["S"]]
        data_out.loc[:, "Error"] = data.loc[:, "BET" + p + "STD"]

    data_out.loc[:, "Value"] = data.loc[:, "BET" + p]
    data_out.loc[:, "Model"] = data.loc[:, "BET" + p + "MDL"]
    return data_out


def get_coupling(path, plane, subnode):
    data = get_tfs(path, "getcouple")
    data_out = data.loc[:, ["S"]]

    plane_map = {'x': '1001', 'y': '1010'}
    error_map = {'x': '1', 'y': '2'}
    node_map = {'amp': 'W', 'real': 'R', 'imaginary': 'I'}

    data_out.loc[:, "Value"] = data.loc[:, "F" + plane_map[plane] + node_map[subnode]]
    data_out.loc[:, "Error"] = data.loc[:, "FWSTD" + error_map[plane]]
    #TODO: Does it make sense to always use FWSTD here?

    return data_out


def get_phase_diff(path, plane, mainnode):
    if mainnode == 'Phase_Total':
        data = get_tfs(path, "getphasetot" + plane)
    elif mainnode == 'Phase':
        data = get_tfs(path, "getphase" + plane)
    p = plane.upper()
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "PHASE" + p] - data.loc[:, "PH" + p + "MDL"]
    data_out.loc[:, "Error"] = data.loc[:, "STDPH" + p]

    return data_out


def get_phase_mdl(path, plane, mainnode):
    if mainnode == 'Phase_Total':
        data = get_tfs(path, "getphasetot" + plane)
    elif mainnode == 'Phase':
        data = get_tfs(path, "getphase" + plane)
    p = plane.upper()
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "PHASE" + p]
    data_out.loc[:, "Error"] = data.loc[:, "STDPH" + p]
    data_out.loc[:, "Model"] = data.loc[:, "PH" + p + "MDL"]
    return data_out


def get_closedorbit(path, plane):
    data = get_tfs(path, "getCO" + plane)
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, plane.upper()]
    data_out.loc[:, "Error"] = data.loc[:, "STD" + plane.upper()]
    return data_out


def get_disperssion_mdl(path, plane):
    data = get_tfs(path, "getD" + plane)
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "D" + plane.upper()]
    data_out.loc[:, "Error"] = data.loc[:, "STDD" + plane.upper()]
    data_out.loc[:, "Model"] = data.loc[:, "D" + plane.upper() + "MDL"]
    return data_out


def get_disperssion_diff(path, plane):
    data = get_tfs(path, "getD" + plane)
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "D" + plane.upper()] -\
                               data.loc[:, "D" + plane.upper() + "MDL"]
    data_out.loc[:, "Error"] = data.loc[:, "STDD" + plane.upper()]
    return data_out


def get_norm_disp_mdl(path, plane):
    if plane == "y":
        LOG.debug("Normal Dispersion Y not implemented. (Just FYI, Taken care of in code)")
        return None

    data = get_tfs(path, "getND" + plane)
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "ND" + plane.upper()]
    data_out.loc[:, "Error"] = data.loc[:, "STDND" + plane.upper()]
    data_out.loc[:, "Model"] = data.loc[:, "ND" + plane.upper() + "MDL"]
    return data_out


def get_norm_disp_diff(path, plane):
    if plane == "y":
        LOG.debug("Normal Dispersion Y not implemented. (Just FYI, Taken care of in code)")
        return None

    data = get_tfs(path, "getND" + plane)
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "ND" + plane.upper()] - \
                               data.loc[:, "ND" + plane.upper() + "MDL"]
    data_out.loc[:, "Error"] = data.loc[:, "STDND" + plane.upper()]
    return data_out


def get_chromatic_amp(path, plane):
    data = get_tfs(path, "chrombeta" + plane)
    data_out = data.loc[:, ["S"]]
    data_out.loc[:, "Value"] = data.loc[:, "W" + plane.upper()]
    data_out.loc[:, "Error"] = data.loc[:, "W" + plane.upper() + "ERR"]
    data_out.loc[:, "Model"] = data.loc[:, "W" + plane.upper() + "M"]
    return data_out


def get_chromatic_coup(path, plane, subnode):
    data = get_tfs(path, "chromcoupling")
    data_out = data.loc[:, ["S"]]

    sn = subnode.lower().replace("chromaticcoupling", "")

    plane_map = {'x': '1001', 'y': '1010'}
    node_map = {'real': 'r', 'imaginary': 'i'}

    if sn in node_map:
        name = "Cf" + plane_map[plane] + node_map[sn]
        data_out.loc[:, "Value"] = data.loc[:, name]
        data_out.loc[:, "Error"] = data.loc[:, name + "ERR"]
    else:
        # sn == 'amp'
        name_r = "Cf" + plane_map[plane] + node_map['real']
        name_i = "Cf" + plane_map[plane] + node_map['imaginary']
        data_out.loc[:, "Value"] = np.sqrt(
            np.square(data.loc[:, name_r]) + np.square(data.loc[:, name_i]))
        data_out.loc[:, "Error"] = np.sqrt(
            np.square(data.loc[:, name_r + "ERR"]) + np.square(data.loc[:, name_i + "ERR"]))
        # TODO: Is this amp-error correct? Before only the imaginary error was used.
    return data_out


"""
=============================   Getter Collector   =============================
"""


def get_data(path, mainnode, subnode):
    if subnode == "Beta_beat":
        data_x = get_beta(path, "x", mainnode)
        data_y = get_beta(path, "y", mainnode)

    elif subnode == "Diff_Phase_PhMdl":
        data_x = get_phase_diff(path, "x", mainnode)
        data_y = get_phase_diff(path, "y", mainnode)

    elif mainnode == "Coupling":
        data_x = get_coupling(path, "x", subnode)
        data_y = get_coupling(path, "y", subnode)

    elif subnode == "CO":
        data_x = get_closedorbit(path, "x")
        data_y = get_closedorbit(path, "y")

    elif subnode == "diff_Disp_DMdl":
        data_x = get_disperssion_diff(path, "x")
        data_y = get_disperssion_diff(path, "y")

    elif mainnode == "ChromaticCoupling":
        data_x = get_chromatic_coup(path, "x", subnode)
        data_y = get_chromatic_coup(path, "y", subnode)

    elif subnode == "diff_NDisp_NDMdl":
        data_x = get_norm_disp_diff(path, "x")
        data_y = get_norm_disp_mdl(path, "x")

    elif subnode == "Beta_BMdl":
        data_x = get_beta_mdl(path, "x", mainnode)
        data_y = get_beta_mdl(path, "y", mainnode)

    elif subnode == "Phase_PhMdl":
        data_x = get_phase_mdl(path, "x", mainnode)
        data_y = get_phase_mdl(path, "y", mainnode)

    elif subnode == "Disp_DMdl":
        data_x = get_disperssion_mdl(path, "x")
        data_y = get_disperssion_mdl(path, "y")

    elif subnode == "ChromaticAmplitude":
        data_x = get_chromatic_amp(path, "x")
        data_y = get_chromatic_amp(path, "y")

    elif subnode == "NDisp_NDMdl":
        data_x = get_norm_disp_mdl(path, "x")
        data_y = get_norm_disp_mdl(path, "y")

    else:
        raise IOError("Subnode '{:s}' unknown".format(subnode))

    return data_x, data_y


"""
=============================   Main Plotting functions   =============================
"""


def plot_analysis(opt):
    """ Plots the specified results

    For required opt-structure look in parse_options.
    """
    LOG.debug("Plotting GetLLM analysis.")
    mdl_analysis = opt.subnode in mdl_subnodes

    ps.set_style("standard", MANUAL_STYLE)
    xmin = min(opt.xplot_xmin, opt.yplot_xmin)
    xmax = max(opt.xplot_xmax, opt.yplot_xmax)

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    ax_x = plt.subplot(gs[0])
    ax_y = None
    ir_pos = None

    paths = opt.path.split(',')

    if opt.label == 'None':
        if mdl_analysis:
            labels = ["mo_" + opt.path.rsplit('/', 1)[-1], "me_" + opt.path.rsplit('/', 1)[-1]]
        else:
            labels = paths
    else:
        labels = opt.label.split(',')

    for idx, path in enumerate(paths):
        data_x, data_y = get_data(path, opt.mainnode, opt.subnode)
        plot_data(ax_x, data_x, labels, idx)

        if ir_pos is None:
            ir_pos = get_irpos(data_x, opt.accel)

        if data_y is not None:
            if ax_y is None:
                ax_x.axes.get_xaxis().set_visible(False)
                ax_y = plt.subplot(gs[1])
            plot_data(ax_y, data_y, labels, idx)

    ax_x.set_xlim(xmin, xmax)
    ax_x.set_ylim(opt.xplot_ymin, opt.xplot_ymax)
    set_yaxis_label(ax_x, 'x', opt.subnode)

    if ax_y is not None:
        ax_y.set_xlim(xmin, xmax)
        ax_y.set_ylim(opt.yplot_ymin, opt.yplot_ymax)
        set_yaxis_label(ax_y, 'y', opt.subnode)
        ps.set_xaxis_label(ax_y)
        if ir_pos:
            ps.show_ir(ir_pos, ax_y, mode='outside')
            ps.show_ir(ir_pos, ax_x, mode='lines')
    else:
        ax_x.axes.get_xaxis().set_visible(True)
        ps.set_xaxis_label(ax_x)
        if ir_pos:
            ps.show_ir(ir_pos, ax_x, mode='outside')

    if int(opt.legendh) > 12:
        show_legend(ax_x, int(opt.legendx), int(opt.legendy))
    return gs


def plot_data(ax, data, labels, idx):
    """ Actual plotting """
    if "Model" in data.columns:
        if idx > 0:
            raise NotImplementedError("Only single model comparison implemented.")

        ax.plot(data["S"], data["Model"],
                color=get_color(idx), markeredgecolor=get_color(idx, True), label=labels[idx])
        idx = idx + 1
        label = labels[idx]
    else:
        label = labels[idx].rsplit('/', 1)[-1]  # TODO: Works only on linux

    ax.errorbar(data["S"], data["Value"], yerr=data["Error"],
                fmt=rcParams['lines.marker'],
                color=get_color(idx), markeredgecolor=get_color(idx, True),
                label=label)


"""
=============================   Labels, Legend, Text   =============================
"""


def set_yaxis_label(plot, axis, subnode):
    labels_map = {  # mapping for labels known by plot_style.py
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

    coupling_map = {'x': r'f_{{1001}}', 'y': r'f_{{1010}}'}

    delta = False
    chromcoup = False

    if subnode.lower().startswith('diff') or subnode.lower() == "co":
        delta = True

    if subnode.lower().startswith('chromaticcoupling'):
        delta = True
        chromcoup = True
        axis = coupling_map[axis]

    if subnode.lower() in ['amp', 'real', 'imaginary']:
        axis = coupling_map[axis]

    ps.set_yaxis_label(param=labels_map[subnode],
                       plane=axis,
                       ax=plot,
                       delta=delta,
                       chromcoup=chromcoup,
                       )


def show_legend(p, legendx, legendy, frameon=False, numpoints=1, ncol=1):
    # handles, labels = p.get_legend_handles_labels()
    # p.legend(handles, labels,
    #          loc='upper left', bbox_to_anchor=(0.02, 1.35),
    #          frameon=frameon, numpoints=numpoints, ncol=ncol)

    p.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),
             fancybox=True, shadow=True, ncol=3)


def get_color(idx, marker=False):
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
    markeredgecolors = [SkyBlue3, ScarletRed3, Butter3, Aluminium3]

    if not marker:
        return colors[idx % len(colors)]
    else:
        return markeredgecolors[idx % len(colors)]


def get_irpos(data, accel):
    # try to load from model_elements in data
    try:
        model_path = re.search(r"(?<=--model=)[^']+", data.Command)  # DoNormal
        if not model_path:
            model_path = re.search(r"(?<=--twissfile=)[^']+", data.Command)  # DoWAnalysis
        model_path = model_path.group()
    except (KeyError, AttributeError):
        pass
    else:
        f, e = os.path.splitext(model_path)
        elements_path = f + '_elements' + e
        try:
            ir_pos = ps.get_ip_positions(elements_path)
        except (IOError, KeyError):
            pass
        else:
            LOG.debug("Loaded IP positions from '{:s}'".format(elements_path))
            return ir_pos

    # loading failed, use defaults
    if accel in IR_POS_DEFAULT:
        LOG.debug("Loaded IP positions from default values.")
        return IR_POS_DEFAULT[accel]
    else:
        LOG.debug("Could not find appropriate IP positions.")
        return {}


if __name__ == "__main__":
    pdf_export()
