#!/afs/cern.ch/work/o/omc/anaconda/bin/python
import argparse
import os
import sys


analysis = ['Diff_Phase_PhMdl', 'Beta_beat', 'amp', 'real', 'imaginary', 'diff_Disp_DMdl', 'diff_NDisp_NDMdl', 'CO', 'ChromaticCouplingReal', 'ChromaticCouplingImaginary', 'ChromaticCouplingAmp']
mdlAnalysis = ['Phase_PhMdl', 'Beta_BMdl', 'Disp_DMdl', 'NDisp_NDMdl', 'ChromaticAmplitude']


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


def main(args=None):
    opt = _parse_options(args)

    # SimplePlotting.setParams()
    mpdf = PdfPages(os.path.abspath(opt.folder))
    if opt.subnode in mdlAnalysis:
        SimplePlotting.plotMdlAnalysis(opt)
    elif opt.subnode in analysis:
        SimplePlotting.plotAnalysis(opt)
    else:
        print "Unknown subnode"

    plt.tight_layout()

    try:
        mpdf.savefig(bbox_inches='tight')
    finally:
        mpdf.close()


def _run_with_anaconda():
    from subprocess import call
    if not sys.platform == "darwin":  # This is Mac
        if "win" in sys.platform:
            print "There is not Windows version of Anaconda in OMC. Aborting."
            return
    interpreter = os.path.join("/afs", "cern.ch", "work", "o", "omc",
                               "anaconda", "bin", "python")
    command = sys.argv
    command.insert(0, interpreter)
    call(command)


if __name__ == "__main__":
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import SimplePlotting
        main()
    except ImportError:
        print "Cannot use this version of Python, trying OMC Anaconda..."
        _run_with_anaconda()
