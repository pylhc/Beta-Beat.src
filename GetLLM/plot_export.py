#!/afs/cern.ch/work/o/omc/anaconda/bin/python
import optparse
import os
import sys


analysis = ['Diff_Phase_PhMdl', 'Beta_beat', 'amp', 'real', 'imaginary', 'diff_Disp_DMdl', 'diff_NDisp_NDMdl', 'CO', 'ChromaticCouplingReal', 'ChromaticCouplingImaginary']
mdlAnalysis = ['Phase_PhMdl', 'Beta_BMdl', 'Disp_DMdl', 'NDisp_NDMdl', 'ChromaticAmplitude']


def _parse_options():
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                     help="What accelerator: ESRF LHCB1 LHCB2",
                     metavar="ACCEL", default="LHCB1", dest="accel")
    parser.add_option("-p", "--path",
                     help="Path to plot data",
                     metavar="PATH", default="", dest="paths")
    parser.add_option("-f", "--folder",
                     help="Folder for pdf file",
                     metavar="FOLDER", default="./", dest="folder")
    parser.add_option("-t", "--title",
                     help="Title",
                     metavar="TITLE", default="", dest="title")
    parser.add_option("-x", "--maxx",
                     help="Upper limit x-plane chart",
                     metavar="MAXX", default="", dest="maxx")
    parser.add_option("-z", "--minx",
                     help="Lower limit x-plane chart",
                     metavar="MINX", default="", dest="minx")
    parser.add_option("-y", "--maxy",
                     help="Upper limit y-plane chart",
                     metavar="MAXY", default="", dest="maxy")
    parser.add_option("-q", "--miny",
                     help="Lower limit y-plane chart",
                     metavar="MINY", default="", dest="miny")
    parser.add_option("-u", "--hmaxx",
                     help="Upper limit x-plane chart x-axis",
                     metavar="HMAXX", default="", dest="hmaxx")
    parser.add_option("-o", "--hminx",
                     help="Lower limit x-plane chart x-axis",
                     metavar="HMINX", default="", dest="hminx")
    parser.add_option("-b", "--hmaxy",
                     help="Upper limit y-plane chart x-axis",
                     metavar="HMAXY", default="", dest="hmaxy")
    parser.add_option("-c", "--hminy",
                     help="Lower limit y-plane chart x-axis",
                     metavar="HMINY", default="", dest="hminy")
    parser.add_option("-s", "--plot",
                     help="Selected plot",
                     metavar="PLOT", default="", dest="plot")
    parser.add_option("-n", "--mainnode",
                     help="Main node",
                     metavar="NODE", default="", dest="mainnode")
    options, _ = parser.parse_args()
    return options


def main():
    options = _parse_options()

    accel = options.accel
    paths = options.paths
    folder = options.folder
    plot = options.plot
    maxx = options.maxx
    minx = options.minx
    maxy = options.maxy
    miny = options.miny
    hmaxx = options.hmaxx
    hminx = options.hminx
    hmaxy = options.hmaxy
    hminy = options.hminy
    mainnode = options.mainnode
    title = options.title

    SimplePlotting.setParams()
    mpdf = PdfPages(os.path.abspath(folder))
    if(plot in mdlAnalysis):
        SimplePlotting.plotMdlAnalysis(paths, accel, plot, mainnode, minx, maxx, miny, maxy, hminx, hmaxx, hminy, hmaxy, title)
    elif(plot in analysis):
        SimplePlotting.plotAnalysis(paths, accel, plot, mainnode, minx, maxx, miny, maxy, hminx, hmaxx, hminy, hmaxy, title)
    else:
        print "Unknown subnode"
    plt.tight_layout
    try:
        mpdf.savefig()
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
