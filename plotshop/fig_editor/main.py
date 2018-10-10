"""
Module plotshop.fig_editor.main: Figure Editor Command Line
-------------------------------------------------------------------

Commandline launcher/wrapper for the advanced qt-gui for matpoltlib.

Currently starts a testing figure.
"""
import sys
from PyQt5 import QtWidgets
from main_window import MainWindow


def _quick_figure():
    """ Create a quick testing figure """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    plt.text(.5, .5, "Hallo")
    plt.errorbar(range(3), range(3), yerr=range(3), xerr=range(3), capsize=3)
    plt.plot([0, 3], [0, 5], color="red")
    return fig


# Main #######################################################################


def main(fig=None):
    """ Prepares QT-environment and starts the editor. """
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("fusion")
    form = MainWindow(fig)
    form.show()
    app.exec_()


# Script Mode #################################################################


if __name__ == "__main__":
    main(_quick_figure())
