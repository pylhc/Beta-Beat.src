"""
Launcher for the advanced qt-gui for matpoltlib.
"""
import sys
from PyQt5 import QtWidgets
from main_window import MainWindow


# Main #######################################################################


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("fusion")
    form = MainWindow()
    form.show()
    app.exec_()


# Script Mode #################################################################


if __name__ == "__main__":
    main()
