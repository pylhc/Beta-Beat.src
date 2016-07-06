import sys
from PyQt4 import QtGui
from sbs_gui_main import SbSGuiMainController


def main(match_path=None):
    app = QtGui.QApplication(sys.argv)
    main_controller = SbSGuiMainController()
    if match_path is None:
        match_path = main_controller.ask_for_match_path()
        if match_path is None:
            return
    main_controller.set_match_path(match_path)
    main_controller.show_view()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
