import os
from matcher_model_default import MatcherModelDefault
from matcher_model_default import MatcherPlotterDefault
from matchers.phase_matcher import PhaseMatcher
from Utilities import tfs_pandas


class MatcherModelPhase(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = super(MatcherModelPhase, self).get_matcher_dict()
        matcher_dict["type"] = "phase"
        return matcher_dict

    def create_matcher(self, lhc_mode, match_path):
        self._matcher = PhaseMatcher(
            lhc_mode,
            self._beam,
            self._name,
            self.get_matcher_dict(),
            match_path
        )

    def get_plotter(self, figure):
        if self._plotter is None:
            self._plotter = MatcherPlotterPhase(figure, self)
        return self._plotter


class MatcherPlotterPhase(MatcherPlotterDefault):

    def plot(self):
        self._figure.clear()
        self._axes_to_data = {}
        axes_x = self._figure.add_subplot(2, 1, 1)
        axes_x.set_title(self._matcher_model.get_name())
        file_horizontal = tfs_pandas.read_tfs(os.path.join(
            self._matcher_model.get_output_path(), "sbs",
            "sbsphasext_" + str(self._matcher_model.get_label()) + ".out")
        )
        self._axes_to_data[axes_x] = file_horizontal
        self._plot_match(axes_x, file_horizontal, "X")
        file_vertical = tfs_pandas.read_tfs(os.path.join(
            self._matcher_model.get_output_path(), "sbs",
            "sbsphaseyt_" + str(self._matcher_model.get_label()) + ".out")
        )
        axes_y = self._figure.add_subplot(2, 1, 2)
        self._axes_to_data[axes_y] = file_vertical
        self._plot_match(axes_y, file_vertical, "Y")
        self._figure.tight_layout()
        self._figure.patch.set_visible(False)
        self._figure.canvas.draw()

    def _plot_match(self, axes, sbs_file, plane):
        if self._matcher_model.get_propagation() == "front":
            MatcherPlotterPhase._plot_front(axes, sbs_file, plane)
        elif self._matcher_model.get_propagation() == "back":
            MatcherPlotterPhase._plot_back(axes, sbs_file, plane)

        axes.legend(loc="lower left", prop={'size': 16})
        axes.set_ylabel(r"$\Delta\Phi_{" + plane.lower() + "}$")

    @staticmethod
    def _plot_front(axes, sbs_file, plane):
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, "PROPPHASE" + plane).values,
                      getattr(sbs_file, "ERRPROPPHASE" + plane).values,
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "PROPPHASE" + plane).values,
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, "CORPHASE" + plane).values,
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "CORPHASE" + plane).values,
                  marker="o", markersize=7., color="green")

    @staticmethod
    def _plot_back(axes, sbs_file, plane):
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, "BACKPHASE" + plane).values,
                      getattr(sbs_file, "ERRBACKPHASE" + plane).values,
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "BACKPHASE" + plane).values,
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, "BACKCORPHASE" + plane).values,
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "BACKCORPHASE" + plane).values,
                  marker="o", markersize=7., color="green")
