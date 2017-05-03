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

    def create_matcher(self, match_path):
        self._matcher = PhaseMatcher(
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
        num_beams = len(self._matcher_model._matcher.get_beams())
        if 1 in self._matcher_model._matcher.get_beams():
            axes_b1_x = self._figure.add_subplot(2, num_beams, 1)
            axes_b1_x.set_title("Beam 1")
            file_beam1_horizontal = tfs_pandas.read_tfs(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbsphasext_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._axes_to_data[axes_b1_x] = file_beam1_horizontal
            self._plot_match(axes_b1_x, file_beam1_horizontal, "X")
            file_beam1_vertical = tfs_pandas.read_tfs(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbsphaseyt_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            axes_b1_y = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 3)
            self._axes_to_data[axes_b1_y] = file_beam1_vertical
            self._plot_match(axes_b1_y, file_beam1_vertical, "Y")
        if 2 in self._matcher_model._matcher.get_beams():
            axes_b2_x = self._figure.add_subplot(2, num_beams, 1 if num_beams == 1 else 2)
            axes_b2_x.set_title("Beam 2")
            file_beam2_horizontal = tfs_pandas.read_tfs(os.path.join(
                self._matcher_model.get_beam2_output_path(), "sbs",
                "sbsphasext_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._axes_to_data[axes_b2_x] = file_beam2_horizontal
            self._plot_match(axes_b2_x, file_beam2_horizontal, "X")
            axes_b2_y = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 4)
            file_beam2_vertical = tfs_pandas.read_tfs(os.path.join(
                self._matcher_model.get_beam2_output_path(), "sbs",
                "sbsphaseyt_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._axes_to_data[axes_b2_y] = file_beam2_vertical
            self._plot_match(axes_b2_y, file_beam2_vertical, "Y")
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
