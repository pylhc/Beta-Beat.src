import os
from matcher_model_default import MatcherModelDefault
from matcher_model_default import MatcherPlotterDefault
from matchers.coupling_matcher import CouplingMatcher
from Utilities import tfs_pandas


class MatcherModelCoupling(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = super(MatcherModelCoupling, self).get_matcher_dict()
        matcher_dict["type"] = "coupling"
        return matcher_dict

    def create_matcher(self, match_path):
        self._matcher = CouplingMatcher(
            self._name,
            self.get_matcher_dict(),
            match_path
        )

    def get_plotter(self, figure):
        if self._plotter is None:
            self._plotter = MatcherPlotterCoupling(figure, self)
        return self._plotter


class MatcherPlotterCoupling(MatcherPlotterDefault):

    def plot(self):
        self._figure.clear()
        self._axes_to_data = {}
        num_beams = len(self._matcher_model._matcher.get_beams())
        if 1 in self._matcher_model._matcher.get_beams():
            axes_b1_f1001 = self._figure.add_subplot(2, num_beams, 1)
            axes_b1_f1001.set_title("Beam 1")
            axes_b1_f1010 = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 3)
            file_beam1 = tfs_pandas.read_tfs(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbscouple_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._axes_to_data[axes_b1_f1001] = file_beam1
            self._axes_to_data[axes_b1_f1010] = file_beam1
            self._plot_match(axes_b1_f1001, file_beam1, "f1001")
            self._plot_match(axes_b1_f1010, file_beam1, "f1010")
        if 2 in self._matcher_model._matcher.get_beams():
            axes_b2_f1001 = self._figure.add_subplot(2, num_beams, 1 if num_beams == 1 else 2)
            axes_b2_f1001.set_title("Beam 2")
            axes_b2_f1010 = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 4)
            file_beam2 = tfs_pandas.read_tfs(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbscouple_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._axes_to_data[axes_b2_f1001] = file_beam2
            self._axes_to_data[axes_b2_f1010] = file_beam2
            self._plot_match(axes_b2_f1001, file_beam2, "f1001")
            self._plot_match(axes_b2_f1010, file_beam2, "f1010")
        self._figure.tight_layout()
        self._figure.patch.set_visible(False)
        self._figure.canvas.draw()

    def _plot_match(self, axes, sbs_file, term):
        if self._matcher_model.get_propagation() == "front":
            MatcherPlotterCoupling._plot_front(axes, sbs_file, term)
        elif self._matcher_model.get_propagation() == "back":
            MatcherPlotterCoupling._plot_back(axes, sbs_file, term)

        axes.legend(loc="lower left", prop={'size': 16})
        if term == "f1001":
            axes.set_ylabel(r"$abs(f_{1001})$")
        elif term == "f1010":
            axes.set_ylabel(r"$abs(f_{1010})$")

    @staticmethod
    def _plot_front(axes, sbs_file, term):
        up_term = term.upper()
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, up_term + "ABSMEAS").values,
                      getattr(sbs_file, "ERR" + up_term + "ABSMEAS").values,
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSMEAS").values,
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSCOR").values,
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSCOR").values,
                  marker="o", markersize=7., color="green")

    @staticmethod
    def _plot_back(axes, sbs_file, term):
        up_term = term.upper()
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, up_term + "ABSBACK").values,
                      getattr(sbs_file, "ERR" + up_term + "ABSBACK").values,
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSBACK").values,
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSBACKCOR").values,
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSBACKCOR").values,
                  marker="o", markersize=7., color="green")
