import sys
import os
from matcher_model_default import MatcherModelDefault
from matchers.coupling_matcher import CouplingMatcher
from Python_Classes4MAD import metaclass


class MatcherModelCoupling(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = super(MatcherModelCoupling, self).get_matcher_dict()
        matcher_dict["type"] = "coupling"
        return matcher_dict

    def create_matcher(self, match_path):
        self._matcher = CouplingMatcher(self._name, self.get_matcher_dict(), match_path)

    def get_plotter(self, figures):
        return MatcherPlotterCoupling(figures, self)


class MatcherPlotterCoupling(object):

    def __init__(self, figures, matcher_model):
        assert type(matcher_model) is MatcherModelCoupling
        self._figures = figures
        self._matcher_model = matcher_model

    def plot(self):
        if 1 in self._matcher_model._matcher.get_beams():
            file_beam1 = metaclass.twiss(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbscouple_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            figure_b1_x = self._figures[0][0]
            figure_b1_x.clear()
            self._plot_match(figure_b1_x, file_beam1, "f1001")
            figure_b1_y = self._figures[0][1]
            figure_b1_y.clear()
            self._plot_match(figure_b1_y, file_beam1, "f1010")
        if 2 in self._matcher_model._matcher.get_beams():
            file_beam2 = metaclass.twiss(os.path.join(
                self._matcher_model.get_beam2_output_path(), "sbs",
                "sbscouple_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            figure_b2_x = self._figures[1][0]
            figure_b2_x.clear()
            self._plot_match(figure_b2_x, file_beam2, "f1001")
            figure_b2_y = self._figures[1][1]
            figure_b2_y.clear()
            self._plot_match(figure_b2_y, file_beam2, "f1010")

    def _plot_match(self, figure, sbs_file, term):
        ax = figure.add_subplot(1, 1, 1)

        if self._matcher_model.get_propagation() == "front":
            MatcherPlotterCoupling._plot_front(ax, sbs_file, term)
        elif self._matcher_model.get_propagation() == "back":
            MatcherPlotterCoupling._plot_back(ax, sbs_file, term)

        ax.legend(loc="lower left", prop={'size': 16})
        if term == "f1001":
            ax.set_ylabel(r"$abs(f_{1001})$")
        elif term == "f1010":
            ax.set_ylabel(r"$abs(f_{1010})$")
        figure.patch.set_visible(False)
        figure.canvas.draw()

    @staticmethod
    def _plot_front(axes, sbs_file, term):
        up_term = term.upper()
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, up_term + "ABSMEAS"),
                      getattr(sbs_file, "ERR" + up_term + "ABSMEAS"),
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSMEAS"),
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSCOR"),
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSCOR"),
                  marker="o", markersize=7., color="green")

    @staticmethod
    def _plot_back(axes, sbs_file, term):
        up_term = term.upper()
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, up_term + "ABSBACK"),
                      getattr(sbs_file, "ERR" + up_term + "ABSBACK"),
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSBACK"),
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSBACKCOR"),
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, up_term + "ABSBACKCOR"),
                  marker="o", markersize=7., color="green")


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
