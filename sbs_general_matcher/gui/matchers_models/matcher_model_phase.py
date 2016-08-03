import sys
import os
from matcher_model_default import MatcherModelDefault
from matchers.phase_matcher import PhaseMatcher
from Python_Classes4MAD import metaclass


class MatcherModelPhase(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = super(MatcherModelPhase, self).get_matcher_dict()
        matcher_dict["type"] = "phase"
        return matcher_dict

    def create_matcher(self, match_path):
        self._matcher = PhaseMatcher(self._name, self.get_matcher_dict(), match_path)

    def get_plotter(self, figures):
        return MatcherPlotterPhase(figures, self)


class MatcherPlotterPhase(object):

    def __init__(self, figures, matcher_model):
        assert type(matcher_model) is MatcherModelPhase
        self._figures = figures
        self._matcher_model = matcher_model

    def plot(self):
        if 1 in self._matcher_model._matcher.get_beams():
            figure_b1_x = self._figures[0][0]
            figure_b1_x.clear()
            file_beam1_horizontal = metaclass.twiss(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbsphasext_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._plot_match(figure_b1_x, file_beam1_horizontal, "X")
            figure_b1_y = self._figures[0][1]
            figure_b1_y.clear()
            file_beam1_vertical = metaclass.twiss(os.path.join(
                self._matcher_model.get_beam1_output_path(), "sbs",
                "sbsphaseyt_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._plot_match(figure_b1_y, file_beam1_vertical, "Y")
        if 2 in self._matcher_model._matcher.get_beams():
            figure_b2_x = self._figures[1][0]
            figure_b2_x.clear()
            file_beam2_horizontal = metaclass.twiss(os.path.join(
                self._matcher_model.get_beam2_output_path(), "sbs",
                "sbsphasext_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._plot_match(figure_b2_x, file_beam2_horizontal, "X")
            figure_b2_y = self._figures[1][1]
            figure_b2_y.clear()
            file_beam2_vertical = metaclass.twiss(os.path.join(
                self._matcher_model.get_beam2_output_path(), "sbs",
                "sbsphaseyt_IP" + str(self._matcher_model.get_ip()) + ".out")
            )
            self._plot_match(figure_b2_y, file_beam2_vertical, "Y")

    def _plot_match(self, figure, sbs_file, plane):
        ax = figure.add_subplot(1, 1, 1)

        if self._matcher_model.get_propagation() == "front":
            MatcherPlotterPhase._plot_front(ax, sbs_file, plane)
        elif self._matcher_model.get_propagation() == "back":
            MatcherPlotterPhase._plot_back(ax, sbs_file, plane)

        ax.legend(loc="lower left", prop={'size': 16})
        ax.set_ylabel(r"$\Delta\Phi$")
        figure.patch.set_visible(False)
        figure.canvas.draw()

    @staticmethod
    def _plot_front(axes, sbs_file, plane):
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, "PROPPHASE" + plane),
                      getattr(sbs_file, "ERRPROPPHASE" + plane),
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "PROPPHASE" + plane),
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, "CORPHASE" + plane),
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "CORPHASE" + plane),
                  marker="o", markersize=7., color="green")

    @staticmethod
    def _plot_back(axes, sbs_file, plane):
        axes.errorbar(sbs_file.S,
                      getattr(sbs_file, "BACKPHASE" + plane),
                      getattr(sbs_file, "ERRBACKPHASE" + plane),
                      label=r"$\Delta\Phi$ measured", color="blue")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "BACKPHASE" + plane),
                  marker="o", markersize=7., color="blue")

        axes.plot(sbs_file.S,
                  getattr(sbs_file, "BACKCORPHASE" + plane),
                  label=r"$\Delta\Phi$ model", color="green")
        axes.plot(sbs_file.S,
                  getattr(sbs_file, "BACKCORPHASE" + plane),
                  marker="o", markersize=7., color="green")


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
