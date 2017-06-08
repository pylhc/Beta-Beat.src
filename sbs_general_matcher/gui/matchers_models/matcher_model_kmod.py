import sys
import os
from matcher_model_default import MatcherModelDefault
from matcher_model_default import MatcherPlotterDefault
from matchers.kmod_matcher import KmodMatcher
from Utilities import tfs_pandas


class MatcherModelKmod(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = super(MatcherModelKmod, self).get_matcher_dict()
        matcher_dict["type"] = "kmod"
        return matcher_dict

    def create_matcher(self, match_path):
        self._matcher = KmodMatcher(self._name, self.get_matcher_dict(), match_path)

    def get_plotter(self, figure):
        if self._plotter is None:
            self._plotter = MatcherPlotterKmod(figure, self)
        return self._plotter


class MatcherPlotterKmod(MatcherPlotterDefault):

    def plot(self):
        self._figure.clear()
        self._axes_to_data = {}
        num_beams = len(self._matcher_model._matcher.get_beams())
        if 1 in self._matcher_model._matcher.get_beams():
            axes_b1_x = self._figure.add_subplot(2, num_beams, 1)
            axes_b1_x.set_title("Beam 1")
            kmod_bb_beam1_horizontal = self._get_kmod_beta_beat_file(1, "X")
            bb_beam1_horizontal = self._get_beta_beat_file(1, "X")
            self._axes_to_data[axes_b1_x] = bb_beam1_horizontal
            self._plot_match(axes_b1_x, kmod_bb_beam1_horizontal, bb_beam1_horizontal, "X")

            axes_b1_y = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 3)
            kmod_bb_beam1_vertical = self._get_kmod_beta_beat_file(1, "Y")
            bb_beam1_vertical = self._get_beta_beat_file(1, "Y")
            self._axes_to_data[axes_b1_y] = bb_beam1_vertical
            self._plot_match(axes_b1_y, kmod_bb_beam1_vertical, bb_beam1_vertical, "Y")

        if 2 in self._matcher_model._matcher.get_beams():
            axes_b2_x = self._figure.add_subplot(2, num_beams, 1 if num_beams == 1 else 2)
            axes_b2_x.set_title("Beam 2")
            kmod_bb_beam2_horizontal = self._get_kmod_beta_beat_file(2, "X")
            bb_beam2_horizontal = self._get_beta_beat_file(2, "X")
            self._axes_to_data[axes_b2_x] = bb_beam2_horizontal
            self._plot_match(axes_b2_x, kmod_bb_beam2_horizontal, bb_beam2_horizontal, "X")

            axes_b2_y = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 4)
            kmod_bb_beam2_vertical = self._get_kmod_beta_beat_file(2, "Y")
            bb_beam2_vertical = self._get_beta_beat_file(2, "Y")
            self._axes_to_data[axes_b2_y] = bb_beam2_vertical
            self._plot_match(axes_b2_y, kmod_bb_beam2_vertical, bb_beam2_vertical, "Y")
        self._figure.tight_layout()
        self._figure.patch.set_visible(False)
        self._figure.canvas.draw()

    def _get_kmod_beta_beat_file(self, beam, plane):
        kmod_bb_path = getattr(self._matcher_model, "get_beam" + str(beam) + "_output_path")()
        file_data = tfs_pandas.read_tfs(os.path.join(
            kmod_bb_path, "sbs",
            "sbskmodbetabeat" + plane.lower() + "_IP" + str(self._matcher_model.get_ip()) + ".out")
        )
        return file_data

    def _get_beta_beat_file(self, beam, plane):
        bb_path = getattr(self._matcher_model, "get_beam" + str(beam) + "_output_path")()
        file_data = tfs_pandas.read_tfs(os.path.join(
            bb_path, "sbs",
            "sbsbetabeat" + plane.lower() + "_IP" + str(self._matcher_model.get_ip()) + ".out")
        )
        return file_data

    def _plot_match(self, axes, kmod_bb_data, bb_data, plane):
        if self._matcher_model.get_propagation() == "front":
            MatcherPlotterKmod._plot_front(axes, kmod_bb_data, bb_data, plane)
        elif self._matcher_model.get_propagation() == "back":
            MatcherPlotterKmod._plot_back(axes, kmod_bb_data, bb_data, plane)

        axes.legend(loc="upper left", prop={'size': 16})
        axes.set_ylabel(r"$\Delta\beta_{" + plane.lower() + r"} / {\beta_{model}}$")
        axes.set_xlabel("S along the segment [m]")

    @staticmethod
    def _plot_front(axes, kmod_bb_data, bb_data, plane):
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCOR" + plane).values,
                  label=r"$\Delta\beta / {\beta_{model}}$ model", color="green")
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCOR" + plane).values,
                  marker="o", markersize=7., color="green")

        axes.errorbar(kmod_bb_data.S,
                      getattr(kmod_bb_data, "BETABEAT" + plane).values,
                      getattr(kmod_bb_data, "ERRBETABEAT" + plane).values,
                      fmt='o', markersize=8.,
                      label=r"$\Delta\beta / {\beta_{model}}$ k-mod",
                      color="blue")

    @staticmethod
    def _plot_back(axes, kmod_bb_data, bb_data, plane):
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCORBACK" + plane).values,
                  label=r"$\Delta\beta / {\beta_{model}}$ model", color="green")
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCORBACK" + plane).values,
                  marker="o", markersize=7., color="green")

        axes.errorbar(kmod_bb_data.S,
                      getattr(kmod_bb_data, "BETABEATBACK" + plane).values,
                      getattr(kmod_bb_data, "ERRBETABEATBACK" + plane).values,
                      fmt='o', markersize=8.,
                      label=r"$\Delta\beta / {\beta_{model}}$ k-mod",
                      color="blue")


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
