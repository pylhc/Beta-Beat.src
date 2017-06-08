import os
from matcher_model_default import MatcherModelDefault
from matcher_model_default import MatcherPlotterDefault
from matchers.amp_matcher import AmpMatcher
from Utilities import tfs_pandas


class MatcherModelAmp(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = super(MatcherModelAmp, self).get_matcher_dict()
        matcher_dict["type"] = "amp"
        return matcher_dict

    def create_matcher(self, match_path):
        self._matcher = AmpMatcher(self._name, self.get_matcher_dict(), match_path)

    def get_plotter(self, figure):
        if self._plotter is None:
            self._plotter = MatcherPlotterAmp(figure, self)
        return self._plotter


class MatcherPlotterAmp(MatcherPlotterDefault):

    def plot(self):
        self._figure.clear()
        self._axes_to_data = {}
        num_beams = len(self._matcher_model._matcher.get_beams())
        if 1 in self._matcher_model._matcher.get_beams():
            axes_b1_x = self._figure.add_subplot(2, num_beams, 1)
            axes_b1_x.set_title("Beam 1")
            amp_bb_beam1_horizontal = self._get_amp_beta_beat_file(1, "X")
            bb_beam1_horizontal = self._get_beta_beat_file(1, "X")
            self._axes_to_data[axes_b1_x] = bb_beam1_horizontal
            self._plot_match(axes_b1_x, amp_bb_beam1_horizontal, bb_beam1_horizontal, "X")

            axes_b1_y = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 3)
            amp_bb_beam1_vertical = self._get_amp_beta_beat_file(1, "Y")
            bb_beam1_vertical = self._get_beta_beat_file(1, "Y")
            self._axes_to_data[axes_b1_y] = bb_beam1_vertical
            self._plot_match(axes_b1_y, amp_bb_beam1_vertical, bb_beam1_vertical, "Y")

        if 2 in self._matcher_model._matcher.get_beams():
            axes_b2_x = self._figure.add_subplot(2, num_beams, 1 if num_beams == 1 else 2)
            axes_b2_x.set_title("Beam 2")
            amp_bb_beam2_horizontal = self._get_amp_beta_beat_file(2, "X")
            bb_beam2_horizontal = self._get_beta_beat_file(2, "X")
            self._axes_to_data[axes_b2_x] = bb_beam2_horizontal
            self._plot_match(axes_b2_x, amp_bb_beam2_horizontal, bb_beam2_horizontal, "X")

            axes_b2_y = self._figure.add_subplot(2, num_beams, 2 if num_beams == 1 else 4)
            amp_bb_beam2_vertical = self._get_amp_beta_beat_file(2, "Y")
            bb_beam2_vertical = self._get_beta_beat_file(2, "Y")
            self._axes_to_data[axes_b2_y] = bb_beam2_vertical
            self._plot_match(axes_b2_y, amp_bb_beam2_vertical, bb_beam2_vertical, "Y")
        self._figure.tight_layout()
        self._figure.patch.set_visible(False)
        self._figure.canvas.draw()

    def _get_amp_beta_beat_file(self, beam, plane):
        amp_bb_path = getattr(self._matcher_model, "get_beam" + str(beam) + "_output_path")()
        file_data = tfs_pandas.read_tfs(os.path.join(
            amp_bb_path, "sbs",
            "sbsampbetabeat" + plane.lower() + "_IP" + str(self._matcher_model.get_ip()) + ".out")
        )
        return file_data

    def _get_beta_beat_file(self, beam, plane):
        bb_path = getattr(self._matcher_model, "get_beam" + str(beam) + "_output_path")()
        file_data = tfs_pandas.read_tfs(os.path.join(
            bb_path, "sbs",
            "sbsbetabeat" + plane.lower() + "_IP" + str(self._matcher_model.get_ip()) + ".out")
        )
        return file_data

    def _plot_match(self, axes, amp_bb_data, bb_data, plane):
        if self._matcher_model.get_propagation() == "front":
            MatcherPlotterAmp._plot_front(axes, amp_bb_data, bb_data, plane)
        elif self._matcher_model.get_propagation() == "back":
            MatcherPlotterAmp._plot_back(axes, amp_bb_data, bb_data, plane)

        axes.legend(loc="upper left", prop={'size': 16})
        axes.set_ylabel(r"$\Delta\beta_{" + plane.lower() + r"} / {\beta_{model}}$")
        axes.set_xlabel("S along the segment [m]")

    @staticmethod
    def _plot_front(axes, amp_bb_data, bb_data, plane):
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCOR" + plane).values,
                  label=r"$\Delta\beta / {\beta_{model}}$ model",
                  color="green")
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCOR" + plane).values,
                  marker="o", markersize=7., color="green")

        axes.errorbar(amp_bb_data.S,
                      getattr(amp_bb_data, "BETABEATAMP" + plane).values,
                      getattr(amp_bb_data, "ERRBETABEATAMP" + plane).values,
                      fmt='o', markersize=8.,
                      label=r"$\Delta\beta / \beta_{model}$ amplitude",
                      color="blue")

    @staticmethod
    def _plot_back(axes, amp_bb_data, bb_data, plane):
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCORBACK" + plane).values,
                  label=r"$\Delta\beta / {\beta_{model}}$ model", color="green")
        axes.plot(bb_data.S, getattr(bb_data, "BETABEATCORBACK" + plane).values,
                  marker="o", markersize=7., color="green")

        axes.errorbar(amp_bb_data.S,
                      getattr(amp_bb_data, "BETABEATAMPBACK" + plane).values,
                      getattr(amp_bb_data, "ERRBETABEATAMPBACK" + plane).values,
                      fmt='o', markersize=8.,
                      label=r"$\Delta\beta / \{beta_{model}}$ amplitude",
                      color="blue")
