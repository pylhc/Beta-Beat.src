from __future__ import print_function
import sys
import os


CURRENT_DIR = os.path.dirname(__file__)


class LhcExcitationMode(object):
    FREE, ACD, ADT = range(3)


class Lhc(object):
    INT_TUNE_X = 64.
    INT_TUNE_Y = 59.

    def __init__(self):
        self._beam = None
        self.optics_file = None
        self.nat_tune_x = None
        self.nat_tune_y = None
        self._excitation = None
        self.drv_tune_x = None
        self.drv_tune_y = None
        self.energy = None
        self.dpp = 0.0

    def verify_object(self):  # TODO: Maybe more checks?
        if self.excitation is None:
            raise LhcError("Excitation mode not set.")
        if (self.excitation == LhcExcitationMode.ACD or
                self.excitation == LhcExcitationMode.ADT):
            if self.drv_tune_x is None or self.drv_tune_y is None:
                raise LhcError("Driven tunes not set.")

    def get_nominal_tmpl(self):
        return self._get_file("nominal.madx")

    def get_best_knowledge_tmpl(self):
        return self._get_file("best_knowledge.madx")

    def _get_file(self, filename):
        return os.path.join(CURRENT_DIR, "lhc", filename)

    @property
    def beam(self):
        return self._beam

    @beam.setter
    def beam(self, beam_num):
        if beam_num not in (None, 1, 2):
            raise LhcError("Beam should be 1 or 2")
        self._beam = beam_num

    @property
    def excitation(self):
        return self._excitation

    @excitation.setter
    def excitation(self, excitation_mode):
        if excitation_mode not in (LhcExcitationMode.FREE,
                                   LhcExcitationMode.ACD,
                                   LhcExcitationMode.ADT):
            raise ValueError("Wrong excitation mode.")
        self._excitation = excitation_mode


class LhcAts(Lhc):
    INT_TUNE_X = 62.
    INT_TUNE_Y = 60.


class LhcRunI(Lhc):
    MACROS_NAME = "lhc_runI"


class LhcRunII2015(Lhc):
    MACROS_NAME = "lhc_runII"


class LhcRunII2016(Lhc):
    MACROS_NAME = "lhc_runII_2016"


class LhcRunII2016Ats(LhcAts):
    MACROS_NAME = "lhc_runII_2016_ats"


class HlLhc(LhcAts):
    MACROS_NAME = "hllhc"


class LhcError(Exception):
    """
    Raised when an LHC creation error occurs.
    """
    pass


if __name__ == "__main__":
    print("Import this module.", file=sys.stderr)
