import os

from tfs_files.tfs_collection import TfsCollection, Tfs


class OpticsMeasurement(TfsCollection):
    """Class to hold and load the measurements from GetLLM.

    The class will try to load the _free file, then the _free2 and then the
    normal file, if none of them if present an IOError will be raised.

    Arguments:
        directory: The path to the measurement directory, usually a GetLLM
            output directory.
    """
    beta = Tfs("getbeta")
    amp_beta = Tfs("getampbeta")
    kmod_beta = Tfs("getkmodbeta")
    phase = Tfs("getphase")
    phasetot = Tfs("getphasetot")
    disp = Tfs("getD")
    coupling = Tfs("getcouple", two_planes=False)
    norm_disp = Tfs("getNDx", two_planes=False)

    def get_filename(self, prefix, plane=""):
        templ = prefix + "{}{}.out"
        for filename in (templ.format(plane, "_free"),
                         templ.format(plane, "_free2"),
                         templ.format(plane, "")):
            if os.path.isfile(os.path.join(self.directory, filename)):
                return filename
        raise IOError("No file name found for prefix {} in {}."
                      .format(prefix, self.directory))

    def write_to(self, value, prefix, plane=""):
        data_frame, suffix = value
        templ = prefix + "{}{}.out"
        filename = templ.format(plane, suffix)
        return filename, data_frame