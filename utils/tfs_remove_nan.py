from utils import tfs_pandas as tfs
from utils import logging_tools as logtools
import os

LOG = logtools.get_logger(__name__)


def clean_files(list_of_files, replace=False):
    for filepath in list_of_files:
        try:
            df = tfs.read_tfs(filepath)
            LOG.info("Read file {:s}".format(filepath))
        except (IOError, tfs.TfsFormatError):
            LOG.info("Skipped file {:s}".format(filepath))
        else:
            df = df.dropna(axis='index')
            if not replace:
                filepath += ".dropna"
            tfs.write_tfs(filepath, df)


if __name__ == "__main__":
    dirs = ["/media/jdilly/Storage/Repositories/Gui_Output/2017-12-06/LHCB1/Results/15-45-45_WANALYSIS_SUSSIX_2"]

    # dirs = [os.path.join("/", "afs", "cern.ch", "user", "j", "jdilly", dir)
    #         for dir in ["w_afterGlobal", "w_afterGlobalCorrection"]]


    for dir_path in dirs:
        all_files = os.listdir(dir_path)
        clean_files([os.path.join(dir_path, f) for f in all_files], replace=True)
