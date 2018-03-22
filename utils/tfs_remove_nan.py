from utils import tfs_pandas as tfs
from utils import logging_tools as logtools

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


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))

