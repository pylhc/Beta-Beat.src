import os
import logging
import cPickle
from Utilities import tfs_pandas
from sdds_files import turn_by_turn_reader

LOGGER = logging.getLogger("__name__")

RAW_SUFFIX = ".raw"
CLEAN_SUFFIX = ".clean"


def write_raw_file(tbt_file, main_input):
    LOGGER.debug("Writing raw sdds")
    tbt_file.write_to_ascii(
        main_input.model,
        get_outpath_with_suffix(main_input.file,
                                main_input.outputdir,
                                RAW_SUFFIX)
    )


class CleanedAsciiWritter(object):

    def __init__(self, main_input, date):
        self._main_input = main_input
        self._date = date

        self.bpm_names_x = None
        self.samples_matrix_x = None
        self.bpm_names_y = None
        self.samples_matrix_y = None
        self.dpp = None

    def write(self):
        LOGGER.debug("Writing clean sdds")
        headers_dict = {}
        if self.dpp is not None:
            headers_dict["dpp"] = self.dpp
        turn_by_turn_reader.write_ascii_file(
            self._main_input.model,
            get_outpath_with_suffix(self._main_input.file,
                                    self._main_input.outputdir,
                                    CLEAN_SUFFIX),
            list(self.bpm_names_x), self.samples_matrix_x,
            list(self.bpm_names_y), self.samples_matrix_y,
            self._date, headers_dict,
        )


def write_bad_bpms(bin_path, bad_bpms_with_reasons, output_dir, plane):
    bad_bpms_file = get_outpath_with_suffix(
        bin_path,
        output_dir,
        ".bad_bpms_" + plane
    )
    with open(bad_bpms_file, 'w') as bad_bpms_writer:
        for line in bad_bpms_with_reasons:
            bad_bpms_writer.write(line + '\n')


def write_harpy_output(main_input, harpy_data_frame, headers, spectrum, plane):
    output_file = get_outpath_with_suffix(
        main_input.file, main_input.outputdir, ".lin" + plane
    )
    tfs_pandas.write_tfs(harpy_data_frame, headers, output_file)
    _dump(get_outpath_with_suffix(
          main_input.file, main_input.outputdir, ".spec" + plane), spectrum)


def _dump(output_file, content):
    with open(output_file, 'wb') as output_data:
        cPickle.Pickler(output_data, -1).dump(content)


def get_outpath_with_suffix(path, output_dir, suffix):
    return os.path.join(output_dir, os.path.basename(path) + suffix)
