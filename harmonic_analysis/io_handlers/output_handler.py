from __future__ import print_function
import os
import logging
from tfs_files import tfs_pandas
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

        self.samples_matrix_x = None
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
            self.samples_matrix_x.index, self.samples_matrix_x,
            self.samples_matrix_y.index, self.samples_matrix_y,
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
    tfs_pandas.write_tfs(output_file, harpy_data_frame, headers)
    if not main_input.skip_files:
        _write_full_spectrum(main_input, spectrum, plane)


def _write_full_spectrum(main_input, spectrum, plane):
    spectr_amps_files = get_outpath_with_suffix(
        main_input.file, main_input.outputdir, ".amps" + plane
    )
    amps_df = spectrum["COEFS"].abs().T
    tfs_pandas.write_tfs(spectr_amps_files, amps_df)
    spectr_freqs_files = get_outpath_with_suffix(
        main_input.file, main_input.outputdir, ".freqs" + plane
    )
    freqs_df = spectrum["FREQS"].T
    tfs_pandas.write_tfs(spectr_freqs_files, freqs_df)


def get_outpath_with_suffix(path, output_dir, suffix):
    return os.path.join(output_dir, os.path.basename(path) + suffix)
