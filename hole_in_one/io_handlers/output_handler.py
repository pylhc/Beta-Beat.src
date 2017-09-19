import os
import logging
from sdds_files import turn_by_turn_reader


LOGGER = logging.getLogger("__name__")

RAW_SUFFIX = ".raw"
CLEAN_SUFFIX = ".clean"

HEADERS = {"X": ["NAME", "S", "BINDEX", "SLABEL", "TUNEX", #"TUNEZ"
                 "NOISE", "PK2PK", "CO", "CORMS", "AMPX",
                 "MUX", "AVG_AMPX", "AVG_MUX", "BPM_RES"],
           "Y": ["NAME", "S", "BINDEX", "SLABEL", "TUNEY", #"TUNEZ"
                 "NOISE", "PK2PK", "CO", "CORMS",
                 "AMPY", "MUY", "AVG_AMPY", "AVG_MUY", "BPM_RES"]}

SPECTR_COLUMN_NAMES = ["FREQ", "AMP"]


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
    bad_bpms_file = get_outpath_with_suffix(bin_path, output_dir, ".bad_bpms_" + plane)
    with open(bad_bpms_file, 'w') as bad_bpms_writer:
        for line in bad_bpms_with_reasons:
            bad_bpms_writer.write(line + '\n')


def write_harpy_output(harpy_data_frame):
    pass


def _create_lin_files(self):
    file_name = self._get_outfile_name(self._plane)
    lin_outfile = tfs_file_writer.TfsFileWriter(
        os.path.join(self._output_dir, file_name)
    )
    headers = HEADERS[self._plane]
    for resonance in RESONANCE_LISTS[self._plane]:
        if resonance == MAIN_LINES[self._plane]:
            continue
        x, y, z = resonance
        if z == 0:
            resstr = (str(x) + str(y)).replace("-", "_")
        else:
            resstr = (str(x) + str(y) + str(z)).replace("-", "_")
        headers.extend(["AMP" + resstr, "PHASE" + resstr])
    headers.extend(["NATTUNE" + self._plane,
                    "NATAMP" + self._plane])
    lin_outfile.add_column_names(headers)
    lin_outfile.add_column_datatypes(
        ["%s"] + ["%le"] * (len(headers) - 1))
    self._lin_outfile = lin_outfile


def _write_single_bpm_results(self, lin_outfile, bpm_results):
    row = [bpm_results.name, bpm_results.position, 0, 0, bpm_results.tune,
           0, bpm_results.peak_to_peak, bpm_results.closed_orbit,
           bpm_results.closed_orbit_rms, bpm_results.amplitude,
           bpm_results.phase, bpm_results.amp_from_avg,
           bpm_results.phase_from_avg, bpm_results.bpm_resolution]
    resonance_list = RESONANCE_LISTS[self._plane]
    main_resonance = MAIN_LINES[self._plane]
    for resonance in resonance_list:
        if resonance != main_resonance:
            if resonance in bpm_results.resonances:
                _, coefficient = bpm_results.resonances[resonance]
                row.append(np.abs(coefficient) / bpm_results.amplitude)
                row.append(np.angle(coefficient) / (2 * np.pi))
            else:
                row.append(0.0)
                row.append(0.0)

    col_name = "NAT" + self._plane.upper()
    try:
        natural_freq, natural_coef = bpm_results.resonances[col_name]
        row.append(natural_freq)
        row.append(np.abs(natural_coef) / bpm_results.amplitude)
    except KeyError:
        row.append(0.0)
        row.append(0.0)
    lin_outfile.add_table_row(row)


def get_outpath_with_suffix(path, output_dir, suffix):
    return os.path.join(output_dir, os.path.basename(path) + suffix)
