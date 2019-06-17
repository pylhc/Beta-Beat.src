"""
Module online_model.extractor_wrapper
---------------------------------------
Python wrapper to use the online model extractor with some conversion functionality.

.. warning::
    Needs to be run from the technical network or using a cwd with server and local access.
"""
import os
import subprocess

import online_model.constants as const
import online_model.data_converter as dc
from tfs_files import tfs_pandas as tfs
from utils import logging_tools
from plotshop.plot_tfs import plot


LOG = logging_tools.get_logger(__name__)


def extract_knob_value_and_definition(knob_names, time, cwd="./", server=None):
    """ Extract knob values from online model using the knob name and the time.

    Example Call::

        extract_knob_values("LHCBEAM/2018_global_ats_flat_b1_for_ip5_waist",
                            "2018-10-30 15:00:00.0",
                            cwd="/afs/cern.ch/user/j/jdilly/extractor/",
                            server="cs-ccr-dev3")

    Args:
        knob_names (list): List of knob names to extract
        time: UTC time in ISOformat to extract
        cwd: output directory for results and log (default: current directory)
        server: server to run on (default: runs local)
    """
    knob_file = os.path.join(cwd, const.get_extractor_knobs_filename())

    # first get the trim value and madx_changefile
    _run_knob_extraction(knob_names, time, cwd, server, "trim")
    if not os.path.exists(knob_file):
        raise IOError("Something went wrong while extracting data form online model. See log file.")

    trim, _ = dc.post_trim_extract(knob_file)

    # run again and get the values and deltas
    _run_knob_extraction(knob_names, time, cwd, server, "k")
    df_list = dc.knobs_k_to_tfs(knob_file, trim)

    for knob, df in zip(knob_names, df_list):
        filename = const.get_knob_tfs_filename(knob)
        tfs.write_tfs(os.path.join(cwd, filename), df)


def extract_overview(knob_names, time=None, cwd="./", server=None, show_plot=False):
    """ Extract overview-data consisting of
        - beamprocess
        - optics used
        - values of given Knobs
        - orbit plot

        Args:
            knob_names (list): List of knob names to extract
            time: UTC time in ISO format (default: now)
            cwd: output directory for results and log (default: current directory)
            server: server to run on (default: runs local)
    """

    if time is None:
        time = dc.get_utc_now()

    if knob_names is None or len(knob_names) == 0:
        raise NotImplementedError("Knob names need to be provided, due to bug in extractor.")
        # knob_names = extract_all_knob_names(time=time, only_active=False, cwd=cwd, server=server)

    # extraction
    _run_overview_extraction(knob_names, time, cwd, server)

    # knobs
    knobs_file = os.path.join(cwd, const.get_extractor_knobs_filename())
    if not os.path.exists(knobs_file):
        raise RuntimeError("Knobs output file '{:s}' ".format(knobs_file) +
                           "not found!")
    df = dc.knobs_kvalues_to_tfs(knobs_file)

    # optics
    optics_file = os.path.join(cwd, const.get_extractor_output_filename())
    if not os.path.exists(optics_file):
        raise RuntimeError("Optics output file '{:s}' ".format(optics_file) +
                           "not found!")
    df.headers[const.get_optics_header()] = dc.get_optics(optics_file)

    # orbit
    filenames = [os.path.join(cwd, const.get_default_orbit_filename(b)) for b in [1, 2]]
    if not all([os.path.exists(fn) for fn in filenames]):
        raise RuntimeError("Orbit files not found in '{:s}'.".format(cwd))
    df.headers[const.get_fill_header()] = dc.get_fill_from_orbitfile(filenames[0])

    tfs.write_tfs(os.path.join(cwd, const.get_overview_filename()), df)
    _log_df(df)
    _plot_orbits(cwd, show_plot)


def extract_all_knob_names(time=None, only_active=True, cwd="./", server=None):
    """
    DOES NOT WORK BECAUSE OF ONLINE-MODEL-EXTRACTOR BUGS....
    Extract all knobs present at given time.

     Args:
        time: UTC time in ISO format (default: now)
        only_active (bool): If set, returns only knobs with non-zero value.
        cwd: output directory for results and log (default: current directory)
        server: server to run on (default: runs local)
     """
    if time is None:
        time = dc.get_utc_now()

    _run_for_log(time, cwd, server)
    trims = dc.get_knobs_and_trims_from_log(
        os.path.join(cwd, const.get_extractor_log_filename())
    )
    if only_active:
        for key, value in trims.items():
            if value == 0:
                del trims[key]
    LOG.debug("Extracted Knob-Names: '{}'".format(",".join(trims.keys())))
    return trims.keys()


# Private Functions ############################################################


# Build command

def _run_knob_extraction(knob_names, time, cwd, server, ktype):
    cmd = " ".join([const.get_om_extractor(),
                    "-time", '"{:s}"'.format(dc.convert_utc_time_to_local(time)),
                    "-Knames", '"{:s}"'.format(",".join(knob_names)),
                    "-Ktype", const.get_ktype(ktype)])
    _run_command(cmd, server, cwd)


def _run_overview_extraction(knob_names, time, cwd, server):
    cmd = " ".join([const.get_om_extractor(),
                    "-e",  # extract optics
                    "-oe",  # extract orbit
                    "-time", '"{:s}"'.format(dc.convert_utc_time_to_local(time)),
                    "-Knames", '"{:s}"'.format(",".join(knob_names)),
                    "-Ktype", const.get_ktype("trim")])
    _run_command(cmd, server, cwd)


def _run_for_log(time, cwd, server):
    cmd = " ".join([const.get_om_extractor(),
                    "-time", '"{:s}"'.format(dc.convert_utc_time_to_local(time)),
                    "-oe",  # extract orbit
                    "-Knames", "dummy",  # makes it crash, but knobs will be in logs
                    "-Ktype", const.get_ktype("trim")
                    ])
    _run_command(cmd, server, cwd)


# Run command

def _run_command(cmd, server, cwd):
    if server:
        cmd = "ssh -t -X -o 'StrictHostKeyChecking no' {:s} 'cd {:s} ; {:s}'".format(
            server, cwd, cmd)
    LOG.debug("Running command:\n{}".format(cmd))
    process = subprocess.Popen(cmd, shell=True, cwd=cwd)
    process.wait()


# Output

def _log_df(df):
    LOG.debug("Extracted data summary:")
    s = "{:<20s}: {:s}"
    k = "  {:<30s} = {:g}"
    with logging_tools.unformatted_console_logging():
        LOG.info("")
        LOG.info(s.format("Time (UTC)", df.headers[const.get_time_header()]))
        LOG.info(s.format("Beamprocess", df.headers[const.get_beamprocess_header()]))
        LOG.info(s.format("Optics", df.headers[const.get_optics_header()]))
        LOG.info(s.format("Fill", str(df.headers[const.get_fill_header()])))
        LOG.info(s.format("Trims", ""))
        for idx in df.index:
            LOG.info(k.format(*df.loc[idx, :].values))
        LOG.info("")


def _plot_orbits(cwd, show_plot):
    beams = [1, 2]
    filenames = [os.path.join(cwd, const.get_default_orbit_filename(b)) for b in beams]
    output = os.path.splitext(filenames[0])[0]
    file_labels = ["Beam {:d}".format(b) for b in beams]

    plot(files=filenames,
         y_cols=["X", "Y"],
         y_labels=["Orbit [m]"] * 2,
         file_labels=file_labels,
         figure_per_file=True,
         output=output,
         no_show=not show_plot,
         )
