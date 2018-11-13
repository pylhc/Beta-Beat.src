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


def extract_knob_values(knob_name, time, cwd="./", server=None):
    """ Extract knob values from online model using the knob name and the time.

    Example Call::

        extract_knob_values("LHCBEAM/2018_global_ats_flat_b1_for_ip5_waist",
                            "2018-10-30 15:00:00.0",
                            cwd="/afs/cern.ch/user/j/jdilly/extractor/",
                            server="cs-ccr-dev3")

    Args:
        knob_name: Name of the knob
        time: UTC time in ISOformat to extract
        cwd: output directory for the results and the log (optional)
        server: server to run on (optional)
    """
    knob_file = os.path.join(cwd, const.get_default_knob_filename())

    # first get the trim value and madx_changefile
    _run_process(knob_name, time, cwd, server, "trim")
    trim, _ = dc.post_trim_extract(knob_file)

    # run again and get the values and deltas
    _run_process(knob_name, time, cwd, server, "k")
    df = dc.knob_k_to_tfs(knob_file, trim)
    tfs.write_tfs(os.path.join(cwd, knob_file.replace(".madx", ".tfs")), df)


# Private Functions ############################################################


def _run_process(knob_name, time, cwd, server, ktype):
    cmd = " ".join([const.get_om_extractor(),
                    "-time", '"{:s}"'.format(dc.convert_utc_time_to_local(time)),
                    "-Knames", '"{:s}"'.format(knob_name),
                    "-Ktype", const.get_ktype(ktype)])

    if server:
        cmd = "ssh -t -X -o 'StrictHostKeyChecking no' {:s} 'cd {:s} ; {:s}'".format(
            server, cwd, cmd)
    process = subprocess.Popen(cmd, shell=True, cwd=cwd)
    process.wait()


# Script Mode ##################################################################


if __name__ == '__main__':
    pass
