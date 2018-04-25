from __future__ import print_function
import os
import sys
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)))
from utils.contexts import temporary_dir, silence
import madx_wrapper


def test_with_macro():
    """ Checks:
         - Output_file is created.
         - Macros resolve correctly.
    """
    content = (
        "!@require lhc_runII\n"
        "!@require lhc_runII_ats.macros.madx\n"
    )
    resolved_lines = [
        'call,file="{}";'.format(os.path.join(madx_wrapper.LIB, "lhc_runII.macros.madx")),
        'call,file="{}";'.format(os.path.join(madx_wrapper.LIB, "lhc_runII_ats.macros.madx"))
    ]

    with temporary_dir() as tmpdir:
        outfile = os.path.join(tmpdir, "job.with_macro.madx")
        with silence():
            madx_wrapper.resolve_and_run_string(content, output_file=outfile, cwd=tmpdir)
        assert os.path.exists(outfile)
        with open(outfile, "r") as of:
            out_lines = of.read().split("\n")
        out_lines = [ol.replace(" ", "") for ol in out_lines]
        assert all([r in out_lines for r in resolved_lines])


def test_with_nonexistent_file():
    """ Checks:
         - Madx crashes when tries to call a non-existent file
         - Logfile is created
         - Error message is read from log
    """
    call_file = "does_not_exist.madx"
    content = "call, file ='{:s}';".format(call_file)
    with temporary_dir() as tmpdir:
        log_file = os.path.join(tmpdir, "tmp_log.log")
        with pytest.raises(madx_wrapper.MadxError) as e:
            madx_wrapper.resolve_and_run_string(content, log_file=log_file, cwd=tmpdir)
        assert os.path.isfile(log_file)
        assert call_file in str(e)


if __name__ == "__main__":
    # test_with_macro()
    # test_with_nonexistent_file()
    pass