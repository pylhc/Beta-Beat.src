import sys
from os import listdir
from os.path import isfile, join, dirname, abspath

sys.path.append(abspath(join(dirname(__file__), "..", "..")))

from utils import tfs_pandas, stats, iotools, outliers
from tests.test_utils.twiss_to_lin import optics_measurement_test_files
import measure_optics
from optics_measurements import optics_input
from model import manager
import pytest


LIMITS = {'P': 1.5e-4, 'B': 1e-2, 'D': 1e-2}
DEFAULT_LIMIT = 5e-3


def test_measure_optics(_create_input):
    input_files, measure_input = _create_input
    measure_optics.measure_optics(input_files, measure_input)
    for f in _find_measurement_files(measure_input.outputdir):
        if not f.startswith("get"):
            continue
        a = tfs_pandas.read_tfs(join(measure_input.outputdir, f))
        cols = [column for column in a.columns.values if column.startswith('DELTA')]
        for col in cols:
            mask = outliers.get_filter_mask(a.loc[:, col].values, limit=0.1)
            rms = stats.weighted_rms(a.loc[mask, col].values)
            if col[5] in LIMITS.keys():
                assert rms < LIMITS[col[5]], "\nFile: {:25}  Column: {:15}   RMS: {:.6f}".format(f, col, rms)
            else:
                assert rms < DEFAULT_LIMIT, "\nFile: {:25}  Column: {:15}   RMS: {:.6f}".format(f, col, rms)
            print("\nFile: {:25}  Column: {:15}   RMS:    {:.6f}".format(f, col, rms))


@pytest.fixture()
def _create_input():
    sourcespath = abspath(join(dirname(__file__), "..", "measurements"))
    iotools.create_dirs(sourcespath)
    accelerator = {"model_dir": join(dirname(__file__), "..", "inputs", "models", "25cm_beam1"),
                   "accel": "lhc", "lhc_mode": "lhc_runII_2018", "beam": 1}
    files_to_load = optics_measurement_test_files(accelerator["model_dir"], sourcespath)
    measure_input = optics_input.OpticsInput()
    measure_input.outputdir = abspath(join(dirname(__file__), "..", "results"))
    measure_input.accelerator = manager.get_accel_instance(accelerator)
    input_files = measure_optics.InputFiles(files_to_load)
    return input_files, measure_input


def _find_measurement_files(meas_path):
    return [f for f in listdir(meas_path) if isfile(join(meas_path, f))]

