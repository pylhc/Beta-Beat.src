from os import listdir
from os.path import isfile, join
from utils import tfs_pandas, stats
from tests.test_utils.twiss_to_lin import optics_measurement_test_files

#def test_measure_optics():
#    create_model()
#    files_to_load = optics_measurement_test_files(modeldir, source_outpath)
#    run_measure_optics()
#    test_output(meas_path)


def test_output(meas_path):
    for f in _find_measurement_files(meas_path):
        try:
            a = tfs_pandas.read_tfs(join(meas_path, f))
        except:
            continue
        cols = [column for column in a.columns.values if column.startswith('DELTA')]
        for col in cols:
            rms = stats.weighted_rms(a.loc[:, col].values)
            print("File: {}  Column: {}   RMS: {}".format(f, col, rms))


def _find_measurement_files(meas_path):
    return [f for f in listdir(meas_path) if isfile(join(meas_path, f))]