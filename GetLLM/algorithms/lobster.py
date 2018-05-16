"""
..module: lobster

"""
import numpy as np
import pandas as pd

from beta import tilt_slice_matrix
from utils import logging_tools

LOGGER = logging_tools.get_logger(__name__)

NEARPI = 1.0e-3
NEARDIST = 6  # maximal separation of BPMs until next pi

def get_local_observable(phase_d, free_model, files_dict):
    """Calculates local observable for pi separated BPMs
    """

    LOGGER.info("Computing local observable")

    modl_phases = phase_d.phase_advances_free_x["MODEL"]
    meas_phases = phase_d.phase_advances_free_x["MEAS"]
    err_phases = phase_d.phase_advances_free_x["ERRMEAS"]

    diff_pi = modl_phases.values - .5
    near_pi = abs(diff_pi) < NEARPI

    tfs_file = files_dict["getlobster.out"]
    tfs_file.add_column_names(['NAME', 'S', 'NEXT', 'DELTAPHI_NEAR', 'ERRPHI_NEAR', 'MDLPHI_NEAR'])
    tfs_file.add_column_datatypes(['%s', "%le", "%s", "%le", "%le", "%le"])

    data = np.empty((len(modl_phases.columns)),
                    dtype=[('NAME', 'S64'),
                           ('S', float),
                           ('NEXT', 'S64'),
                           ('DELTAPHI', float)])

    collected = []
    for i, column in enumerate(modl_phases.columns):
        # search only downstream
        mask = near_pi[i]
        mask[:i] = False
        mask[i + NEARDIST:] = False

        modl_column = modl_phases[column][mask]
        meas_column = meas_phases[column][mask]
        erro_column = err_phases[column][mask]
        nextname = "---"
        delta_phi = 0
        err = 0
        s = free_model.loc[column, "S"]

        if len(meas_column) > 0:
            mdl_phi = modl_column[0]
            delta_phi = meas_column[0] - mdl_phi
            if delta_phi > 1.0:
                delta_phi -= 2.0
            if delta_phi < -1.0:
                delta_phi += 2.0
            nextname = meas_column.index[0]
            err = erro_column[0]

            collected.append(column)
            tfs_file.add_table_row([column, s, nextname, delta_phi, err, mdl_phi])

        data[i] = (column, s, nextname, delta_phi)

    df = pd.DataFrame(data, index=modl_phases.columns)

    LOGGER.info("X: collected local observables: {}".format(len(collected)))

    return df.loc[collected]

