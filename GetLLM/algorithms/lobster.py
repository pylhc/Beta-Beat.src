import numpy as np
import pandas as pd

from beta import tilt_slice_matrix

NEARPI = 1.0e-3

def get_local_observable(phase_d, getllm_d, files_dict):
    """Calculates local observable for pi separated BPMs
    """

    free_model = getllm_d.accelerator.get_model_tfs()
    modl_phases = phase_d.phase_advances_free_x["MODEL"]
    meas_phases = phase_d.phase_advances_free_x["MEAS"]
    err_phases = phase_d.phase_advances_free_x["ERRMEAS"]

    diff_pi = modl_phases.values - 1
    near_pi = abs(diff_pi) < NEARPI

    tfs_file = files_dict["getlobster.out"]
    tfs_file.add_column_names(['NAME', 'S', 'NEXT', 'DELTAPHI', 'ERRPHI', 'MDLPHI'])
    tfs_file.add_column_datatypes(['%s', "%le", "%s", "%le", "%le", "%le"])

    data = np.empty((len(modl_phases.columns)),
                    dtype=[('NAME', 'S64'),
                           ('S', float),
                           ('NEXT', 'S64'),
                           ('DELTAPHI', float)])

    for i, column in enumerate(modl_phases.columns):
        mask = near_pi[:, i]
        mask[:i] = False
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
            if delta_phi > .5:
                delta_phi -= 1.0
            nextname = meas_column.index[0]
            err = erro_column[0]

            tfs_file.add_table_row([column, s, nextname, delta_phi, err, mdl_phi])

        data[i] = (column, s, nextname, delta_phi)

    df = pd.DataFrame(data, index=modl_phases.columns)



