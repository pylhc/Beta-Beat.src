"""
..module: lobster

"""
import numpy as np
from numpy import sin
import pandas as pd

from beta import tilt_slice_matrix
from utils import logging_tools

LOGGER = logging_tools.get_logger(__name__)

NEARPI = 4.0e-4
NEARDIST = 6  # maximal separation of BPMs until next pi
ROB = 7  # range of BPMs

def get_local_observable(phase_d, free_model, tune):
    """Calculates local observable for pi separated BPMs
    """

    LOGGER.info("Computing local observable")

    modl_phases = phase_d.phase_advances_free_x["MODEL"]
    meas_phases = phase_d.phase_advances_free_x["MEAS"]
    err_phases = phase_d.phase_advances_free_x["ERRMEAS"]

    diff_pi = modl_phases.values - .5
    near_pi = abs(diff_pi) < NEARPI

    data = np.empty((len(modl_phases.columns)),
                    dtype=[('NAME', 'S64'),
                           ('S', float),
                           ('NEXT', 'S64'),
                           ('DELTAPHI_NEAR', float)])

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
         #   if delta_phi > 1.0:
         #       delta_phi -= 2.0
         #   if delta_phi < -1.0:
         #       delta_phi += 2.0
            nextname = meas_column.index[0]
            err = erro_column[0]

            collected.append(column)

        data[i] = (column, s, nextname, delta_phi)

    df = pd.DataFrame(data, index=modl_phases.columns)
                      #columns=["NAME", "S", "NEXT", "DELTAPHI"])

    LOGGER.info("collected local observables: {}".format(len(collected)))

    # regional
    LOGGER.info("calculating regional observable")
    LOGGER.debug("range of BPMs: {}".format(ROB))

    combos = [(x,y,z,w)
              for x in range(1, ROB)
              for y in range(x + 1, ROB)
              for z in range(y + 1, ROB)
              for w in range(z + 1, ROB)]

    tilted_meas = (tilt_slice_matrix(meas_phases.values, 0, len(meas_phases), tune))
    tilted_modl = (tilt_slice_matrix(modl_phases.values, 0, len(modl_phases), tune))

    reg_obs = pd.DataFrame(index=modl_phases.index,
                               columns=[get_comb_pattern(co) for co in combos])

    for (ix, iy, iz, iw) in combos:

        patt = get_comb_pattern((ix, iy, iz, iw))
        sinii = sin([
        tilted_modl[iz] - tilted_modl[ix],
        tilted_modl[iy] - tilted_modl[ix],
        tilted_modl[iw] - tilted_modl[ix],
        tilted_modl[iy] - tilted_modl[ix]]) / sin([
            tilted_modl[iz] - tilted_modl[iy],
            tilted_modl[iz] - tilted_modl[iy],
            tilted_modl[iw] - tilted_modl[iy],
            tilted_modl[iw] - tilted_modl[iy]]) / sin([
                tilted_modl[iy] - tilted_modl[ix],
                tilted_modl[iz] - tilted_modl[ix],
                tilted_modl[iy] - tilted_modl[ix],
                tilted_modl[iw] - tilted_modl[ix]])

        phadvbeat = np.array([
                tilted_meas[iy] - tilted_meas[ix],
                -tilted_meas[iz] + tilted_meas[ix],
                -tilted_meas[iy] + tilted_meas[ix],
                tilted_meas[iw] - tilted_meas[ix]])

        regobsmatr = np.dot(np.transpose(sinii), phadvbeat)
        LOGGER.info(regobsmatr.shape)
        reg_obs.loc[:, patt] = np.sum(regobsmatr, axis=1)

    reg_obs.loc[:, "S"] = free_model.loc[modl_phases.index, "S"]

    return df.loc[collected], reg_obs

def get_comb_pattern(comb):
    """Returns a string representing the pattern of the combination
    ex.: (1,2,4,5) --> -ij-kl
    """
    combstring = ["-"] * (ROB - 1)
    combstring[comb[0]-1] = "i"
    combstring[comb[1]-1] = "j"
    combstring[comb[2]-1] = "k"
    combstring[comb[3]-1] = "l"
    return "".join(combstring)
