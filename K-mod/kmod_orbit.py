import sys
import os
import numpy as np
import time
import datetime
sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    ".."
)))
from Utilities import tfs_pandas


BEAMS = BEAM1, BEAM2 = (0, 1)
SIDES = LEFT, RIGHT = (0, 1)
PLANES = HOR, VER = (0, 1)

PLANE_STR = {HOR: "X", VER: "Y"}
BEAM_STR = {BEAM1: "1", BEAM2: "2"}
SIDE_STR = {LEFT: "L", RIGHT: "R"}

CUTOFF = 0.15  # BPMs we ignor
# CUTOFF = 1.5 * np.max(v[0,:])


def compute_offset(model1_path, model2_path,
                   orbit_path_left, orbit_path_right,
                   kleft_path, kright_path,
                   ip):

    ks, orbits = _collect_orbit_data(orbit_path_left, orbit_path_right,
                                     kleft_path, kright_path)

    model1 = tfs_pandas.read_tfs(model1_path)
    model1.set_index("NAME", inplace=True)
    model2 = tfs_pandas.read_tfs(model2_path)
    model2.set_index("NAME", inplace=True)

    bpm_names = _get_bpm_names(orbit_path_left, orbit_path_right)
    models = {BEAM1: model1, BEAM2: model2}

    results = _apply_to_beam_side_plane(
        lambda beam, side, plane: _compute_and_clean(
            ip, beam, side, plane, models, bpm_names, ks, orbits
        )
    )
    print(results)
    return results


def _collect_orbit_data(orbit_path_left, orbit_path_right,
                        kleft_path, kright_path):
    k_file_left = tfs_pandas.read_tfs(kleft_path)
    k_file_left["TIME"] = k_file_left["TIME"].astype(long)
    k_file_left.set_index("TIME", inplace=True)

    k_file_right = tfs_pandas.read_tfs(kright_path)
    k_file_right["TIME"] = k_file_right["TIME"].astype(long)
    k_file_right.set_index("TIME", inplace=True)

    orbit_paths = {LEFT: orbit_path_left, RIGHT: orbit_path_right}
    k_files = {LEFT: k_file_left, RIGHT: k_file_right}

    k = {}
    _apply_to_beam_side_plane(
        lambda beam, side, plane: k.__setitem__((beam, side, plane), [])
    )

    bpm = {}
    _apply_to_beam_side_plane(
        lambda beam, side, plane: bpm.__setitem__((beam, side, plane), [])
    )
    for side in SIDES:
        orbit_path = orbit_paths[side]
        k_file = k_files[side]
        for filename in os.listdir(orbit_path):
            orbit_data = tfs_pandas.read_tfs(os.path.join(
                orbit_path,
                filename
            ))
            timestamp = _date_str_to_timestamp(orbit_data.headers["OrbitDate"])
            if filename.startswith('LHCB1'):
                beam = BEAM1
            elif filename.startswith('LHCB2'):
                beam = BEAM2
            try:
                k[(beam, side, HOR)].append(_get_sign(beam, HOR) *
                                            k_file.loc[timestamp, "K"])
                k[(beam, side, VER)].append(_get_sign(beam, VER) *
                                            k_file.loc[timestamp, "K"])
                orbit_x = orbit_data.loc[:, "X"]
                bpm[(beam, side, HOR)].append(orbit_x - np.mean(orbit_x))
                orbit_y = orbit_data.loc[:, "Y"]
                bpm[(beam, side, VER)].append(orbit_y - np.mean(orbit_y))
            except KeyError:
                continue
    return k, bpm


def _get_bpm_names(orbit_path_left, orbit_path_right):
    return {
        (BEAM1, LEFT): tfs_pandas.read_tfs(
            os.path.join(orbit_path_left, "LHCB1.orbit1.tfs"))["NAME"],
        (BEAM1, RIGHT): tfs_pandas.read_tfs(
            os.path.join(orbit_path_right, "LHCB1.orbit1.tfs"))["NAME"],
        (BEAM2, LEFT): tfs_pandas.read_tfs(
            os.path.join(orbit_path_left, "LHCB2.orbit1.tfs"))["NAME"],
        (BEAM2, RIGHT): tfs_pandas.read_tfs(
            os.path.join(orbit_path_right, "LHCB2.orbit1.tfs"))["NAME"],
    }


def _compute_and_clean(ip, beam, side, plane, models, bpm_names, ks, orbits):
    this_orbit = orbits[(beam, side, plane)]
    this_model = models[beam]
    this_bpm_names = bpm_names[(beam, side)]
    this_bpm_model = this_model.loc[this_bpm_names, :]

    quadname = "MQXA.1" + SIDE_STR[side] + str(ip)
    this_k = _compute_kl(ks[(beam, side, plane)], this_model, quadname, plane)

    orb = np.array(this_orbit)

    # Work out transfer matrix for given beam and quad
    model_data = _compute_transfer_matrix(this_model, this_bpm_model,
                                          quadname, this_k, plane, beam)

    # SVD on model and data
    u, s, v = np.linalg.svd(model_data)
    ud, sd, vd = np.linalg.svd(orb, full_matrices=False)

    # Remove offset
    sd[0] = 0
    orb = np.dot(np.dot(ud, np.diag(sd)), vd)
    ud, sd, vd = np.linalg.svd(orb)

    # SVD clean on data

    while max(vd[0, :]) > CUTOFF or min(vd[0, :]) < -CUTOFF:
        for j in reversed(range(0, len(vd[0, :]))):
            if abs(vd[0, j]) > CUTOFF:
                orb = np.delete(orb, (j), axis=1)
                model_data = np.delete(model_data, (j), axis=1)
        u, s, v = np.linalg.svd(model_data)
        ud, sd, vd = np.linalg.svd(orb)

    # Work out offset using matrix multiplication
    offset = np.dot(np.dot(np.transpose(u[:, 0]), orb),
                    np.transpose(v[0, :])) / s[0]
    N = orb - offset * model_data
    errorbars = np.abs(np.dot(np.dot(np.transpose(u[:, 0]), N),
                              np.transpose(v[0, :])) / s[0])
    return offset, errorbars


def _compute_kl(ks, model, quadname, ip):
    quad_length = (model.loc[quadname, "S"] -
                   model.ix[model.index.get_loc(quadname) - 1, "S"]),

    avg_k = np.mean(ks)
    return (ks - avg_k) * quad_length


def _compute_transfer_matrix(model, bpm_model, quadname, ks, plane, beam):
    mu_m = bpm_model.loc[:, "MU" + PLANE_STR[plane]]
    b_m = bpm_model.loc[:, "BET" + PLANE_STR[plane]]

    mu_q = model.loc[quadname, "MU" + PLANE_STR[plane]]
    b_q = model.loc[quadname, "BET" + PLANE_STR[plane]]
    a_q = model.loc[quadname, "ALF" + PLANE_STR[plane]]

    tune = model.headers["Q" + BEAM_STR[beam]]

    m11 = (np.sqrt(b_m / b_q) *
           (np.cos((tune - np.abs(mu_m - mu_q)) * 2 * np.pi) +
            a_q * np.sin((tune - np.abs(mu_m - mu_q)) * 2 * np.pi)))
    m12 = (np.sqrt(b_m * b_q) *
           np.sin((tune - np.abs(mu_m - mu_q)) * 2 * np.pi))

    # F = change in orbit
    F = -(b_q * ks / np.tan(tune * np.pi) /
          (2 + b_q * ks / np.tan(tune * np.pi)))

    # G = kick for each k value in the data (for 1m misalignment)
    G = (-(1 - a_q / np.tan(tune * np.pi)) * ks /
         (2 + b_q * ks / np.tan(tune * np.pi)))

    # Work out model data at each bpm for 1m misalignment
    return np.outer(F, m11) + np.outer(G, m12)


def _apply_to_beam_side_plane(function):
    results = {}
    for beam in BEAMS:
        for side in SIDES:
            for plane in PLANES:
                results[(beam, side, plane)] = function(beam, side, plane)
    return results


# TODO: Replace with dict
def _get_sign(beam, plane):
    if beam == BEAM1 and plane == HOR:
        return 1
    elif beam == BEAM1 and plane == VER:
        return -1
    elif beam == BEAM2 and plane == HOR:
        return -1
    elif beam == BEAM2 and plane == VER:
        return 1


def _date_str_to_timestamp(str_date):
    timestamp = time.mktime(datetime.datetime.strptime(
        str_date,
        "%Y-%m-%d %H:%M:%S.%f"
    ).timetuple())
    timestamp = timestamp * 1000
    return long(timestamp)
