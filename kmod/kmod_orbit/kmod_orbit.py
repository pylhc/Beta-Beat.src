from __future__ import print_function
import sys
import os
import numpy as np
np.seterr(all="raise")
import time
import datetime
import argparse
sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
)))
from utils import tfs_pandas


BEAMS = BEAM1, BEAM2 = (0, 1)
SIDES = LEFT, RIGHT = (0, 1)
PLANES = HOR, VER = (0, 1)

PLANE_STR = {HOR: "X", VER: "Y"}
BEAM_STR = {BEAM1: "1", BEAM2: "2"}
SIDE_STR = {LEFT: "L", RIGHT: "R"}

CUTOFF = 0.15  # BPMs we ignor
# CUTOFF = 1.5 * np.max(v[0,:])


def _parse_args():
    parser = argparse.ArgumentParser(description="")  # TODO

    parser.add_argument(
        "--modelb1",
        help="Model path for beam 1 (with drifts).",
        required=True,
        dest="model1_path"
    )
    parser.add_argument(
        "--modelb2",
        help="Model path for beam 2 (with drifts).",
        required=True,
        dest="model2_path"
    )
    parser.add_argument(
        "--orbitleft",
        help="Path to the directory containing left quad orbits",
        required=True,
        dest="orbit_path_left"
    )
    parser.add_argument(
        "--orbitright",
        help="Path to the directory containing right quad orbits",
        required=True,
        dest="orbit_path_right"
    )
    parser.add_argument(
        "--kmodksleft",
        help="Path to the kmod output file containt the left quad K samples",
        required=True,
        dest="kleft_path"
    )
    parser.add_argument(
        "--kmodksright",
        help="Path to the kmod output file containt the right quad K samples",
        required=True,
        dest="kright_path"
    )
    parser.add_argument(
        "--ip",
        help="Interaction point to use.",
        required=True,
        type=int,
        dest="ip",
    )
    options = parser.parse_args()
    input = (options.model1_path, options.model2_path,
             options.orbit_path_left, options.orbit_path_right,
             options.kleft_path, options.kright_path,
             options.ip)
    return input


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
            for plane in (HOR, VER):
                closest = np.argmin(np.abs(timestamp - k_file.index.values))
                val = _get_sign(beam, plane) * k_file.K.iloc[closest]
                k[(beam, side, plane)].append(val)
                orbit = orbit_data.loc[:, PLANE_STR[plane]]
                bpm[(beam, side, plane)].append(orbit - np.mean(orbit))
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
    this_k = _compute_kl(ks[(beam, side, plane)], this_model, quadname)

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


def _compute_kl(ks, model, quadname):
    quad_length = (model.loc[quadname, "S"] -
                   model.ix[model.index.get_loc(quadname) - 1, "S"])
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
    return {(beam, side, plane): function(beam, side, plane)
            for beam in BEAMS for side in SIDES for plane in PLANES}


def _get_sign(beam, plane):
    return {(BEAM1, HOR): 1,
            (BEAM1, VER): -1,
            (BEAM2, HOR): -1,
            (BEAM2, VER): 1, }[(beam, plane)]


def _date_str_to_timestamp(str_date):
    timestamp = time.mktime(datetime.datetime.strptime(
        str_date,
        "%Y-%m-%d %H:%M:%S.%f"
    ).timetuple())
    timestamp = timestamp * 1000
    return long(timestamp)


def _terminal_start():
    _input = _parse_args()
    results = compute_offset(*_input)
    _apply_to_beam_side_plane(
        lambda beam, side, plane: print(
            "Beam" + BEAM_STR[beam],
            "Side " + SIDE_STR[side],
            "Plane " + PLANE_STR[plane] + ":",
            str(results[(beam, side, plane)])
        )
    )


if __name__ == "__main__":
    _terminal_start()
