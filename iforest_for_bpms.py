# coding: utf-8
'''
Created on May 4, 2018
Application of Isolation Forest algorithm for the diagnostics of BPM signal.
Original paper: Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. “Isolation forest.” Data Mining, 2008. ICDM‘08. Eighth IEEE International Conference.

@author: Elena Fol
'''
from __future__ import print_function
import os
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm
import pandas
import argparse
from sklearn.ensemble import IsolationForest
from utils import tfs_pandas, logging_tools
from model.accelerators import lhc


LOGGER = logging_tools.get_logger(__name__)
ARCS_CONT = 0.01
IRS_CONT = 0.01
FEATURES = "TUNE{0},NOISE_SCALED,AMP{0}"
FEATURES_WITH_NAME = "NAME,TUNE{0},NOISE_SCALED,AMP{0}"
PLANE = ("x", "y")


def get_bad_bpms(files, remove_bpms):
    files_list = files.split(',')
    files_x, files_y = separate_by_plane(files_list)
    bad_dfs = {"x": [], "y": []}
    bad_bpms_to_write = {"x": None, "y": None}
    for files, plane in ((files_x, "x"), (files_y, "y")):
        bad_dfs[plane].append(get_bad_bpms_from_measurement(files, plane))
    bad_bpms_to_write["x"] = pandas.concat(bad_dfs["x"])
    bad_bpms_to_write["y"] = pandas.concat(bad_dfs["y"])
    write_bad_bpms(files_list[0], bad_bpms_to_write)
    if remove_bpms:
        remove_bpms_from_file(files_x, set(bad_bpms_to_write["x"].NAME), "x")
        remove_bpms_from_file(files_y, set(bad_bpms_to_write["y"].NAME), "y")


def write_bad_bpms(first_file, bad_bpms_to_write):
    meas_dir = os.path.abspath(os.path.join(first_file, os.pardir))
    for plane in PLANE:
        bad_bpms_summary_path = os.path.join(meas_dir, "bad_bpms_iforest_{}.tfs".format(plane))
        tfs_pandas.write_tfs(bad_bpms_summary_path, bad_bpms_to_write[plane])
    LOGGER.info("Bad BPMs summary from Isolation Forest written to: %s", meas_dir)


def get_bad_bpms_from_measurement(files, plane):
    uplane = plane.upper()
    bpm_tfs_data = _create_tfs_data(files, plane)
    arc_bpm_data, ir_bpm_data = get_data_for_clustering(bpm_tfs_data, plane)
    dataframes = []
    for data_for_clustering, cont, title in ((arc_bpm_data, ARCS_CONT, "Arcs " + uplane),
                                      (ir_bpm_data, IRS_CONT, "IRs " + uplane)):
        bad_bpms, good_bpms, all_bpms_scores, bad_bpms_scores =\
             detect_anomalies(cont, data_for_clustering, uplane)
        bpm_tfs_data, data_for_clustering, bad_bpms, good_bpms =\
            [reassign_index(data)
             for data in (bpm_tfs_data, data_for_clustering, bad_bpms, good_bpms)]
        signif_feature = get_significant_features(bpm_tfs_data, data_for_clustering, bad_bpms, good_bpms, plane)
        signif_feature.loc[:, "SCORE"] = bad_bpms_scores
        dataframes.append(signif_feature)
        if plot:
            plot_scores_threshold(all_bpms_scores, cont, title)
            plot_bpms_3d(good_bpms, bad_bpms, uplane, title, bad_bpms_scores)
    return pandas.concat(dataframes)


def reassign_index(data):
    data["NEW_INDEX"] = range(len(data.NAME))
    return data.set_index("NEW_INDEX")


def get_significant_features(bpm_tfs_data, data_for_clustering, bad_bpms, good_bpms, plane):
    features_df = pandas.DataFrame(index=bad_bpms.index)
    columns = FEATURES.format(plane.upper()).split(",")
    for index in bad_bpms.index:
        max_dist = max([(abs(data_for_clustering.loc[index, col] -
                             good_bpms.loc[:, col].mean()), col)
                        for col in columns])
        max_dist, sig_col = max_dist
        features_df.loc[index, "NAME"] = bad_bpms.loc[index, "NAME"]
        features_df.loc[index, "FEATURE"] = sig_col
        features_df.loc[index, "VALUE"] = bpm_tfs_data.loc[index, sig_col]
        features_df.loc[index, "AVG"] = np.mean(bpm_tfs_data.loc[good_bpms.index][sig_col])
    return features_df


def detect_anomalies(contamination, data, uplane):
    iforest = IsolationForest(n_estimators=100, max_samples='auto',
                              contamination=contamination, max_features=1.0,
                              bootstrap=False)
    features = data[FEATURES.format(uplane).split(",")]
    iforest.fit(features)
    labels = iforest.predict(features)
    bad_bpms = data.iloc[np.where(labels == -1)].copy()
    good_bpms = data.iloc[np.where(labels != -1)].copy()
    all_bpms_scores = iforest.decision_function(features)
    bad_bpms_scores = iforest.decision_function(features.iloc[np.where(labels == -1)])
    return bad_bpms, good_bpms, all_bpms_scores, bad_bpms_scores


def get_data_for_clustering(bpm_tfs_data, plane):
    columns = FEATURES.format(plane.upper()).split(",")
    ir_bpm_data_for_clustering = bpm_tfs_data.iloc[~lhc.Lhc.get_arc_bpms_mask(bpm_tfs_data.NAME)].copy()
    arc_bpm_data_for_clustering = bpm_tfs_data.iloc[lhc.Lhc.get_arc_bpms_mask(bpm_tfs_data.NAME)].copy()
    for col in columns:
        ir_bpm_data_for_clustering.loc[:, col] = _normalize_parameter(ir_bpm_data_for_clustering.loc[:, col])
        arc_bpm_data_for_clustering.loc[:, col] = _normalize_parameter(arc_bpm_data_for_clustering.loc[:, col])
    return arc_bpm_data_for_clustering, ir_bpm_data_for_clustering


def _create_tfs_data(filepaths, plane):
    bpm_data_rows = []
    for filepath in filepaths:
        bpm_tfs_file = tfs_pandas.read_tfs(filepath)
        bpms_tfs_data = bpm_tfs_file[FEATURES_WITH_NAME.format(plane.upper()).split(",")]
        bpm_data_rows.append(bpms_tfs_data)
    return pandas.concat(bpm_data_rows)


def _normalize_parameter(column_data):
    return (column_data - column_data.min()) / (column_data.max() - column_data.min())


def separate_by_plane(files_list):
    files_x = []
    files_y = []
    for file_in in files_list:
        if os.path.basename(file_in).endswith("linx"):
            files_x.append(file_in)
        elif os.path.basename(file_in).endswith("liny"):
            files_y.append(file_in)
        else:
            print("Given file is not a measurement!")
    return files_x, files_y


def remove_bpms_from_file(paths, bad_bpm_names, plane):
    """
    Writes a backup of the original .lin files (e.g .linx --> .linx.notcleaned)
    and removes the BPNs identified by iForest as bad.
    :param paths: original lin files
    :param bad_bpm_names: list of the names of bad BPMs identified by iForest
    """
    for path in paths:
        src_dir = os.path.abspath(os.path.join(path, os.pardir))
        filename = os.path.basename(path)
        new_filename = os.path.join(src_dir, filename + ".notcleaned")
        os.rename(path, new_filename)
        original_file_tfs = tfs_pandas.read_tfs(new_filename).set_index("NAME", drop=False)
        original_file_tfs = original_file_tfs.loc[~original_file_tfs.index.isin(bad_bpm_names)]
        pln_num = "1" if plane == "x" else "2"
        original_file_tfs.headers["Q{}".format(pln_num)] =\
            original_file_tfs["TUNE" + plane.upper()].mean()
        original_file_tfs.headers["Q{}RMS".format(pln_num)] =\
            np.std(original_file_tfs["TUNE" + plane.upper()])
        tfs_pandas.write_tfs(path, original_file_tfs, original_file_tfs.headers)


def revert_forest_cleaning(files):
    """
    Reverts the cleaning. The backup files are renamed back to the original names (e.g .linx.notcleaned --> .linx)
    :param paths: list of files, where bad bpms identified by iForest are removed
    """
    files_list = files.split(',')
    for path in files_list:
        src_dir = os.path.abspath(os.path.join(path, os.pardir))
        filename = os.path.basename(path)
        notcleaned_file = os.path.join(src_dir, filename + ".notcleaned")
        original_file_tfs = tfs_pandas.read_tfs(notcleaned_file).set_index("NAME", drop=False)
        os.remove(path)
        lin_file = os.path.join(src_dir, notcleaned_file.replace(".notcleaned",""))
        os.rename(notcleaned_file, lin_file)
        tfs_pandas.write_tfs(lin_file, original_file_tfs)


def plot_scores_threshold(scores, cont, title):
    threshold = stats.scoreatpercentile(scores,100*ARCS_CONT)
    plt.hist(scores, bins=100, range=(-0.5,0.5), edgecolor='black', linewidth=1)
    plt.axvline(x=threshold, color='r', linestyle='-', label='Learned Threshold')
    plt.text(threshold, 0, str(threshold)[:-10], color='r', fontsize=8, verticalalignment='bottom', horizontalalignment='right')
    plt.title(title)
    plt.xlabel("Anomaly score (the lower, the more abnormal)")
    plt.ylabel("Number of BPMs")
    plt.legend()
    plt.show()
    

def plot_bpms_3d(good_bpms, bad_bpms, plane, title, scores):
    columns = FEATURES_WITH_NAME.format(plane.upper()).split(",")
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.plot3D(good_bpms.loc[:, columns[1]], good_bpms.loc[:, columns[2]], good_bpms.loc[:, columns[3]], 'o', markerfacecolor="black",
                      markeredgecolor='black', markersize=14, label = "good")
    ax.plot3D(bad_bpms.loc[:, columns[1]], bad_bpms.loc[:, columns[2]], bad_bpms.loc[:, columns[3]], '^', markerfacecolor="red",
                      markeredgecolor='black', markersize=14, label = "faulty")
    for index, score in zip(bad_bpms.index, scores):
        ax.text(bad_bpms.loc[index, columns[1]], bad_bpms.loc[index, columns[2]], bad_bpms.loc[index, columns[3]], bad_bpms.loc[index,"NAME"] + " {" + str(score)[:-10] + "}")
    ax.set_xlabel('Tune', fontsize = 25, linespacing=3.2)
    ax.set_ylabel('Amplitude', fontsize = 25, linespacing=3.2)
    ax.set_zlabel('Noise', fontsize = 25, linespacing=3.2)
    for axis in ('x', 'y', 'z'):
        ax.tick_params(axis=axis, labelsize=15)
    plt.legend(fontsize = 25)
    plt.title(title)
    plt.show()


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--files",
        dest="files", type=str,
    )
    parser.add_argument(
        "--remove_bpms",
        dest="remove_bpms",
        action="store_true",
    )
    parser.add_argument(
        "--revert",
        dest="revert",
        action="store_true",
    )
    parser.add_argument(
        "--plot",
        dest="plot",
        action="store_true",
    )
    options = parser.parse_args()
    return options.files, options.remove_bpms, options.revert, options.plot


if __name__ == '__main__':
    _files, _remove_bpms, _revert, plot = _parse_args()
    if(not _revert):
        get_bad_bpms(_files, _remove_bpms)
    else:
        revert_forest_cleaning(_files)
