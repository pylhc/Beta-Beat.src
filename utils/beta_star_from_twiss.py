"""
Module utils.beta_star_from_twiss
------------------------------------

Functions to calculate beta* and waist information from twiss files.


.. rubric:: References

..  [#CarlierAccuracyfeasibilitybeta2017]
    F. Carlier et al.,
    Accuracy and feasibility of the beta* measurement for LHC
    and High Luminosity LHC using k modulation,
    Phys. Rev. Accel. Beams, vol. 20, no. 1, p. 011005, Jan. 2017.
    doi: [10.1103/PhysRevAccelBeams.20.011005](https://doi.org/10.1103/PhysRevAccelBeams.20.011005)
"""

from tfs_files import tfs_pandas as tfs

POSITION_COLUMN = "S"
BETA_COLUMN = "BET{plane:s}"
ALPHA_COLUMN = "ALF{plane:s}"
PLANES = "XY"
RES_COLUMNS = ["LABEL", "S", "BETASTAR", "WAIST", "BETAWAIST"]
RESULTS_PREFIX = "betaip_"
RESULTS_SUFFIX = ".tfs"


def get_beta_star_and_waist_from_ip(tfs_df, beam, locations, labels=None):
    """ Get beta* and waist* info directly from IP.

    Args:
        tfs: pandas dataframe or path to tfs file
        locations: ip-labels to get the info from
        labels: labels to use (default: location names)

    Returns:
         table with beta* and waist info
    """

    try:
        tfs_df = tfs.read_tfs(tfs_df, index="NAME")
    except TypeError:
        pass

    if labels and len(labels) != len(locations):
        raise ValueError("Length of locations and labels need to be equal.")

    if not labels:
        labels = locations

    tfs_out = tfs.TfsDataFrame(columns=RES_COLUMNS).set_index(RES_COLUMNS[0])
    for loc, lab in zip(locations, labels):
        for plane in PLANES:
            pos = tfs_df.loc[loc, POSITION_COLUMN]
            beta_star = tfs_df.loc[loc, BETA_COLUMN.format(plane=plane)]
            alpha_star = tfs_df.loc[loc, ALPHA_COLUMN.format(plane=plane)]
            beta_waist = beta_star / (1 + alpha_star)
            waist = alpha_star * beta_waist
            tfs_out.loc[get_full_label(lab, beam, plane)] = [pos, beta_star, waist, beta_waist]
    return tfs_out


def get_full_label(label, beam, plane):
    """ Full label name """
    return "{:s}.B{:d}.{:s}".format(label, beam, plane)



