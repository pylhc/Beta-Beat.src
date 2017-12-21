from Utilities import tfs_pandas as tfs
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    rdt = "f0040_AMP"
    twissfile = 'twiss_2oct_elements.dat'
    ptcfile = 'ptc_normal_2octupoles.dat'

    twiss = tfs.read_tfs(twissfile).set_index('NAME')
    ptc = tfs.read_tfs(ptcfile).set_index('NAME')
    idx_rdt = np.where(ptc.columns == rdt)[0][0]

    d = ptc.iloc[1:, idx_rdt] - ptc.iloc[:-1, idx_rdt].values
    bpms_end = d.index[np.abs(d) > .1].values
    bpms_start = d.index[[np.where(d.index == b)[0][0] - 1 for b in bpms_end]]

    [twiss.loc[bpms_start[i]:bpms_end[i], ['K1L', 'K2L', 'K3L']].sum() for i in range(len(bpms_start))]

    cm = plt.get_cmap('tab20c')
    n = len(bpms_start)

    ax = ptc.plot(x="S", y=rdt)
    y = [-1000, 1000]
    [ax.plot([ptc.loc[b, "S"], ptc.loc[b, "S"]], y, linestyle='--', color=cm(1.*i/n))
     for i, b in enumerate(bpms_start)]
    [ax.text(ptc.loc[b, "S"], 550, b, color=cm(1.*i/n))
     for i, b in enumerate(bpms_start)]

    [ax.plot([ptc.loc[b, "S"], ptc.loc[b, "S"]], y, linestyle='-', color=cm(1.*i/n))
     for i, b in enumerate(bpms_end)]
    [ax.text(ptc.loc[b, "S"], 500, b, color=cm(1.*i/n))
     for i, b in enumerate(bpms_end)]
    plt.show()
