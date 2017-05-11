import os
import sys

sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    ".."
)))
from Utilities import outliers
from Utilities import tfs_pandas as tfs


SUFFIXES = {'x' : ".linx", 'y': ".liny"}
PLANES=['x','y']
TUNES = {'x' : "TUNEX", 'y' : "TUNEY"}


def clean_tunes(filename):
    for plane in PLANES:
        file_name=filename + SUFFIXES[plane]
        f = tfs.read_tfs(file_name)
        mask = outliers.get_filter_mask(f.loc[:,TUNES[plane]], limit=0.00001)
        f = f.loc[mask,:]
        tfs.write_tfs(f,f.headers,file_name + ".new")
        os.rename(file_name, file_name + ".raw")
        os.rename(file_name + ".new", file_name) 
        

if __name__ == "__main__":
    clean_tunes(*sys.argv[1:])
