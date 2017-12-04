from __future__ import print_function
import sys
import json
import pandas as pd
from collections import OrderedDict

sys.path.append("/afs/cern.ch/work/j/jcoellod/public/Beta-Beat.src")

from Utilities import tfs_pandas


def main(elems_file, sequences, output):
    sequences = sequences.split(",")
    elems_data = tfs_pandas.read_tfs(elems_file)
    s = elems_data.S
    names = elems_data.NAME
    length = elems_data.L

    elem_to_vars = get_elem_to_vars(names, sequences)
    variables = []
    for name in names:
        variables.append("".join(elem_to_vars[name]))

    all_data = OrderedDict((
        ("S", s),
        ("L", length),
        ("NAME", names),
        ("VARS", variables),
    ))
    data_frame = pd.DataFrame(
        data=all_data,
        columns=all_data.keys(),
    ).sort_values(by="S")
    tfs_pandas.write_tfs(output, data_frame)


def get_elem_to_vars(elems, sequences):
    elem_to_vars = OrderedDict({})
    for sequence in sequences:
        with open(sequence, "r") as seq_data:
            for line in seq_data:
                for elem in elems:
                    if elem in line and ":=" in line and "Kmax" not in line:
                        parts = line.split()
                        elem_name = parts[0].replace(",", "")
                        for part in parts:
                            if "kq" in part:
                                var = part
                                break
                        print(var)
                        try:
                            elem_name = elem_name[:elem_name.index(":")]
                        except ValueError:
                            pass
                        if elem_name not in elem_to_vars:
                            elem_to_vars[elem_name] = []
                        if var not in elem_to_vars[elem_name]:
                            elem_to_vars[elem_name].append(var)
    return elem_to_vars


if __name__ == "__main__":
    _, _elems_file, _sequences, _output = sys.argv
    main(_elems_file, _sequences, _output)
