from __future__ import print_function
import sys
import json
import pandas as pd
from collections import OrderedDict

sys.path.append("/afs/cern.ch/work/j/jcoellod/public/Beta-Beat.src")

from Utilities import tfs_pandas


def main(all_listss, sequence, model, output):
    var_to_classes = get_var_to_classes(all_listss.split(","))
    elem_to_vars = get_elem_to_vars(var_to_classes.keys(), sequence)
    model_data = tfs_pandas.read_tfs(model).set_index("NAME")

    s = []
    elements = []
    vars = []
    classes = []
    for elem in elem_to_vars:
        try:
            model_data.loc[elem, "S"]
        except KeyError:
            continue
        elements.append(elem)
        vars.append(",".join(elem_to_vars[elem]))
        this_classes = []
        for var in elem_to_vars[elem]:
            for cls in var_to_classes[var]:
                if cls not in this_classes:
                    this_classes.append(cls)
        classes.append(",".join(this_classes))
        s.append(model_data.loc[elem, "S"])

    all_data = OrderedDict((
        ("S", s),
        ("ELEMENT", elements),
        ("CLASSES", classes),
        ("VARS", vars),
    ))
    data_frame = pd.DataFrame(
        data=all_data,
        columns=all_data.keys(),
    ).sort_values(by="S")
    tfs_pandas.write_tfs(data_frame, {}, output)


def get_var_to_classes(all_listss):
    var_to_classes = OrderedDict({})
    for all_lists in all_listss:
        with open(all_lists, "r") as all_list_reader:
            all_list_data = json.load(all_list_reader)
            for var_class in all_list_data:
                for var in all_list_data[var_class]:
                    if type(var) == dict:
                        continue
                    if not var.startswith("k"):
                        continue
                    if var not in var_to_classes:
                        var_to_classes[var] = []
                    var_to_classes[var].append(var_class)
    return var_to_classes


def get_elem_to_vars(vars, sequence):
    elem_to_vars = OrderedDict({})
    with open(sequence, "r") as seq_data:
        for line in seq_data:
            for var in vars:
                if var in line and ":=" in line:
                    parts = line.split()
                    elem_name = parts[0].replace(",", "")
                    if elem_name not in elem_to_vars:
                        elem_to_vars[elem_name] = []
                    if var not in elem_to_vars[elem_name]:
                        elem_to_vars[elem_name].append(var)
    return elem_to_vars


if __name__ == "__main__":
    _, _all_listss, _sequence, _model, _output = sys.argv
    main(_all_listss, _sequence, _model, _output)
