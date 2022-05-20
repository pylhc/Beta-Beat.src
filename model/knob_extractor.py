import argparse
import math
from datetime import datetime, timedelta
import re
import os

KNOBS_TXT_MDLDIR = "acc-models-lhc/operation/knobs.txt"
KNOBS_TXT_FALLBACK =  "/afs/cern.ch/eng/acc-models/lhc/current/operation/knobs.txt"

KNOB_NAMES = {
    "sep": [
        "LHCBEAM:IP1-SEP-H-MM",
        "LHCBEAM:IP1-SEP-V-MM",
        "LHCBEAM:IP5-SEP-H-MM",
        "LHCBEAM:IP5-SEP-V-MM",
    ],
    "xing": [
        "LHCBEAM:IP1-XING-V-MURAD",
        "LHCBEAM:IP1-XING-H-MURAD",
        "LHCBEAM:IP5-XING-V-MURAD",
        "LHCBEAM:IP5-XING-H-MURAD",
    ],
    "chroma": [
        "LHCBEAM1:QPH",
        "LHCBEAM1:QPV",
        "LHCBEAM2:QPH",
        "LHCBEAM2:QPV",
    ],
    "ip_offset": [
        "LHCBEAM:IP1-OFFSET-V-MM",
        "LHCBEAM:IP2-OFFSET-V-MM",
        "LHCBEAM:IP5-OFFSET-H-MM",
        "LHCBEAM:IP8-OFFSET-H-MM",
    ],
    "disp": [
        "LHCBEAM:IP1-SDISP-CORR-SEP",
        "LHCBEAM:IP1-SDISP-CORR-XING",
        "LHCBEAM:IP5-SDISP-CORR-SEP",
        "LHCBEAM:IP5-SDISP-CORR-XING",
    ],
    "mo": [
        "LHCBEAM1:LANDAU_DAMPING"
        "LHCBEAM2:LANDAU_DAMPING"
    ]
}


def main():
    parser = argparse.ArgumentParser("Knob extraction tool.")

    parser.add_argument("knobs", type=str, nargs='*',
                        help="a list of knob categories to extract",
                        choices=[key for key in KNOB_NAMES.keys()]
                        )
    parser.add_argument("--extract", type=str, nargs='+',
                        help=("defines the time at which to extract the knob settings. "
                              "'now': extracts the current knob setting. "
                              "<time>: extracts the knob setting for a given time. "
                              "<time> <timedelta>: extracts the knob setting for a given time, "
                              "<timedelta> is a timedelta specification (7m = 7 minutes, 1d = 1day, 7m30s = 7 min 30 secs) "
                              "a prefix '_' specifies a negative timedelta"
                              )
                        )
    parser.add_argument("--state", action='store_true',
                        help="prints the state of the statetracker")
    parser.add_argument("--output", type=str, default='knobs.madx',
                        help="specify user-defined output path. This should probably be `model_dir/knobs.madx`"
                        )

    args = parser.parse_args()

    if args.extract is not None:
        end = args.extract[1] if len(args.extract) > 1 else None
        _extract(args.knobs, args.extract[0], end, args.output)

    if args.state:
        print("---- STATE ------------------------------------")
        print(ldb.get("LhcStateTracker:State", t1))
        print(ldb.get("LhcStateTracker/State", t1))


def _extract(knobs, start, end = None, output="./knobs.madx"):
    import pytimber
    print("---- EXTRACTING KNOBS -------------------------")

    t1 = _time_from_str(start)

    if end is not None:
        t1 = _add_delta(t1, end)

    print(f"extracting knobs for {t1}")

    knobdict = _get_knobs_dict()

    ldb = pytimber.LoggingDB(source="nxcals")
    print("---- KNOBS ------------------------------------")
    with open(output, "w") as outfile:
        for knob in knobs:
            outfile.write(f"!! --- {knob:10} --------------------\n")
            for knobname in KNOB_NAMES[knob]:
                print(f"looking for {knobname}")
                knobkey = f"LhcStateTracker:{knobname}:target"
                knobvalue = ldb.get(knobkey, t1)
                print(knobvalue)
                if not knobkey in knobvalue:
                    outfile.write(f"! no value for {knobname}\n")
                    continue
                (timestamps, values) = knobvalue[knobkey]
                value = values[-1]
                (madxname, scaling) = knobdict[knobname.replace(":", "/")]
                print(f"{madxname}: {value*scaling}")
                if not math.isnan(value):
                    outfile.write(f"{madxname} := {value*scaling};\n")
            outfile.write("\n")


def _time_from_str(pattern):
    if pattern == "now":
        return datetime.now()
    try:
        return datetime.fromisoformat(pattern)
    except:
        pass
    try:
        return datetime.fromtimestamp(int(pattern))
    except:
        pass

    print(f"couldn't read datetime '{pattern}'")
    return None


def _add_delta(t1, pattern):
    is_negative = pattern.startswith('_')
    sign = 1

    if is_negative:
        sign = -1

    deltare = re.compile("(\\d+)(\\w)")
    all_deltas = deltare.findall(pattern)

    for delta in all_deltas:
        unit = delta[1]
        value = sign*int(delta[0])
        if unit == 's':
            t1 = t1 + timedelta(seconds=value)
        if unit == 'm':
            t1 = t1 + timedelta(minutes=value)
        if unit == 'h':
            t1 = t1 + timedelta(hours=value)
        if unit == 'd':
            t1 = t1 + timedelta(days=value)
        if unit == 'w':
            t1 = t1 + timedelta(days=7*value)
        if unit == 'M':
            t1 = t1 + timedelta(months=value)

    return t1


def _get_knobs_dict(user_defined = None):
    # if all fails, fall back to lhc acc-models
    if user_defined is not None:
        filename = user_defined
        print(f"take user defined knobs.txt: '{filename}")
    elif os.path.isfile(KNOBS_TXT_MDLDIR):
        filename = KNOBS_TXT_MDLDIR
        print(f"take model folder's knobs.txt: '{filename}")
    else:
        # if all fails, fall back to lhc acc-models
        filename = KNOBS_TXT_FALLBACK
        print(f"take fallback knobs.txt: '{filename}'")

    knobdict = {}
    with open(filename) as knobtxt:
        next(knobtxt)
        for line in knobtxt:
            knob = line.split(',')
            if len(knob) == 4:
                knobdict[knob[1].strip()] = (knob[0].strip(), float(knob[2]))

    return knobdict


if __name__ == "__main__":
    main()
