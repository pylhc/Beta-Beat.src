import sys
import numpy as np

USAGE_STR = """
Usage:

To check if the file contains bad BPM combinations:
python read_bet_deviations.py check input_file

To try to remove bad BPM combinations from the file:
python read_bet_deviations.py fix input_file output_file
"""

NO_ALLOWED_WORDS = ["DOROS", "ITLK"]

def main(mode, input_file, output_file=None):
    assert mode in ["fix", "check"]

    loaded_file = np.load(input_file)
    horizontal = loaded_file[0]
    new_horizontal = {}
    for key, value in horizontal.iteritems():
        is_good = True
        for word in NO_ALLOWED_WORDS:
            if word in key:
                is_good = False
                if mode == "check":
                    print "Found bad BPM combination:", key
                    return
                break
        if is_good:
            new_horizontal[key] = value

    vertical = loaded_file[1]
    new_vertical = {}
    for key, value in vertical.iteritems():
        is_good = True
        for word in NO_ALLOWED_WORDS:
            if word in key:
                is_good = False
                if mode == "check":
                    print "Found bad BPM combination:", key
                    return
                break
        if is_good:
            new_vertical[key] = value
    if mode == "check":
        print "The file doesn't containt bad BPMs."
    else:
        new_file_content = np.array([new_horizontal, new_vertical])
        np.save(output_file, new_file_content)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print USAGE_STR
