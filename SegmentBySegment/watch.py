import os
import time
import sys


def main():
    path_to_watch1 = sys.argv[1]  # This must be the corrections file
    path_to_watch2 = sys.argv[2]  # This must be gplot file
    sbs_command = sys.argv[3].strip("\"")  # Full SbS command to run

    gplot = path_to_watch2.find("gplot")  # string index for gplot in gplot path

    path = path_to_watch2[:gplot]     # Extract path from gplot filename
    label = path_to_watch2[gplot + 6:]  # Extract label from gplot filename

    epsfile = path + "plot_" + label + ".phasetotal.eps"

    os.system("emacs " + path_to_watch1 + " &")
    os.system("gv --watch " + epsfile + " &")

    a = os.stat(path_to_watch1)
    before1 = a.st_mtime
    a = os.stat(path_to_watch2)
    before2 = a.st_mtime

    while True:
        time.sleep(1)
        a = os.stat(path_to_watch1)
        after1 = a.st_mtime
        a = os.stat(path_to_watch2)
        after2 = a.st_mtime
        if after1 > before1:  # MAD file changed, run everything
            print "Changed: ", path_to_watch1
            os.system("python " + sbs_command + "; gnuplot " + path_to_watch2)
        if  after2 > before2:        # Only gplot file changed
            print "Changed: ", path_to_watch2
            os.system("gnuplot " + path_to_watch2)

        before1 = after1
        before2 = after2


if __name__ == "__main__":
    main()
