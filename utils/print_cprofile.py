
import sys
import pstats


def main():
    """
    Little snipet to print the output of cProfile in
    a human redable way.
    Usage: python print_cprofile.py input_file.dat >> ouput_file.txt
    """
    _, file_path = sys.argv
    p = pstats.Stats(file_path)
    p.sort_stats('cumtime')
    p.print_stats()


if __name__ == "__main__":
    main()
