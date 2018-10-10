import os
import matplotlib
from PyQt5 import QtGui

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
IMAGE_DIR = os.path.join(THIS_DIR, 'images')


# Images #####################################################################


def get_image_dir():
    return IMAGE_DIR


def get_image_path(filename):
    return os.path.join(IMAGE_DIR, filename)


def get_icon(name):
    if os.path.splitext(name)[1] == "":
        name += ".png"
    img_file = get_image_path(name)
    if not os.path.isfile(img_file):
        img_file = os.path.join(matplotlib.rcParams['datapath'], 'images', name)
        if not os.path.isfile(img_file):
            raise IOError("Image '{}' not found.".format(img_file))
    return QtGui.QIcon(img_file)

