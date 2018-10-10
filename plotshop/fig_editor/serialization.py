from utils.contexts import suppress_exception
from utils.dict_tools import DotDict




def line_to_struct(line):
    """ Converts a line to a dictionary structure of plotly-standard.
    See https://plot.ly/python/reference/#scatter
    """
    ls = DotDict()
    ls.name = line.get_name()
    ls.visible = line.get_visible()
    ls.x = line.get_xdata()
    ls.y = line.get_ydata()
    ls.mode = _get_mode(line)

    ls.line = DotDict()
    ls.line.color = line.get_color()
    ls.line.width = line.get_width()
    ls.line.dash = line.get_linestyle()


def struct_to_line(ax, struct):
    pass



def linestyle_to_dash(ls):
    """ Map linestyle to dash_names """
    return {
        "-": "solid",
        "--": "dash",
        ".": "dot",
        "-.": "dashdot",
    }[ls]


def dash_to_linestyle(dash):
    try:
        return {
            "solid": "-",
            "dash": "--",
            "dot": ".",
            "dashdot": "-.",
            "longdash": "--",
            "longdashdot": "-.",
        }[dash]
    except KeyError:
        return "--"


def _get_mode(line):
    """ Returns 'lines', 'markers' or 'lines+markers'.
    """
    mode = []
    ls = line.get_linestyle()
    ms = line.get_marker()
    if "none" != ls.lower():
        mode.append("lines")
    if "none" != ms.lower():
        mode.append("markers")
    return "+".join(mode)


def _set_mode(line, mode):
    """ Sets lines or markers to None if not in mode. """
    if "lines" not in mode:
        line.set_linestyle("None")
    if "markers" not in mode:
        line.set_marker("None")