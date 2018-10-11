import pytz


# Global 'Parameters' for easy editing #########################################
def get_planes():
    """ Names for the planes."""
    return "XY"


def get_experiment_timezone():
    """ Get time zone for measurement data. """
    return pytz.timezone("Europe/Zurich")


def get_time_col():
    """ Label for the TIME column."""
    return "TIME"


def get_tstart_head():
    """ Label for fill start time from header. """
    return "START_TIME"


def get_tend_head():
    """ Label for fill end time from header. """
    return "END_TIME"


def get_bbq_col(plane):
    """ Label for the BBQ column """
    return 'BBQ_{:s}'.format(plane.upper())


def get_mav_col(plane):
    """ Label for the moving average BBQ column. """
    return "{:s}_MAV".format(get_bbq_col(plane))


def get_used_in_mav_col(plane):
    """ Label for the column showing if BBQ value was used in moving average. """
    return "{:s}_IN_MAV".format(get_bbq_col(plane))


def get_natq_col(plane):
    """ Label for the natural tune column. """
    return 'NATQ{:s}'.format(plane.upper())


def get_natq_corr_col(plane):
    """ Label for the corrected natural tune column. """
    return "{:s}_CORRECTED".format(get_natq_col(plane))


def get_natq_err_col(plane):
    """ Label for the natural tune error column. """
    return "{:s}RMS".format(get_natq_col(plane))


def get_action_col(plane):
    """ Label for the action column. """
    return "2J{:s}RES".format(plane.upper())


def get_action_err_col(plane):
    """ Label for the action error column. """
    return get_action_col("{:s}STD".format(plane))


def get_paired_columns(plane, other_plane):
    """ Four-Tupel of action/tune pair columns + errors. """
    return {"x": get_action_col(plane), "x_err": get_action_err_col(plane),
            "y": get_natq_corr_col(other_plane), "y_err": get_natq_err_col(other_plane)}


def get_paired_lables(plane, other_plane):
    """ Labels for the action/tune plots. """
    return (r'$2J_{:s} \quad [\mu m]$'.format(plane.lower()),
            r'$\Delta Q_{:s}$'.format(other_plane.lower()))


def get_timber_bbq_key(plane, beam):
    """ Key to extract bbq from timber. """
    return 'lhc.bofsu:eigen_freq_{:d}_b{:d}'.format({"X": 1, "Y": 2}[plane], beam)


# Script Mode #################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
