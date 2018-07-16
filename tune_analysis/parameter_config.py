# Global 'Parameters' for easy editing #########################################
def get_planes():
    return "XY"


def get_time_col():
    return "TIME"


def get_tstart_head():
    return "START_TIME"


def get_tend_head():
    return "END_TIME"


def get_bbq_col(plane):
    return 'BBQ_{:s}'.format(plane)


def get_mav_col(plane):
    return "{:s}_MAV".format(get_bbq_col(plane))


def get_used_in_mav_col(plane):
    return "{:s}_IN_MAV".format(get_bbq_col(plane))


def get_natq_col(plane):
    return 'NATQ{:s}'.format(plane)


def get_corrected_col(plane):
    return "{:s}_CORRECTED".format(get_natq_col(plane))


def get_timber_key(plane, beam):
    return 'lhc.bofsu:eigen_freq_{:d}_b{:d}'.format({"X": 1, "Y": 2}[plane], beam)
