def can_str_be_parsed_to_number(str_to_test):
    """ Checks if the given string is a float.
    "str" --> False,
    "23.45" --> True,
    "2E3" --> True,
    "2e-3" --> True
    """
    try:
        float(str_to_test)
        return True
    except ValueError:
        return False
