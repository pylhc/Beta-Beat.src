def can_str_be_parsed_to_number(str_to_test):
    """ Checks if the given string is a float. """
    try:
        float(str_to_test)
        return True
    except ValueError:
        return False
