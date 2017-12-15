"""
    Tools to handle errors easier
    (and to not to have to look it up on StackOverflow all the time)
"""


def replace_error_message(err, message):
    """ Replace the message of error 'err' with 'message' """
    if not err.args or len(err.args) == 1:
        err.args = (message,)
    else:
        err.args = (message,) + err.args[1:]
    return err


def append_error_message(err, message):
    """ Append 'message' to the message of error 'err' """
    if not err.args:
        err.args = (message,)
    elif len(err.args) == 1:
        err.args = (err.args[0] + message,)
    else:
        err.args = (err.args[0] + message,) + err.args[1:]
    return err
