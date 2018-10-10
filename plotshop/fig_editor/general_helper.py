from contextlib import contextmanager

@contextmanager
def suppress_exception(exception):
    """ Catch exception and ignore it. """
    try:
        yield
    except exception:
        pass
