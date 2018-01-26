import sys
import os

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

from Utilities import tfs_pandas


class _MetaTfsCollection(type):
    def __new__(mcs, cls_name, bases, dct):
        new_dict = dict(dct)
        for name in dct:
            value = dct[name]
            try:
                args = value.args
                kwargs = value.kwargs
            except AttributeError:
                continue
            new_props = _define_property(dct, args, kwargs)
            try:
                prop_x, prop_y = new_props
                new_dict.pop(name)
                new_dict[name + "_x"] = prop_x
                new_dict[name + "_y"] = prop_y
            except TypeError:
                new_dict[name] = new_props
        return super(_MetaTfsCollection, mcs).__new__(mcs, cls_name, bases, new_dict)


class TfsCollection(object):
    """
    TODO
    """
    __metaclass__ = _MetaTfsCollection

    def __init__(self, directory):
        self.directory = directory
        self._buffer = {}

    def get_filename(*args, **kwargs):
        """
        TODO
        """
        raise NotImplementedError(
            "This is an abstract method, it should be implemented in subclasses."
        )

    def clear(self):
        """
        TODO
        """
        self._buffer = {}

    def _load_tfs(self, filename):
        try:
            return self._buffer[filename]
        except KeyError:
            self._buffer[filename] =\
                tfs_pandas.read_tfs(os.path.join(self.directory, filename))
            return self._buffer[filename]

    def _write_tfs(self, filename, data_frame):
        tfs_pandas.write_tfs(os.path.join(self.directory, filename), data_frame)
        self._buffer[filename] = data_frame


class Tfs(object):
    """
    TODO
    """
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs


# Private methods to define the properties ##################################

def _define_property(dct, args, kwargs):
    if "two_planes" not in kwargs:
        prop_x, prop_y = _define_property_two_planes(dct, args, kwargs)
        return prop_x, prop_y
    elif kwargs["two_planes"]:
        kwargs.pop("two_planes")
        return _define_property_two_planes(dct, args, kwargs)
    else:
        kwargs.pop("two_planes")
        def getter_funct(self):
            return _getter(self, *args, **kwargs)
        def setter_funct(self, tfs_data):
            return _setter(self, tfs_data, *args, **kwargs)
        return property(fget=getter_funct, fset=setter_funct)


def _define_property_two_planes(dct, args, kwargs):
    x_kwargs = dict(kwargs)
    y_kwargs = dict(kwargs)
    x_kwargs["plane"] = "x"
    y_kwargs["plane"] = "y"
    def x_getter_funct(self):
        return _getter(self, *args, **x_kwargs)
    def x_setter_funct(self, tfs_data):
        return _setter(self, tfs_data, *args, **x_kwargs)
    def y_getter_funct(self):
        return _getter(self, *args, **y_kwargs)
    def y_setter_funct(self, tfs_data):
        return _setter(self, tfs_data, *args, **y_kwargs)
    property_x = property(fget=x_getter_funct, fset=x_setter_funct)
    property_y = property(fget=y_getter_funct, fset=y_setter_funct)
    return property_x, property_y


def _getter(self, *args, **kwargs):
    filename = self.get_filename(*args, **kwargs)
    return self._load_tfs(filename)


def _setter(self, data_frame, *args, **kwargs):
    # TODO: Get filename for writing.
    filename = self.get_filename(*args, **kwargs)
    self._write_tfs(filename, data_frame)
