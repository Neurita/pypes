# -*- coding: utf-8 -*-
"""
Configuration manager for workflow definitions.
This works over Kaptan:https://github.com/emre/kaptan
"""
import os.path as op

from   kaptan import Kaptan, HANDLER_EXT

from .utils import get_extension


def _handler_for_ext(ext):
    """

    Parameters
    ----------
    ext: str
        The file extension without the dot.

    Returns
    -------
    handler: a Kaptan.BaseHandler
    """
    if ext not in HANDLER_EXT:
        raise ValueError("Expected a file with an extension "
                         "in {}, got {}.".format(HANDLER_EXT.keys(), ext))

    return HANDLER_EXT[ext]


def _handler_for_file(file_path):
    fname = op.basename(file_path)
    fext  = get_extension(fname)[1:]
    return _handler_for_ext(fext)


def _load_config(file_path):
    cpt = Kaptan()
    return cpt.import_config(file_path)


def _check_file(file_path):
    fpath = op.abspath(op.expanduser(file_path))

    if not op.isfile(fpath):
        raise IOError("Could not find file {}.".format(fpath))


def _update_kaptan(kptn, new_values):
    for k, v in new_values.items():
        kptn.upsert(k, v)


class Borg:
    __shared_state = {}

    def __init__(self):
        self.__dict__ = self.__shared_state
    # and whatever else you want in your class -- that's all!


class Config(object):
    """ This class has a shared state like a Borg.

    This is a Kaptan class that infers the file handler
    from the file extension using the `from_file` function.

    Parameters
    ----------
    handler: str or kaptan.BaseHandler
        This parameter is to keep compatibility with Kaptan's instantiation protocol.
        See its documentation on how to use it: http://emre.github.io/kaptan/
        I personally prefer to use the `from_file` function.
    """
    __shared_state = {}

    def __init__(self, handler=None):
        self.__dict__ = self.__shared_state
        self._cpt = Kaptan(handler)

    @classmethod
    def from_file(cls, file_path):
        """ Returns a Config instance with the data from
        `file_path`.

        Parameters
        ----------
        file_path: str
            Path to a configuration file.
            Its extension can be any from config.HANDLER_EXT, i.e.,
            {'conf': 'ini',
             'ini': 'ini',
             'json': 'json',
             'py': 'file',
             'yaml': 'yaml',
             'yml': 'yaml'}

        Returns
        -------
        cfg: Config
        """
        _check_file(file_path)
        cfg = Config()
        cfg._cpt = _load_config(file_path)
        return cfg

    def update_from_file(self, file_path):
        """ Updates the config parameters with the data from
        `file_path`.

        Parameters
        ----------
        file_path: str
            Path to a configuration file.
            Its extension can be any from config.HANDLER_EXT, i.e.,
            {'conf': 'ini',
             'ini': 'ini',
             'json': 'json',
             'py': 'file',
             'yaml': 'yaml',
             'yml': 'yaml'}

        Returns
        -------
        cfg: Config
        """
        _check_file(file_path)
        cpt = _load_config(file_path)

        params = cpt.configuration_data
        _update_kaptan(self._cpt, params)

    def check_file(self, item):
        """ This is a __getitem__ operator for file path values, if the file does not exist an IOError is raised."""
        try:
            fpath = self.__getitem__(item)
        except KeyError:
            raise
        else:
            if not op.exists(fpath):
                raise IOError('Could not find file for key {} in the {}.'.format(item, fpath))
            return fpath

    def keys(self):
        return self._cpt.configuration_data.keys()

    def items(self):
        return self._cpt.configuration_data.items()

    def __getitem__(self, item):
        if item not in self._cpt.configuration_data:
            raise KeyError('Could not find key {} in configuration content.'.format(item))

        return self._cpt.configuration_data[item]

    def __setitem__(self, key, value):
        return self._cpt.upsert(key, value)

    def __repr__(self):
        return '<config.Config> ({})'.format(list(self.items()))