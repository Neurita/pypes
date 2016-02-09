# -*- coding: utf-8 -*-
"""
Configuration manager for workflow definitions.
This works over Kaptan:https://github.com/emre/kaptan
"""
import os
import os.path as op

from   kaptan import Kaptan, HANDLER_EXT

from .utils import remove_ext, get_extension


class Borg:
    __shared_state = {}

    def __init__(self):
        self.__dict__ = self.__shared_state
    # and whatever else you want in your class -- that's all!


class Config(Kaptan):
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

    def __init__(self, handler):
        self.__dict__ = self.__shared_state

        super(Config, self).__init__(handler=handler)

    @staticmethod
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

    def from_file(self, file_path):
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
        fpath = op.abspath(op.expanduser(file_path))

        if not op.isfile(fpath):
            raise IOError("Could not find file {}.".format(fpath))

        fname   = op.basename(file_path)
        fext    = get_extension(fname)[1:]
        fpyname = remove_ext(fname)
        handler = self._handler_for_ext(fext)

        cpt = Kaptan(handler=handler)

        cwd = op.abspath(op.curdir)
        is_neighbour = op.samefile(op.dirname(fpath), cwd)
        if is_neighbour:
            _ = cpt.import_config(fpyname)

        else:
            os.chdir(op.dirname(fpath))
            _ = cpt.import_config(fpyname)
            os.chdir(cwd)

        self.__dict__.update(cpt.__dict__)

        return self

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

    def __getitem__(self, item):
        if item not in self.configuration_data:
            raise KeyError('Could not find key {} in configuration content.'.format(item))

        return self.configuration_data[item]

    def __setitem__(self, key, value):
        return self.setitem(key, value)
