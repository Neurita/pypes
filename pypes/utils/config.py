# -*- coding: utf-8 -*-
"""
Configuration manager for workflow definitions.
This works over Kaptan:https://github.com/emre/kaptan

The global configuration registry is declared in the bottom of this file.
"""
import os.path as op

from   kaptan import Kaptan, HANDLER_EXT

from pypes.utils import get_extension


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


# the global configuration registry
PYPES_CFG = Config()


def node_settings(node_name):
    #TODO: items from PYPES_CFG that startwith node_name
    global PYPES_CFG
    for k, v in PYPES_CFG.keys():
        if k.startswith(node_name):
            yield k, v


def update_config(file_path):
    global PYPES_CFG
    PYPES_CFG.update_from_file(file_path)


def _set_node_inputs(node, params):
    for k, v in params.items():
        setattr(node.inputs, k, v)


def _get_params_for(node_name):
    pars = {}
    for k, v in node_settings(node_name):
        nuk = '.'.join(k.split('.')[1:]) if '.' in k else k
        pars[nuk] = v

    return pars


def setup_node(node_element, name, inputs=None):
    """ Create a pe.Node from `node_element` with a given name.
    Check in the global configuration if there is any value for the node name and will set it.

    Parameters
    ----------
    node_element: nipype.interface

    name: str

    inputs: dict
        Dictionary with values for the pe.Node inputs.
        These will have higher priority than the ones in the global Configuration.

    Returns
    -------
    node: nipype.Node
    """
    node = pe.Node(node_element, name=name)

    params = _get_params_for(name)
    if inputs is not None:
        params.update(inputs)

    _set_node_inputs(node, params)

    return node