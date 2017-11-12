# -*- coding: utf-8 -*-
"""
Configuration manager for workflow definitions.
This works over Kaptan:https://github.com/emre/kaptan

The global configuration registry is declared in the bottom of this file.
"""
import os.path as op

from   nipype.pipeline.engine import Node, MapNode, JoinNode
from   nipype.interfaces.base import isdefined
from   kaptan import Kaptan


def _load_config(file_path):
    cpt = Kaptan()
    return cpt.import_config(file_path)


def _check_file(file_path):
    fpath = op.abspath(op.expanduser(file_path))

    if not op.isfile(fpath):
        raise IOError("Could not find configuration file {}.".format(fpath))


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

    def update(self, adict):
        for k, v in adict.items():
            self._cpt.upsert(k, v)

    def keys(self):
        return self._cpt.configuration_data.keys()

    def items(self):
        return self._cpt.configuration_data.items()

    def get(self, item, default=None):
        return self._cpt.configuration_data.get(item, default)

    def __getitem__(self, item):
        if item not in self._cpt.configuration_data:
            raise KeyError('Could not find key {} in configuration content.'.format(item))
        return self._cpt.configuration_data[item]

    def __setitem__(self, key, value):
        return self._cpt.upsert(key, value)

    def __repr__(self):
        return '<config.Config> ({})'.format('\n'.join([str(i) for i in self.items()]))


# the global configuration registry
PYPES_CFG = Config()


def node_settings(node_name):
    global PYPES_CFG
    for k, v in PYPES_CFG.items():
        if k.startswith(node_name):
            yield k, v


def update_config(value):
    """ Value can be a configuration file path or a dictionary with
    configuration settings."""
    global PYPES_CFG
    if isinstance(value, str):
        PYPES_CFG.update_from_file(value)
    elif isinstance(value, dict):
        PYPES_CFG.update(value)
    else:
        raise NotImplementedError('Cannot update the configuration with {}.'.format(value))


def _set_node_inputs(node, params, overwrite=False):
    for k, v in params.items():
        try:
            if not isdefined(getattr(node.inputs, k)):
                setattr(node.inputs, k, v)
            else:
                if overwrite:
                    setattr(node.inputs, k, v)
        except AttributeError as ate:
            raise AttributeError('Error in configuration settings: node `{}` '
                                 'has no attribute `{}`.'.format(node, k)) from ate


def _get_params_for(node_name):
    pars = {}
    for k, v in node_settings(node_name):
        nuk = '.'.join(k.split('.')[1:]) if '.' in k else k
        pars[nuk] = v

    return pars


def check_mandatory_inputs(node_names):
    """ Raise an exception if any of the items in the List[str] `node_names` is not
    present in the global configuration settings."""
    for name in node_names:
        if name not in PYPES_CFG:
            raise AttributeError('Could not find a configuration parameter for {}. '
                                 'Please set it in the an input configuration file.'.format(name))


def get_config_setting(param_name, default=''):
    """ Return the value for the entry with name `param_name` in the global configuration."""
    return PYPES_CFG.get(param_name, default=default)


def setup_node(interface, name, settings=None, overwrite=True, **kwargs):
    """ Create a pe.Node from `interface` with a given name.
    Check in the global configuration if there is any value for the node name and will set it.

    Parameters
    ----------
    interface: nipype.interface

    name: str

    settings: dict
        Dictionary with values for the pe.Node inputs.
        These will have higher priority than the ones in the global Configuration.

    overwrite: bool
        If True will overwrite the settings of the node if they are already defined.
        Default: True

    kwargs: keyword arguments
        type: str or None.
            choices: 'map', 'join, or None.
            If 'map' will return a MapNode.
            If 'join' will return a JoinNode.
            If None will return a Node.

        Extra arguments to pass to nipype.Node __init__ function.


    Returns
    -------
    node: nipype.Node
    """
    typ = kwargs.pop('type', None)
    if typ == 'map':
        node_class = MapNode
    elif typ == 'join':
        node_class = JoinNode
    else:
        node_class = Node
    node = node_class(interface=interface, name=name, **kwargs)

    params = _get_params_for(name)
    if settings is not None:
        params.update(settings)

    _set_node_inputs(node, params, overwrite=overwrite)

    return node


# ---------------------------------------------------------------------------
# Helper functions for specific parameters of config
# ---------------------------------------------------------------------------
def check_atlas_file():
    """ Return True and the path to the atlas_file if the configuration settings
    `normalize_atlas` is True and `atlas_file` points to an existing file.
    If `normalize_atlas` is False will return False and an empty string.
    Otherwise will raise a FileNotFoundError.

    Returns
    -------
    do_atlas: bool

    atlas_file: str
        Existing file path to the atlas file

    Raises
    ------
    FileNotFoundError
        If the `normalize_atlas` option is True but the
        `atlas_file` is not an existing file.
    """
    normalize_atlas = get_config_setting('normalize_atlas', default=False)
    if not normalize_atlas:
        return False, ''

    atlas_file = get_config_setting('atlas_file', default='')
    if not op.isfile(atlas_file):
        raise FileNotFoundError('Could not find atlas file in {}. '
                                'Please set `normalize_atlas` to False '
                                'or give an existing atlas image.'.format(atlas_file))
    return True, atlas_file

