# -*- coding: utf-8 -*-
"""
Helper functions for joining, merging, managing the workflow nodes.
"""

from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink, DataGrabber
from nipype.interfaces.base import traits, isdefined
import nipype.interfaces.fsl as fsl

from ..crumb  import DataCrumb


def get_trait_value(traitspec, value_name, default=None):
    """ Return the attribute `value_name` from traitspec if it is defined.
    If not will return the value of `default`.
    Parameters
    ----------
    traitspec: TraitedSpec

    value_name: str
        Name of the `traitspect` attribute.

    default: any
        A default value in case the attribute does not exist or is not defined.

    Returns
    -------
    trait_value: any
    """
    val = getattr(traitspec, value_name, default)
    return default if not isdefined(val) else val


def selectindex(files, idx, flatten=True):
    """ Select the items in the list `files` in the indexes given by the list of
    integers `idx`."""
    import numpy as np
    from nipype.utils.filemanip import filename_to_list, list_to_filename
    from neuro_pypes._utils import flatten_list

    if flatten:
        files = flatten_list(files)

    if isinstance(idx, list):
        return list_to_filename(np.array(filename_to_list(files))[idx].tolist())
    else:
        return files[idx]


def get_node(wf, node_types, name=''):
    """ Return the first node found in `wf` of any of the types in `node_types`.

    Parameters
    ----------
    wf: nipype Workflow

    node_types: sequence of node types

    Returns
    -------
    nodes: nipype Node
        The DataSink in `wf`.
    """
    node = None
    for ionode in node_types:
        node = find_wf_node(wf, ionode, name=name)
        if node is not None:
            return node

    if node is None:
        raise KeyError('Could not find a node of type {} in the worflow.'.format(node_types))


def get_datasink(wf, name=''):
    """ Return the DataSink node in `wf`

    Parameters
    ----------
    wf: nipype Workflow

    Returns
    -------
    nodes: nipype Node
        The DataSink in `wf`.
    """
    return get_node(wf, (DataSink, ), name=name)


def get_input_node(wf, name=''):
    """ Return the first node of type: (DataCrumb, SelectFiles, DataGrabber) in wf."""
    return get_node(wf, (DataCrumb, SelectFiles, DataGrabber), name=name)


def get_interface_node(wf, name):
    """ Return the first node found of type IdentityInterface in `wf` with name `name`."""
    return get_node(wf, (IdentityInterface, ), name=name)


def wf_nodes(wf, iface_type):
    """ Return the nodes found in the list of node names in `wf` that
    has an interface of type `iface_type`.

    Parameters
    ----------
    wf: nipype Workflow

    iface_type: nipype class

    Returns
    -------
    nodes: list of nipype Node
    """
    nodes = []
    for name in wf.list_node_names():
        if isinstance(getattr(wf.get_node(name), 'interface', None), iface_type):
            nodes.append(wf.get_node(name))
    return nodes


def find_wf_node(wf, iface_type, name=''):
    """ Return the first node found in the list of node names in `wf` that
    has an interface of type `iface_type`.

    Parameters
    ----------
    wf: nipype Workflow

    iface_type: nipype class

    Returns
    -------
    node: nipype Node
    """
    for node_nom in wf.list_node_names():
        if isinstance(getattr(wf.get_node(node_nom), 'interface', None), iface_type):
            if not name:
                return wf.get_node(node_nom)
            elif name in node_nom.split('.'):
                return wf.get_node(node_nom)
    return None


def extend_trait_list(trait_list, extension):
    """ Extend or initialize `trait_list` with `extension`.

    Parameters
    ----------
    trait_list: traits.List

    extension: list or tuple

    Returns
    -------
    trait_list: traits.List
        The extended list.
    """
    if isdefined(trait_list):
        if not isinstance(extension, (list, tuple)):
            raise ValueError('Expected `extension` to be list or tuple, got {}.'.format(type(extension)))

        trait_list.extend(extension)
    else:
        trait_list = extension

    return trait_list


def sum2args():
    """ Return a nipype function that sums up two args: `arg1` and `arg2` and
    leaves the result in `out`.

    Returns
    -------
    fi: nipype.interfaces.Function
    """
    func = 'def func(arg1, arg2): return arg1 + arg2'
    fi = Function(input_names=['arg1', 'arg2'], output_names=['out'])
    fi.inputs.function_str = func
    return fi


def joinstrings(n_args=2):
    """ Return a nipype function that joins up to `n_args` parts. and
    leaves the result in `out`.
    Parameters
    ----------
    n_args: int
        Number of `argX` parameters the function might have.
        Example: If `n_args` == 2, then the node will have `arg1` and `arg2` nodes.

    Returns
    -------
    fi: nipype.interfaces.Function
    """
    arg_names = ['arg{}'.format(n) for n in range(1, n_args+1)]

    func = '''def func({0}): import os; return os.path.join({0})'''.format(', '.join(arg_names))
    fi = Function(input_names=arg_names, output_names=['out'])
    fi.inputs.function_str = func
    return fi


def fsl_merge(in_files=traits.Undefined, dimension='t'):
    """ Merge the NifTI files in `in_files` in the given `dimension`.
    This uses `fslmerge`.

    Parameters
    ----------
    in_files: list of str.
        Paths to the files to merge.

    dimension: str
        Character indicating the merging dimension.
        Choices: 't', 'x', 'y', 'z'

    Returns
    -------
    merger: fsl.Merge
    """
    merger = fsl.Merge()
    merger.inputs.dimension = dimension
    merger.inputs.output_type = "NIFTI_GZ"
    merger.inputs.in_files = in_files

    return merger


def get_input_file_name(input_node, fname_key):
    """ Return the name of the file given by the node key `fname_key` in `input_node`.

    Parameters
    ----------
    input_node: nipype Node
        a node with a file input interface (SelectFiles, DataCrumb), for now.

    fname_key: str
        The key that is used to access the file path using `input_node`.
        Example: 'anat'

    Returns
    -------
    filename: str
        The path template of the input file `fname_key` path from the `input_node`.
    """
    if isinstance(input_node.interface, SelectFiles):
        try:
            fname = input_node.interface._templates[fname_key]
        except AttributeError:
            raise
        except KeyError:
            raise KeyError('Could not find file key {} in node {}({}).'.format(fname_key, input_node.name,
                                                                           type(input_node.interface)))
        else:
            return fname

    if isinstance(input_node.interface, DataCrumb):
        try:
            crumb_args = input_node.interface._templates[fname_key]
            incrumb    = input_node.interface._crumb.replace(**dict(crumb_args))
        except AttributeError:
            raise
        except KeyError:
            raise KeyError('Could not find file key {} in node {}({}).'.format(fname_key, input_node.name,
                                                                           type(input_node.interface)))
        else:
            return incrumb.path

    else:
        raise NotImplementedError('`get_input_file_name` has not been implemented for nodes'
                                  ' of type {}.'.format(type(input_node.interface)))