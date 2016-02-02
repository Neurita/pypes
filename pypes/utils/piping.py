# -*- coding: utf-8 -*-
"""
Helper functions for joining, merging, managing the workflow nodes.
"""
from   nipype import Function
from   nipype.interfaces.traits_extension import isdefined
from   nipype.interfaces.base import traits
from   nipype.interfaces.io import SelectFiles, DataGrabber, DataSink
import nipype.interfaces.fsl as fsl


def get_input_node(wf):
    """ Return the file input node in `wf`

    Parameters
    ----------
    wf: nipype Workflow

    Returns
    -------
    nodes: nipype Node
        The file input nodes in `wf`.
    """
    input_types = (SelectFiles, DataGrabber)
    return find_wf_node(wf, input_types)


def get_datasink(wf):
    """ Return the DataSink node in `wf`

    Parameters
    ----------
    wf: nipype Workflow

    Returns
    -------
    nodes: nipype Node
        The DataSink in `wf`.
    """
    return find_wf_node(wf, DataSink)


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


def find_wf_node(wf, iface_type):
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
    for name in wf.list_node_names():
        if isinstance(getattr(wf.get_node(name), 'interface', None), iface_type):
            return wf.get_node(name)
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
    fi: nipype.interfaces.utility.Function
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
    fi: nipype.interfaces.utility.Function
    """
    arg_names = ['arg{}'.format(n) for n in range(1, n_args+1)]

    func = 'def func({0}): import os; return os.path.join({0})'.format(', '.join(arg_names))
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