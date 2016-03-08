# -*- coding: utf-8 -*-
"""
Helper functions for joining, merging, managing the workflow nodes.
"""
from nipype import Function, Node, SelectFiles, DataSink, DataGrabber
import nipype.interfaces.fsl as fsl
from nipype.interfaces.base import (traits, isdefined)

from ..crumb  import DataCrumb


def iterable_record_node(records, node_name):
    """ A node to iterate over each set of parameters given in records, e.g.
            [[('subjid', '001'), ('diagnosis', 'ad')],
              ('subjid', '002'), ('diagnosis', 'hc')],
              ...
            ]
    Parameters
    ----------
    records: list of list of 2-tuple of str
        The set of

    node_name: str
        A name for the node

    Returns
    -------
    node: nipype.Node

    Raises
    ------
    IndexError
        If the keys in the records are not the same.
    """
    # check if all of them have the same set of keys
    one_fields = set(dict(records[0]).keys())
    for idx, rec in enumerate(records):
        if one_fields != set(dict(rec).keys()):
            raise ValueError('Expected that all `records` were the same, '
                             'got sets of keys {} in index {}, given '
                             'the first set {}.'.format(set(dict(rec).keys()),
                                                        idx, one_fields))

    # defined the getter function
    def get_record(index, items, fields):
        idict = dict(items[index])
        if len(fields) == 1:
            return idict[fields[0]]
        else:
            return tuple([idict[fname] for fname in fields])

    # define the get_record Function interface
    fiface = Function(input_names=['index', 'items', 'fields'],
                      output_names=one_fields,
                      function=get_record)

    # define the node
    node = Node(fiface, name=node_name)
    node.inputs.items = records
    node.inputs.fields = list(one_fields)
    node.iterables = [('index', list(range(len(records)))),]

    return node


def get_node(wf, node_types):
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
        node = find_wf_node(wf, ionode)
        if node is not None:
            return node

    if node is None:
        raise KeyError('Could not find a node of type {} in the worflow.'.format(node_types))


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
    return get_node(wf, (DataSink, ))


def get_input_node(wf):
    """ Return the first node of type: (DataCrumb, SelectFiles, DataGrabber) in wf."""
    return get_node(wf, (DataCrumb, SelectFiles, DataGrabber))


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