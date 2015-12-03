"""
Helper functions for joining, merging, managing the workflow nodes.
"""
from   nipype import Function
from   nipype.interfaces.traits_extension import isdefined
from   nipype.interfaces.base import traits
import nipype.interfaces.fsl as fsl


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


def joinpaths():
    """ Return a nipype function that joins up two path parts: `arg1` and `arg2` and
    leaves the result in `out`.

    Returns
    -------
    fi: nipype.interfaces.utility.Function
    """
    func = 'def func(arg1, arg2): import os; return os.path.join(arg1, arg2)'
    fi = Function(input_names=['arg1', 'arg2'], output_names=['out'])
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