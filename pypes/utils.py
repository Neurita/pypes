"""
Helper nipype nodes.
"""
import os
import os.path as op

import nipype.interfaces.spm as spm
import nipype.interfaces.fsl as fsl
from   nipype.interfaces.base import traits
from   nipype.interfaces.io import isdefined
from   nipype.interfaces.utility import Function


def spm_tpm_priors_path(spm_dir=None):
    """ Return the path to the TPM.nii file from SPM.

    Parameters
    ----------
    spm_dir: str
        Path to SPM.
        If `None` will try nipype to guess it. #TOBEDONE

    Returns
    -------
    tpm_path: str
        Path to the TPM.nii file.

    Raises
    ------
    FileNotFoundError
        If `template` is `None` and can't find the TPM.nii file from SPM.
    """
    if spm_dir is None:
        spm_dir = spm.Info.version().get('path', None)

    if spm_dir is None:
        spm_dir = op.expanduser('~/Software/matlab_tools/spm12')

    tpm_path = op.join(spm_dir, 'tpm', 'TPM.nii')
    if not op.exists(tpm_path):
        raise FileNotFoundError('Could not find TPM.nii file from SPM.')

    return tpm_path


def get_fsl_extension(default='NIFTI_GZ'):
    """ Return the file extension that FSL will output.

    Parameters
    ----------
    default: str
        Default output type for FSL in FSL format.
        Choices: ANALYZE, NIFTI, NIFTI_PAIR,
                ANALYZE_GZ, NIFTI_GZ, NIFTI_PAIR_GZ

    Returns
    -------
    ext: str
        The file extension

    Note
    ----
    In the case of the pair files, this will only return
    the extension of the header file.

    #TODO: deal with pair file renaming with header modification.
    Check boyle for this.
    """
    outtype = os.environ.get('FSLOUTPUTTYPE', default)

    otype_ext = {'NIFTI':           '.nii',
                 'NIFTI_GZ':        '.nii.gz',
                 'ANALYZE':         '.hdr',
                 'NIFTI_PAIR':      '.hdr',
                 'NIFTI_PAIR_GZ':   '.hdr.gz',
                 'ANALYZE_PAIR_GZ': '.hdr.gz',
                }

    return otype_ext[outtype]


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

