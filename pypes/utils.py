"""
Helper nipype nodes.
"""
import os.path as op

import nipype.interfaces.spm as spm
import nipype.interfaces.fsl as fsl
from   nipype.interfaces.base import traits


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
    # TODO use spm.Info.version() to pick the path to SPM
    if spm_dir is None:
        spm_dir = '/home/alexandre/Software/matlab_toolboxes/spm12'

    tpm_path = op.join(spm_dir, 'tpm', 'TPM.nii')
    if not op.exists(tpm_path):
        raise FileNotFoundError('Could not find TPM.nii file from SPM.')

    return tpm_path


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