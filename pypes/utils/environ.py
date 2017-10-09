# -*- coding: utf-8 -*-
"""
Helper functions to check for value configurations in the system and the external tools.
"""
import os
from os import path as op

from nipype.interfaces import spm

from ..config import get_config_setting


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
    NotADirectoryError
        If the SPM path is not found.

    FileNotFoundError
        If `template` is `None` and can't find the TPM.nii file from SPM.
    """
    matlab_cmd = os.environ.get('SPMMCRCMD')
    spm_version = spm.Info.version(matlab_cmd=matlab_cmd,
                                   use_mcr='FORCE_SPMMCR' in os.environ)

    if spm_version is None:
        raise RuntimeError("Nipype could not find a valid Matlab or SPM configuration.")

    if spm_dir is None:
        spm_dir = spm_version.get('path', None)

    if spm_dir is None:
        spm_dir = op.expanduser(get_config_setting('spm_dir'))

    if not spm_dir:
        raise NotADirectoryError('Could not find SPM path.')

    if not op.exists(spm_dir):
        raise NotADirectoryError('The specified SPM path ({}) does not exist.'.format(spm_dir))

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
