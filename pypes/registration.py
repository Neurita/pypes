"""
Nipype registration nodes
"""
import nipype.interfaces.spm as spm
from   nipype.interfaces.base import traits

from .utils import spm_tpm_priors_path


def spm_apply_deformations(in_imgs=traits.Undefined, trans_field=traits.Undefined):
    """Return a Normalize12 interface object.
    SPM12’s new Normalise routine for warping an image to a template.
    For more info:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#normalize12
    """
    norm12 = spm.Normalize12(jobtype='write')
    norm12.inputs.deformation_file  = trans_field
    norm12.inputs.image_to_align    = in_imgs
    norm12.inputs.write_voxel_sizes = [1, 1, 1]

    #norm12.run()
    return norm12


def spm_normalize(in_imgs=traits.Undefined, voxel_size=(1, 1, 1), template=None):
    """Return a Normalize12 interface object.
    SPM12’s new Normalise routine for warping an image to a template.

    Parameters
    ----------
    in_imgs: iterable of str

    voxel_size: iterable of 3 int
        Size of the output image voxel in mm in each dimension.
        Default: `(1, 1, 1)`

    template: str
        Path to the target registration template.
        Default: the SPM TPM file.

    Raises
    ------
    FileNotFoundError
        If `template` is `None` and can't find the TPM.nii file from SPM.

    Notes
    -----
    For more info:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#normalize12
    """
    if template is None:
        template = spm_tpm_priors_path()

    norm12 = spm.Normalize12(jobtype='estwrite',
                             tpm=template,
                             write_voxel_sizes=voxel_size)

    #norm12.run()
    return norm12


def spm_coregister(src_img=traits.Undefined, tgt_img=traits.Undefined, cost_function='mi'):
    """Use spm_coreg for estimating cross-modality rigid body alignment.

    More info: http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=39
    Parameters
    ----------
    src_img: str
        Path to the coregistration source image

    tgt_img: str
        Path to the coregistration target image

    cost_function: ('mi' or 'nmi' or 'ecc' or 'ncc')
        cost function, one of: 'mi' - Mutual Information,
         'nmi' - Normalised Mutual Information,
         'ecc' - Entropy Correlation Coefficient,
         'ncc' - Normalised Cross Correlation

        For more info:
        http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#coregister

    Returns
    -------
    coreg: nipype.interfaces.smp.Coregister
        spm.Coregister interface object
    """
    coreg = spm.Coregister()

    coreg.inputs.source = src_img
    coreg.inputs.target = tgt_img

    coreg.inputs.cost_function = cost_function
    coreg.inputs.jobtype = 'estwrite'

    #coreg.run()
    return coreg


