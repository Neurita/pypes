# -*- coding: utf-8 -*-
"""
Nipype registration nodes and workflows
"""
import nipype.pipeline.engine    as pe
from   nipype.interfaces.utility import IdentityInterface, Function

import nipype.interfaces.spm  as spm
import nipype.interfaces.afni as afni
from   nipype.interfaces.base import traits

from ..nilearn import mean_img, concat_imgs, smooth_img
from ..config  import setup_node
from ..utils   import spm_tpm_priors_path


def spm_apply_deformations(in_imgs=traits.Undefined,
                           trans_field=traits.Undefined,
                           bbox=traits.Undefined):
    """Return a Normalize12 interface object.
    SPM12's new Normalise routine for warping an image to a template.
    For more info:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#normalize12

    Parameters
    ----------
    in_imgs: iterable of str

    trans_field: file
        file y_*.nii containing 3 deformation fields for the deformation in
        x, y and z dimension

    bbox: (a list of from 2 to 2 items which are a list of
           from 3 to 3 items which are a float)
        write_bounding_box option.
        3x2-element list of lists representing the bounding box (in mm) to
        be written
        See the preproc.get_bounding_box function.
        This is important to set when applying deformation fields
        to cropped images.

    Nipype Ouputs
    -------------
    deformation_field: (a list of items which are an existing file name)
        NIfTI file containing 3 deformation fields for the deformation in x,
        y and z dimension

    normalized_files: (a list of items which are an existing file name)
        Normalized other files

    normalized_image: (a list of items which are an existing file name)
        Normalized file that needed to be aligned

    """
    norm12 = spm.Normalize12(jobtype='write')

    norm12.inputs.deformation_file   = trans_field
    norm12.inputs.image_to_align     = in_imgs
    norm12.inputs.write_voxel_sizes  = [1, 1, 1]
    norm12.inputs.write_bounding_box = bbox

    #norm12.run()
    return norm12


def spm_normalize(in_imgs=traits.Undefined, template=None, **kwargs):
    """Return a Normalize12 interface object.
    SPM12's new Normalise routine for warping an image to a template.

    Parameters
    ----------
    in_imgs: iterable of str

    template: str
        Template in form of tissue probablitiy maps to normalize to
        mutually_exclusive: deformation_file.
        Default: the SPM TPM file.

    Nipype Ouputs
    -------------
    deformation_field: (a list of items which are an existing file name)
        NIfTI file containing 3 deformation fields for the deformation in x,
        y and z dimension

    normalized_files: (a list of items which are an existing file name)
        Normalized other files

    normalized_image: (a list of items which are an existing file name)
        Normalized file that needed to be aligned


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
                             **kwargs)

    #norm12.run()
    return norm12


def spm_coregister(src_img=traits.Undefined, tgt_img=traits.Undefined, cost_function='mi'):
    """Use spm_coreg for estimating cross-modality rigid body alignment.
    The write_interp option is set to 0 by default.

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


def afni_deoblique(in_file=traits.Undefined, out_file=traits.Undefined, out_type='NIFTI_GZ'):
    """ Return a nipype interface for AFNI '3dWarp -deoblique'.

    Parameters
    ----------
    in_file: str
        Path to the input file

    out_file: str
        Path to the output file.

    out_type: str
        ('NIFTI_GZ' or 'AFNI' or 'NIFTI')
        AFNI output filetype

    Returns
    -------
    deob: nipype.interfaces.afni.Warp
    """

    deob = afni.Warp()
    deob.inputs.in_file = in_file
    deob.inputs.deoblique = True
    deob.inputs.out_file = out_file
    deob.inputs.outputtype = out_type

    return deob


def spm_group_template(wf_name="spm_group_template"):
    """ Pick all subject files in `grptemplate_input.in_files`, calculate an average
    image and smooth it with `"{}_smooth".format(wf_name)` node (you can configure the smooth `fwhm` from
    a config file.).

    It does:
    - calculate a mean image (across subjects) and
    - smooth it with 8x8x8 gaussian kernel -> this is the template.
    - Finally, warp all PET images again to this template.
    - If tissue data is available from MR, do PVC.

    Parameters
    ----------
    wf_name: str
        Name of the workflow.

    Nipype Inputs
    -------------
    grptemplate_input.in_files: list of traits.File
        The raw NIFTI_GZ PET image files

    Nipype outputs
    --------------
    grptemplate_output.template: existing file
        The common custom PET template file.

    Returns
    -------
    wf: nipype Workflow
    """
    # input
    input = setup_node(IdentityInterface(fields=["in_files"]), name="grptemplate_input",)

    # merge
    concat = setup_node(Function(function=concat_imgs, input_names=["in_files"], output_names=["out_file"],
                                 imports=['from pypes.nilearn import ni2file']),
                        name='merge_time')

    # average
    average = setup_node(Function(function=mean_img, input_names=["in_files"], output_names=["out_file"],
                                 imports=['from pypes.nilearn import ni2file']),
                        name='average')

    # smooth
    smooth = setup_node(Function(function=smooth_img,
                                 input_names=["in_file", "fwhm"],
                                 output_names=["out_file"],
                                 imports=['from pypes.nilearn import ni2file']),
                         name="{}_smooth".format(wf_name))

    # output
    output = setup_node(IdentityInterface(fields=["template"]), name="grptemplate_output",)

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    wf.connect([
                # input
                (input,       concat,  [("in_files", "in_files")]),

                # merge, average and smooth
                (concat,      average, [("out_file", "in_file")]),
                (average,     smooth,  [("out_file", "in_file")]),

                # output
                (smooth,     output,   [("out_file", "template")]),
               ])

    return wf
