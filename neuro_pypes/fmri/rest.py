# -*- coding: utf-8 -*-
"""
Nipype workflows to process resting-state functional MRI.
"""

from .grouptemplate import attach_spm_fmri_grouptemplate_wf
from .clean import attach_fmri_cleanup_wf
from .warp import attach_spm_warp_fmri_wf


def _attach_rest_preprocessing(main_wf, registration_wf_name="spm_warp_fmri", do_group_template=False):
    """ Attach the resting-state MRI pre-processing workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    registration_wf_name: str
         Name of the registration workflow.

    do_group_template: bool
        If True will attach the group template creation and pre-processing pipeline.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.select.anat: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    main_wf = attach_fmri_cleanup_wf(main_wf)
    main_wf = attach_spm_warp_fmri_wf(main_wf,
                                      registration_wf_name=registration_wf_name,
                                      do_group_template=False)

    if do_group_template:
        main_wf = attach_spm_fmri_grouptemplate_wf(main_wf, wf_name="spm_fmri_grptemplate")
        main_wf = attach_spm_warp_fmri_wf(main_wf,
                                          registration_wf_name=registration_wf_name,
                                          do_group_template=True)

        reg_wf       = main_wf.get_node("{}_{}".format(registration_wf_name, 'grptemplate'))
        grp_template = main_wf.get_node("group_template")
        main_wf.connect([(grp_template, reg_wf,  [("fmri_template",  "wfmri_input.epi_template")]),])

    return main_wf


def attach_rest_preprocessing(main_wf, wf_name="spm_warp_fmri"):
    return _attach_rest_preprocessing(main_wf,
                                      registration_wf_name=wf_name,
                                      do_group_template=False)


def attach_rest_grptemplate_preprocessing(main_wf, wf_name="spm_warp_fmri"):
    return _attach_rest_preprocessing(main_wf,
                                      registration_wf_name=wf_name,
                                      do_group_template=True)
