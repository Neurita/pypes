# -*- coding: utf-8 -*-
"""
fMRI group template registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io
from   nipype.interfaces.utility import IdentityInterface

from   ..preproc import (spm_create_group_template_wf,
                         spm_register_to_template_wf)
from   ..config  import setup_node
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def attach_spm_fmri_grouptemplate(main_wf, wf_name="spm_fmri_template"):
    """ Attach a fMRI pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow picks all spm_fmri_preproc outputs 'fmri_output.warped_files' in `main_wf`
    to create a group template.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    fmri_output.warped_files: input node

    datasink: nipype Node

    spm_fmri_preproc: nipype Workflow

    Nipype Outputs
    --------------
    group_template.fmri_template: file
        The path to the PET group template.

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_fmri_preproc
    - spm_anat_preproc if `spm_fmri_template.do_fmripvc` is True.

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    fmri_wf = main_wf.get_node("spm_rest_preproc")

    in_files = get_input_node(main_wf)
    datasink = get_datasink(main_wf, name='datasink')

    # The base name of the 'fmri' file for the substitutions
    fmri_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))

    # the group template datasink
    base_outdir  = datasink.inputs.base_directory
    grp_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,), name="group_datasink")
    grp_datasink.inputs.container = '{}_grouptemplate'.format(fmri_fbasename)

    # the group template workflow
    template_wf = spm_create_group_template_wf(wf_name)

    # the list of the slicetime corrected fmri subjects
    warped_fmris = pe.JoinNode(interface=IdentityInterface(fields=["warped_fmris"]),
                               joinsource="warped_fmris",
                               joinfield="warped_fmris",
                               name="warped_fmris")

    # output node
    output = setup_node(IdentityInterface(fields=["fmri_template"]), name="group_template")

    # group dataSink output substitutions
    regexp_subst = [
                     (r"/wgrptemplate{rest}_merged_mean_smooth.nii$",  "/{rest}_grouptemplate_mni.nii"),
                     (r"/w{rest}_merged_mean_smooth.nii$",             "/{rest}_grouptemplate_mni.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, fmri=fmri_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    grp_datasink.inputs.regexp_substitutions = extend_trait_list(grp_datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # warped fmris file list input
                     (fmri_wf,       warped_fmris, [("fmri_output.warped_files",     "warped_fmris")]),

                     # group template wf
                     (warped_fmris,  template_wf, [("warped_fmris",                 "grptemplate_input.in_files")]),

                     # output node
                     (template_wf, output,       [("grptemplate_output.template", "fmri_template")]),

                     # template output
                     (output,      grp_datasink, [("fmri_template",                "@fmri_group_template")]),
                   ])

    # Now we start with the correction and registration of each subject to the group template
    # add the fmri template to the preproc workflow
    reg_wf = spm_register_to_template_wf(wf_name="spm_fmri_register_to_grouptemplate")
    main_wf.connect([
                     (output,   reg_wf,  [("fmri_template",  "reg_input.template")]),
                     (in_files, reg_wf,  [("fmri",           "reg_input.in_file"),]),

                     (reg_wf,   datasink, [("reg_output.warped",     "mrfmri.@warped"),
                                           ("reg_output.warp_field", "mrfmri.@warp_field"),
                                          ]),
                     ])

    # per-subject datasink output substitutions
    regexp_subst = [
                     (r"/{fmri}_sn.mat$",           "/{fmri}_grptemplate_params.mat"),
                     (r"/wgrptemplate_{fmri}.nii$", "/{fmri}_grptemplate.nii"),
                     (r"/w{fmri}.nii",              "/{fmri}_grptemplate.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, fmri=fmri_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    return main_wf

