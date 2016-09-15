# -*- coding: utf-8 -*-
"""
fMRI group template registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io
from   nipype.interfaces.utility import IdentityInterface

from   ..preproc import spm_create_group_template_wf
from   ..config  import setup_node
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def attach_spm_fmri_grouptemplate_wf(main_wf, wf_name="spm_fmri_template"):
    """ Attach a fMRI pre-processing workflow that uses SPM12 to `main_wf`.
    This workflow picks all spm_fmri_preproc outputs 'fmri_output.warped_files' in `main_wf`
    to create a group template.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow


    Nipype Inputs
    -------------
    rest_input.in_file: traits.File
        The slice time and motion corrected fMRI file.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    rest_output.avg_epi_mni: input node

    datasink: nipype Node

    spm_rest_preproc_mni: nipype Workflow

    Nipype Outputs
    --------------
    group_template.fmri_template: file
        The path to the fMRI group template.

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_fmri_preproc

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    fmri_warp_wf = main_wf.get_node("spm_warp_fmri_mni")

    in_files = get_input_node(main_wf)
    datasink = get_datasink(main_wf, name='datasink')

    # The base name of the 'rest' file for the substitutions
    fmri_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))

    # the group template datasink
    base_outdir  = datasink.inputs.base_directory
    grp_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,), name="group_datasink")
    grp_datasink.inputs.container = '{}_grouptemplate'.format(fmri_fbasename)

    # the list of the average EPIs from each subject
    warped_epis = pe.JoinNode(interface=IdentityInterface(fields=["warped_epis"]),
                              joinsource="infosrc",
                              joinfield="warped_epis",
                              name="warped_epis")

    # the group template workflow
    template_wf = spm_create_group_template_wf(wf_name)

    # output node
    output = setup_node(IdentityInterface(fields=["fmri_template"]), name="group_template")

    # group dataSink output substitutions
    regexp_subst = [
                     (r"/wgrptemplate{fmri}_merged_mean_smooth.nii$",  "/{fmri}_grouptemplate_mni.nii"),
                     (r"/w{fmri}_merged_mean_smooth.nii$",             "/{fmri}_grouptemplate_mni.nii"),
                   ]
    regexp_subst = format_pair_list(regexp_subst, fmri=fmri_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    grp_datasink.inputs.regexp_substitutions = extend_trait_list(grp_datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # warped fmris file list input
                     (fmri_warp_wf,  warped_epis, [("wfmri_output.wavg_epi",        "warped_epis")]),

                     # group template wf
                     (warped_epis,  template_wf,  [("warped_epis",                  "grptemplate_input.in_files")]),

                     # output node
                     (template_wf, output,        [("grptemplate_output.template",  "fmri_template")]),

                     # template output
                     (output,      grp_datasink,  [("fmri_template",                "@fmri_group_template")]),
                   ])

    return main_wf

