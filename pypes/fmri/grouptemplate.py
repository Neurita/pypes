# -*- coding: utf-8 -*-
"""
fMRI group template registration nipype workflow.
"""
import os.path as op

import nipype.pipeline.engine as pe
from   nipype.interfaces.io import DataSink
from   nipype.interfaces import IdentityInterface

from   ..preproc import spm_create_group_template_wf, spm_warp_to_mni
from   ..config  import setup_node
from   .._utils  import format_pair_list
from   ..utils   import (get_datasink,
                         extend_trait_list,
                         get_input_node,
                         remove_ext,
                         get_input_file_name,
                         extension_duplicates)


def attach_spm_fmri_grouptemplate_wf(main_wf, wf_name='spm_epi_grouptemplate'):
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
    - fmri_cleanup: for the `rest_output.avg_epi` output

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    fmri_cleanup_wf = main_wf.get_node('fmri_cleanup')

    in_files = get_input_node(main_wf)
    datasink = get_datasink(main_wf, name='datasink')

    # The base name of the 'rest' file for the substitutions
    fmri_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'rest')))

    # the group template datasink
    base_outdir  = datasink.inputs.base_directory
    grp_datasink = pe.Node(DataSink(parameterization=False,
                                    base_directory=base_outdir,),
                                    name='{}_grouptemplate_datasink'.format(fmri_fbasename))
    grp_datasink.inputs.container = '{}_grouptemplate'.format(fmri_fbasename)

    # the list of the average EPIs from all the subjects
    #avg_epi_map = pe.MapNode(IdentityInterface(fields=['avg_epis']), iterfield=['avg_epis'], name='avg_epi_map')

    avg_epis = pe.JoinNode(IdentityInterface(fields=['avg_epis']), joinsource='infosrc', joinfield='avg_epis',
                           name='avg_epis')

    # directly warp the avg EPI to the SPM standard template
    warp_epis = spm_warp_to_mni("spm_warp_avgepi_to_mni")

    # the group template workflow
    template_wf = spm_create_group_template_wf(wf_name)

    # output node
    output = setup_node(IdentityInterface(fields=['fmri_template']), name='group_template')

    # group dataSink output substitutions
    regexp_subst = [
                     (r'/wgrptemplate{fmri}_merged_mean_smooth.nii$',  '/{fmri}_grouptemplate_mni.nii'),
                     (r'/w{fmri}_merged_mean_smooth.nii$',             '/{fmri}_grouptemplate_mni.nii'),
                   ]
    regexp_subst = format_pair_list(regexp_subst, fmri=fmri_fbasename)
    regexp_subst += extension_duplicates(regexp_subst)
    grp_datasink.inputs.regexp_substitutions = extend_trait_list(grp_datasink.inputs.regexp_substitutions,
                                                                 regexp_subst)

    # Connect the nodes
    main_wf.connect([
                     # the avg EPI inputs
                     (fmri_cleanup_wf, avg_epis,  [('rest_output.avg_epi',         'avg_epis')]),
                     #(fmri_cleanup_wf, avg_epi_map,  [('rest_output.avg_epi',         'avg_epis')]),
                     #(avg_epi_map,     avg_epis,     [('avg_epis',                    'avg_epis')]),

                     # warp avg EPIs to MNI
                     (avg_epis,        warp_epis,    [('avg_epis',                    'warp_input.in_files')]),

                     # group template wf
                     (warp_epis,       template_wf,  [('warp_output.warped_files',    'grptemplate_input.in_files')]),

                     # output node
                     (template_wf,     output,       [('grptemplate_output.template', 'fmri_template')]),

                     # template output
                     (output,          grp_datasink, [('fmri_template',               '@fmri_group_template')]),
                     (warp_epis,       grp_datasink, [('warp_output.warped_files',    'individuals.@warped')]),
                   ])

    return main_wf

