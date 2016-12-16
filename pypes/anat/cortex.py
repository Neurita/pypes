# -*- coding: utf-8 -*-
"""
Nipype workflow for cortical thickness measures using ANTs.
"""
import os.path as op

from   nipype.interfaces.ants import CorticalThickness

from   ..config  import setup_node
from   .._utils  import format_pair_list
from   ..utils   import (remove_ext,
                         extend_trait_list,
                         get_input_node,
                         get_datasink,
                         get_input_file_name,
                         extension_duplicates)


def attach_ants_cortical_thickness(main_wf, wf_name='ants_cortical_thickness'):
    """ Attach the ANTs anatomical cortical thickness workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Not used, this is there only to keep a uniform function signature.

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.anat: input node

    cortical_thickness.brain_template: traits.File
        Anatomical *intensity* template (possibly created using a population
        data set with buildtemplateparallel.sh in ANTs).
        This template is *not* skull-stripped.

    cortical_thickness.segmentation_priors: string
        Pattern for paths to the tissue segmentation priors in the space of the
        given template.

        Specified using c-style formatting, e.g. "labelsPriors%02d.nii.gz".
        We assume that the first four priors are ordered as follows
            1:  csf
            2:  cortical gm
            3:  wm
            4:  deep gm

    corthick_input.brain_probability_mask: traits.File
        Brain probability mask in template space.

    cortical_thickness.t1_registration_template:
        Anatomical *intensity* template(assumed to be skull-stripped). A
        commoncase would be where this would be the sametemplate as
        specified in the -e option whichis not skull stripped.

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # cortical thickness node
    cort_thick = setup_node(CorticalThickness(), name="cortical_thickness")

    # The base name of the 'anat' file for the substitutions
    anat_fbasename = remove_ext(op.basename(get_input_file_name(in_files, 'anat')))

    # dataSink output substitutions
    regexp_subst = [
                   ]
    regexp_subst = format_pair_list(regexp_subst, anat=anat_fbasename)

    # add nii.gz patterns
    regexp_subst += extension_duplicates(regexp_subst)
    datasink.inputs.regexp_substitutions = extend_trait_list(datasink.inputs.regexp_substitutions,
                                                             regexp_subst)

    main_wf.connect([(in_files, cort_thick, [("anat", "anatomical_image")]),

                     # output
                     (cort_thick, datasink, [("BrainExtractionMask",                "anat.cortex.@brain_mask"),
                                             ("BrainSegmentation",                  "anat.cortex.@brain_segmentation"),
                                             ("BrainSegmentationN4",                "anat.cortex.@brain_n4"),
                                             ("BrainSegmentationPosteriors",        "anat.cortex.@brain_seg_posteriors"),
                                             ("BrainVolumes",                       "anat.cortex.@brain_volumes"),
                                             ("CorticalThickness",                  "anat.cortex.@cortical_thickness"),
                                             ("CorticalThicknessNormedToTemplate",  "anat.cortex.@normed_cortical_thickness"),
                                             ("SubjectToTemplate0GenericAffine",    "anat.cortex.@subj_to_template_affine"),
                                             ("SubjectToTemplate1Warp",             "anat.cortex.@subj_to_template_warp"),
                                             ("SubjectToTemplateLogJacobian",       "anat.cortex.@subj_to_template_logjac"),
                                             ("TemplateToSubject0Warp",             "anat.cortex.@template_to_subj_warp"),
                                             ("TemplateToSubject1GenericAffine",    "anat.cortex.@template_to_subj_affine"),
                                            ]),
                    ])

    return main_wf
