# -*- coding: utf-8 -*-
"""
Nipype workflows to co-register anatomical MRI to diffusion MRI.
"""
import nipype.pipeline.engine as pe
from   nipype.interfaces.fsl import MultiImageMaths
from   nipype.interfaces.utility import IdentityInterface, Select, Split
from   nipype.algorithms.misc import Gunzip

from .._utils  import flatten_list
from ..config  import setup_node, check_atlas_file
from ..preproc import spm_coregister


def spm_anat_to_diff_coregistration(wf_name="spm_anat_to_diff_coregistration"):
    """ Co-register the anatomical image and other images in anatomical space to
    the average B0 image.

    This estimates an affine transform from anat to diff space, applies it to
    the brain mask and an atlas.

    Nipype Inputs
    -------------
    dti_co_input.avg_b0: traits.File
        path to the average B0 image from the diffusion MRI.
        This image should come from a motion and Eddy currents
        corrected diffusion image.

    dti_co_input.anat: traits.File
        path to the high-contrast anatomical image.

    dti_co_input.tissues: traits.File
        paths to the NewSegment c*.nii output files, in anatomical space

    dti_co_input.atlas_anat: traits.File
        Atlas in subject anatomical space.

    Nipype Outputs
    --------------
    dti_co_output.anat_diff: traits.File
        Anatomical image in diffusion space.

    dti_co_output.tissues_diff: traits.File
        Tissues images in diffusion space.

    dti_co_output.brain_mask_diff: traits.File
        Brain mask for diffusion image.

    dti_co_output.atlas_diff: traits.File
        Atlas image warped to diffusion space.
        If the `atlas_file` option is an existing file and `normalize_atlas` is True.

    Nipype Workflow Dependencies
    ----------------------------
    This workflow depends on:
    - spm_anat_preproc

    Returns
    -------
    wf: nipype Workflow
    """
    # specify input and output fields
    in_fields  = ["avg_b0", "tissues", "anat"]
    out_fields = ["anat_diff",
                  "tissues_diff",
                  "brain_mask_diff",]

    do_atlas, _ = check_atlas_file()
    if do_atlas:
        in_fields  += ["atlas_anat"]
        out_fields += ["atlas_diff"]

    # input interface
    dti_input = pe.Node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                        name="dti_co_input")

    gunzip_b0 = pe.Node(Gunzip(), name="gunzip_b0")
    coreg_b0  = setup_node(spm_coregister(cost_function="mi"), name="coreg_b0")

    # co-registration
    brain_sel    = pe.Node(Select(index=[0, 1, 2]),            name="brain_sel")
    coreg_split  = pe.Node(Split(splits=[1, 2], squeeze=True), name="coreg_split")

    brain_merge  = setup_node(MultiImageMaths(), name="brain_merge")
    brain_merge.inputs.op_string = "-add '%s' -add '%s' -abs -kernel gauss 4 -dilM -ero -kernel gauss 1 -dilM -bin"
    brain_merge.inputs.out_file = "brain_mask_diff.nii.gz"

    # output interface
    dti_output = pe.Node(IdentityInterface(fields=out_fields),
                         name="dti_co_output")

    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # Connect the nodes
    wf.connect([
                # co-registration
                (dti_input, coreg_b0, [("anat", "source")]),

                (dti_input,     brain_sel,   [("tissues",             "inlist")]),
                (brain_sel,     coreg_b0,    [(("out", flatten_list), "apply_to_files")]),

                (dti_input,     gunzip_b0,   [("avg_b0",   "in_file")]),
                (gunzip_b0,     coreg_b0,    [("out_file", "target")]),

                (coreg_b0,      coreg_split, [("coregistered_files", "inlist")]),
                (coreg_split,   brain_merge, [("out1",               "in_file")]),
                (coreg_split,   brain_merge, [("out2",               "operand_files")]),

                # output
                (coreg_b0,     dti_output,     [("coregistered_source", "anat_diff")]),
                (coreg_b0,     dti_output,     [("coregistered_files",  "tissues_diff")]),
                (brain_merge,  dti_output,     [("out_file",            "brain_mask_diff")]),
              ])

    # add more nodes if to perform atlas registration
    if do_atlas:
        coreg_atlas = setup_node(spm_coregister(cost_function="mi"), name="coreg_atlas")

        # set the registration interpolation to nearest neighbour.
        coreg_atlas.inputs.write_interp = 0
        wf.connect([
                    (dti_input,   coreg_atlas, [("anat",               "source"),
                                                ("atlas_anat",         "apply_to_files"),
                                               ]),
                    (gunzip_b0,   coreg_atlas, [("out_file",           "target")]),
                    (coreg_atlas, dti_output,  [("coregistered_files", "atlas_diff")]),
                  ])

    return wf

