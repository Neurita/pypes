# -*- coding: utf-8 -*-
"""
DICOM to Nifti converter node based on dcm2niix and hansel.Crumb.
"""
import nipype.pipeline.engine as pe
from nipype.interfaces.base    import traits
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.dcm2nii import Dcm2niix, Dcm2nii

from ..config import setup_node
from ..utils  import (get_datasink,
                      get_input_node)


def dcm2niix_wf(wf_name='dcm2niix'):
    """Run dcm2niix over one folder with DICOM files.

    Nipype Inputs
    -------------
    dcm2niix.in_dcmdir: traits.Dir
        path to the DICOM images folder.

    Nipype Outputs
    --------------
    dcm2niix.bids: (a list of items which are an existing file name)

    dcm2niix.bvals: (a list of items which are an existing file name)

    dcm2niix.bvecs: (a list of items which are an existing file name)

    dcm2niix.converted_files: (a list of items which are an existing file name)

    Returns
    -------
    wf: nipype Workflow
    """
    # Create the workflow object
    wf = pe.Workflow(name=wf_name)

    # specify input and output fields
    in_fields  = ["in_dcmdir",]

    out_fields = ["bids",
                  "bvals",
                  "bvecs",
                  "converted_files"]

    # input node
    dcm2niix_input = setup_node(IdentityInterface(fields=in_fields,
                                                  mandatory_inputs=True),
                                name="dcm2niix_input")
    # T1 preprocessing nodes
    dcm2niix = setup_node(Dcm2niix(), name="dcm2niix")

    # output node
    dcm2niix_output = setup_node(IdentityInterface(fields=out_fields),
                                 name="dcm2niix_output")

    # Connect the nodes
    wf.connect([
                # input
                (dcm2niix_input, dcm2niix, [("in_dcmdir",  "source_dir"),
                                           ]),

                # output
                (dcm2niix, dcm2niix_output, [("bids",            "bids"),
                                             ("bvals",           "bvals"),
                                             ("bvecs",           "bvecs"),
                                             ("converted_files", "converted_files"),
                                            ]),
              ])

    return wf


def attach_dcm2niix(main_wf, wf_name="dcm2niix"):
    """ Attach the dcm2niix workflow to the `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the dcm2niix workflow

    Nipype Inputs for `main_wf`
    ---------------------------
    Note: The `main_wf` workflow is expected to have an `input_files` and a `datasink` nodes.

    input_files.dcm_dir: input node

    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    in_files = get_input_node(main_wf)
    datasink = get_datasink  (main_wf)

    # The workflow box
    d2n_wf = dcm2niix_wf(wf_name=wf_name)

    main_wf.connect([(in_files, d2n_wf,   [("dcm_dir",                         "dcm2niix_input.in_dcmdir")]),
                     (d2n_wf,   datasink, [("dcm2niix_output.bids",            "@bids"),
                                           ("dcm2niix_output.bvals",           "@bvals"),
                                           ("dcm2niix_output.bvecs",           "@bvecs"),
                                           ("dcm2niix_output.converted_files", "@converted_files"),
                                          ],),
                    ])

    return main_wf


def dcm2nii_converter(source_names=traits.Undefined):
    """ Create a dcm2nii workflow.

    It does:
    - dcm2nii

    Nipype Inputs
    -------------
    source_names: traits.File
        path to the DICOM files or series folder.

    Returns
    -------
    wf: nipype Workflow

    Note
    ----
    For more info: http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.dcm2nii.html
    # TODO: an attach function for this node.
    """
    dcm2nii = Dcm2nii()

    dcm2nii.inputs.gzip_output     = True
    dcm2nii.inputs.output_dir      = '.'
    dcm2nii.inputs.terminal_output = 'file'
    dcm2nii.inputs.source_names    = source_names

    return dcm2nii

