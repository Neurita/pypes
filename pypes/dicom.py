# -*- coding: utf-8 -*-
"""
DICOM to Nifti workflows.

Here I implement 2 possible options:
one based on `dcm2nii` and the other based on `dcmstack`
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.dcm2nii as dcm2nii
import nipype.interfaces.dcmstack as dcmstack
from .io import build_crumb_workflow


def dcm2nii_wf(data_crumb, output_dir, cache_dir='', **kwargs):
    """ Returns a workflow for the a clinical database.

    Parameters
    ----------
    data_crumb: hansel.Crumb
        The crumb until the dicom files files.
        Example: Crumb('/home/hansel/data/{subject_id}/{session_id}/{acquisition}/{dicom_file})
        The last crumb argument of `data_crumb` must be '{image}', which indicates each of the
        subject/session files. This argument will be replaced by the corresponding image name.

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path

    kwargs: keyword arguments
        Keyword arguments with values for the data_crumb crumb path.
    """
    if kwargs:
        data_crumb = data_crumb.replace(**kwargs)

    if not data_crumb.exists():
        raise IOError("Expected an existing folder for `data_crumb`, got {}.".format(data_crumb))

    wfs = {"spm_anat_preproc": attach_spm_anat_preprocessing,
           "spm_mrpet_preproc": attach_spm_mrpet_preprocessing,
          }

    if wf_name not in wfs:
        raise ValueError("Expected `wf_name` to be in {}, got {}.".format(list(wfs.keys()),
                                                                          wf_name))

    if not cache_dir:
        cache_dir = op.join(op.dirname(output_dir), "wd")

    # generate the workflow
    main_wf = in_out_crumb_wf(work_dir=cache_dir,
                              data_crumb=data_crumb,
                              output_dir=output_dir,
                              crumb_arg_values=dict(**kwargs),
                              files_crumb_args={'anat': [('image', 'anat_hc.nii.gz')],
                                                'pet':  [('image', 'pet_fdg.nii.gz')],
                                                'diff': [('image', 'diff.nii.gz')],
                                                'bval': [('image', 'diff.bval')],
                                                'bvec': [('image', 'diff.bvec')],
                                               },
                              input_wf_name='input_files')

    wf = wfs[wf_name](main_wf=main_wf)

    # move the crash files folder elsewhere
    wf.config["execution"]["crashdump_dir"] = op.join(wf.base_dir, wf.name, "log")

    return wf
