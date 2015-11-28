"""
Functions to create pipelines for public and not so public available datasets.
"""

import os.path as op

from   .run  import in_out_workflow
from   .anat import attach_t1_preprocessing
from   .pet  import attach_pet_preprocessing


def cobre_workflow(wf_name, base_dir, cache_dir, output_dir, year):
    """ Returns a workflow for the COBRE database.

    Parameters
    ----------
    wf_name: str
        A name for the workflow.

    base_dir: str
        The folder path where the raw data is.

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path

    year: str or int
        The year of the subject set.
    """

    data_dir = base_dir
    if not data_dir or not op.exists(data_dir):
        raise IOError("Expected an existing folder for `data_dir`, got {}.".format(data_dir))

    wfs = {"spm_t1_preproc": attach_t1_preprocessing,
           # TODO: "spm_rest_preproc": attach_rest_preprocessing,
          }

    if wf_name not in wfs:
        raise ValueError("Expected `wf_name` to be in {}, got {}.".format(list(wfs.keys()),
                                                                          wf_name))

    # check some args
    if not output_dir:
        output_dir = op.join(op.dirname(data_dir), "out")

    if not cache_dir:
        cache_dir = op.join(op.dirname(data_dir), "wd")

    # generate the workflow
    main_wf = in_out_workflow(work_dir=cache_dir,
                              data_dir=data_dir,
                              output_dir=output_dir,
                              session_names=['session_1'],
                              file_names=['anat_1/mprage.nii.gz', 'rest_1/rest.nii.gz'],
                              subject_ids=None,
                              input_wf_name='input_files')

    wf = wfs[wf_name](main_wf=main_wf,
                      data_dir=data_dir,
                      work_dir=cache_dir,
                      output_dir=output_dir,)

    # move the crash files folder elsewhere
    wf.config["execution"]["crashdump_dir"] = op.join(wf.base_dir, wf.name, "log")

    return wf


def clinical_dataset(wf_name, base_dir, cache_dir, output_dir, year):
    """ Run an specific pipeline.

    Parameters
    ----------
    wf_name: str
        A name for the workflow.

    base_dir: str
        The folder path where the raw data is.

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path

    year: str or int
        The year of the subject set.
    """
    if not year:
        data_dir = base_dir
    else:
        data_dir = op.join(base_dir, year)

    if not data_dir or not op.exists(data_dir):
        raise IOError("Expected an existing folder for `data_dir`, got {}.".format(data_dir))

    wfs = {"spm_t1_preproc": attach_t1_preprocessing,
           "spm_pet_preproc": attach_pet_preprocessing,
          }

    if wf_name not in wfs:
        raise ValueError("Expected `wf_name` to be in {}, got {}.".format(list(wfs.keys()),
                                                                          wf_name))

    # check some args
    if not output_dir:
        output_dir = op.join(op.dirname(data_dir), "out", year)

    if not cache_dir:
        cache_dir = op.join(op.dirname(data_dir), "wd", year)

    # generate the workflow
    main_wf = in_out_workflow(work_dir=cache_dir,
                              data_dir=data_dir,
                              output_dir=output_dir,
                              session_names=['session_0'],
                              file_names=['anat_hc.nii.gz', 'pet_fdg.nii.gz'],
                              subject_ids=None,
                              input_wf_name='input_files')

    wf = wfs[wf_name](main_wf=main_wf,
                      data_dir=data_dir,
                      work_dir=cache_dir,
                      output_dir=output_dir,)

    # move the crash files folder elsewhere
    wf.config["execution"]["crashdump_dir"] = op.join(wf.base_dir, wf.name, "log")

    return wf
