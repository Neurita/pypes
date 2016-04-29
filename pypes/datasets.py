# -*- coding: utf-8 -*-
"""
Functions to create pipelines for public and not so public available datasets.
"""

from collections import OrderedDict

from   .config  import update_config
from   .anat    import attach_spm_anat_preprocessing
from   .dti     import (attach_fsl_dti_preprocessing,
                        attach_camino_tractography)
from   .fmri    import attach_rest_preprocessing
from   .io      import build_crumb_workflow
from   .pet     import (attach_spm_mrpet_preprocessing,
                        attach_spm_pet_preprocessing,
                        attach_spm_pet_grouptemplate)


def _cobre_wf_setup(wf_name):
    """ Return a list of workflow-attach functions and a dict with kwargs for
    the `in_out_crumb_wf` function to run workflows against the COBRE database.

    Parameters
    ----------
    wf_name

    Returns
    -------
    wf_attachers: dict with wf attachers

    in_out_wf_kwargs: dict with kwargs
    """
    attach_functions = {"spm_anat_preproc":         [("spm_anat_preproc",  attach_spm_anat_preprocessing)],

                        "spm_anat_rest_preproc":    [("spm_anat_preproc",  attach_spm_anat_preprocessing),
                                                     ("spm_rest_preproc",  attach_rest_preprocessing),
                                                     ],
                       }

    files_crumb_args = {'anat':  [('modality', 'anat_1'),
                                  ('image',    'mprage.nii.gz')]} #'anat_1/mprage.nii.gz',

    if 'rest' in wf_name:
        files_crumb_args.update({'rest':  [('modality', 'rest_1'),
                                           ('image',    'rest.nii.gz')], # 'rest_1/rest.nii.gz'},
                                })

    params = {'file_templates': files_crumb_args}

    return OrderedDict(attach_functions[wf_name]), params


def _clinical_wf_setup(wf_name):
    """ Return a list of workflow-attach functions and a dict with kwargs for
    the `in_out_crumb_wf` function to run workflows against the clinical database.

    Parameters
    ----------
    wf_name

    Returns
    -------
    wf_attachers: dict with wf attachers

    in_out_wf_kwargs: dict with kwargs
    """
    attach_functions = {"spm_anat_preproc":     [("spm_anat_preproc",  attach_spm_anat_preprocessing)],

                        "spm_pet_preproc":      [("spm_pet_preproc",   attach_spm_pet_preprocessing)],

                        "spm_pet_template":     [("spm_pet_preproc",       attach_spm_pet_preprocessing),
                                                 ("spm_pet_grouptemplate", attach_spm_pet_grouptemplate),
                                                ],

                        "spm_anat_pet_preproc": [("spm_anat_preproc",  attach_spm_anat_preprocessing),
                                                 ("spm_mrpet_preproc", attach_spm_mrpet_preprocessing),
                                                ],

                        "fsl_dti_preproc":      [("spm_anat_preproc",  attach_spm_anat_preprocessing),
                                                 ("fsl_dti_preproc",   attach_fsl_dti_preprocessing),
                                                ],

                        "anat_dti_camino":      [("spm_anat_preproc",  attach_spm_anat_preprocessing),
                                                 ("fsl_dti_preproc",   attach_fsl_dti_preprocessing),
                                                 ("camino_tract" ,     attach_camino_tractography),
                                                ],

                        "anat_pet_dti_camino":  [("spm_anat_preproc",  attach_spm_anat_preprocessing),
                                                 ("spm_mrpet_preproc", attach_spm_mrpet_preprocessing),
                                                 ("fsl_dti_preproc",   attach_fsl_dti_preprocessing),
                                                 ("camino_tract",      attach_camino_tractography),
                                                ],
                       }

    files_crumb_args = {}

    if 'anat' in wf_name:
        files_crumb_args.update({'anat':  [('image', 'anat_hc.nii.gz')]})

    if 'pet' in wf_name:
        files_crumb_args.update({'pet':  [('image', 'pet_fdg.nii.gz')],})

    if 'dti' in wf_name:
        files_crumb_args.update({'diff': [('image', 'diff.nii.gz')],
                                 'bval': [('image', 'diff.bval')],
                                 'bvec': [('image', 'diff.bvec')],
                                })

    params = {'file_templates': files_crumb_args}

    return OrderedDict(attach_functions[wf_name]), params


def cobre_crumb_workflow(wf_name, data_crumb, output_dir, cache_dir='', config_file='', params=None):
    """ Returns a workflow for the COBRE database.

    Parameters
    ----------
    wf_name: str
        A name for the workflow.
        Choices: 'spm_anat_preproc': MPRAGE preprocessing with SPM12
                 'spm_rest_preproc': MPRAGE+rs-fMRI preprocessing with SPM12 (not implemented yet)

    data_crumb: hansel.Crumb
        The crumb until the subject files.
        Example: Crumb('/home/hansel/cobre/raw/{subject_id}/session_1/{modality}/{image})
        The last 2 crumb arguments of `data_crumb` must be '{modality}/{image}', which indicates each of the
        subject/session files. This argument will be replaced by the corresponding image name.

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path

    config_file: str
        Path to a configuration file. Will use anything compatible with Kaptan.
        Have a look at config.py.

    params: dict with arguments
        crumb_replaces

        atlas_file
    """
    attach_funcs, in_out_wf_kwargs = _cobre_wf_setup(wf_name)

    if config_file:
        update_config(config_file)

    if params:
        update_config(params)

    wf = build_crumb_workflow(attach_funcs,
                              data_crumb=data_crumb,
                              in_out_kwargs=in_out_wf_kwargs,
                              output_dir=output_dir,
                              cache_dir=cache_dir,)

    return wf


def clinical_crumb_workflow(wf_name, data_crumb, output_dir, cache_dir='', config_file='', params=None):
    """ Returns a workflow for the a clinical database.

    Parameters
    ----------
    wf_name: str
        A name for the workflow to be created.
        Choices: 'spm_anat_preproc':     MPRAGE preprocessing with SPM12
                 'spm_anat_pet_preproc': MPRAGE+FDG-PET preprocessing with SPM12
                 'spm_pet_preproc':      FDG-PET only preprocessing with SPM12

    data_crumb: hansel.Crumb
        The crumb until the subject files.
        Example: Crumb('/home/hansel/data/{subject_id}/{session_id}/{modality}/{image})
        The last crumb argument of `data_crumb` must be '{image}', which indicates each of the
        subject/session files. This argument will be replaced by the corresponding image name.

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path

    config_file: str
        Path to a configuration file. Will use anything compatible with Kaptan.
        Have a look at config.py.

    params: dict with arguments
        crumb_replaces

        atlas_file

        raise_on_filenotfound
    """
    attach_funcs, in_out_wf_kwargs = _clinical_wf_setup(wf_name)

    if params is None:
        params = {}

    params.update(in_out_wf_kwargs)

    if config_file:
        update_config(config_file)

    if params:
        update_config(params)

    wf = build_crumb_workflow(attach_funcs,
                              data_crumb=data_crumb,
                              in_out_kwargs=in_out_wf_kwargs,
                              output_dir=output_dir,
                              cache_dir=cache_dir,)

    return wf


