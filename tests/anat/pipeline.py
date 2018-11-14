#!/usr/bin/env python

import os
import socket
from configparser import ConfigParser

from hansel import Crumb
from nipype import config as nipype_config

from neuro_pypes.io import build_crumb_workflow
from neuro_pypes.config import update_config
## Import functions
from neuro_pypes.anat import (
    attach_spm_anat_preprocessing,
)
from neuro_pypes.fmri import (
    attach_rest_preprocessing,
)
from neuro_pypes.pet import (
    attach_spm_mrpet_preprocessing,
    attach_spm_pet_preprocessing,
    attach_spm_pet_grouptemplate
)


if __name__ == "__main__":
    STDB_DIR = os.path.expanduser('~/projects/std_brains')
    SPM_DIR  = os.path.expanduser('~/software/spm12')
    BASE_DIR = os.path.expanduser('/data')
    data_dir = os.path.join(BASE_DIR, 'raw')
    cache_dir = os.path.join(BASE_DIR, 'wd')
    output_dir = os.path.join(BASE_DIR, 'out')

    HAMM_DIR = os.path.join(STDB_DIR, 'atlases', 'hammers')
    HAMM_MNI = os.path.join(HAMM_DIR, 'Hammers_mith_atlas_n30r83_SPM5.nii.gz')
    HAMM_LABELS = os.path.join(HAMM_DIR, 'labels.txt')

    SPM_CANONICAL_BRAIN_2MM = os.path.join(STDB_DIR, 'templates', 'spm_canonical', 'single_subj_T1_brain.nii.gz')

    # template files
    PET_MNI  = os.path.join(STDB_DIR, 'templates', 'pet.nii')
    MNI_MASK = os.path.join(STDB_DIR, 'templates', 'MNI152_T1_2mm_brain_mask.nii.gz')

    # data input/output os.path.dirname(__file__) means alongside with this .py file
    settings_file = os.path.join(os.path.dirname(__file__), 'config.yml')
    nipype_cfg_file = os.path.join(os.path.dirname(__file__), 'nipype.cfg')

    # passing configuration to pypes
    update_config(settings_file)
    update_config({'atlas_file': HAMM_MNI})

    # update nipype config
    config_data = ConfigParser()
    config_data.read(nipype_cfg_file)
    nipype_config.update_config(config_data)

    #Crumb path definition in a dictionary
    workflow_functions = dict([
        ("spm_anat_preproc", attach_spm_anat_preprocessing),
    ])

    data_path = os.path.join(os.path.expanduser(data_dir), '{subject_id}', '{session}', '{image}')
    data_crumb = Crumb(data_path, ignore_list=['.*'])

    crumb_modalities = {
        'anat': [('scan', 'T1'), ('image', 'Head_MPRAGE_highContrast.nii.gz')],
    }

    # put pypes workflow with crumb argunments into nipype
    workflow = build_crumb_workflow(
        workflow_functions,
        data_crumb,
        crumb_modalities,
        output_dir,
        cache_dir,
        "spm_anat" # final node name
    )
    workflow.run(plugin=None)
