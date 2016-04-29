# -*- coding: utf-8 -*-

from .pet_utils import petpvc_cmd, petpvc_mask, intensity_norm
from .petpvc import PETPVC
from .realign import nipy_motion_correction
from .denoise import (nlmeans_denoise_img,
                      create_regressors,
                      extract_noise_components,
                      motion_regressors)
from .registration import (spm_apply_deformations,
                           spm_coregister,
                           spm_normalize,
                           afni_deoblique,
                           spm_group_template,
                           spm_tpm_priors_path,
                          )
from .slicetime import (afni_slicetime,
                        spm_slicetime,
                        auto_spm_slicetime,
                        auto_nipy_slicetime)
from .slicetime_params import STCParametersInterface

from .spatial import get_bounding_box