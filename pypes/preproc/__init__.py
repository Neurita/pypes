# -*- coding: utf-8 -*-

from .petpvc import PETPVC
from .realign import nipy_motion_correction
from .denoise import (nlmeans_denoise_img,
                      create_regressors,
                      extract_noise_components,
                      motion_regressors,
                      reslice_img)
from .registration import (spm_apply_deformations,
                           spm_coregister,
                           spm_normalize,
                           spm_warp_to_mni,
                           afni_deoblique,
                           spm_tpm_priors_path,
                           spm_create_group_template_wf,
                           spm_register_to_template_wf,
                           )
from .slicetime import (afni_slicetime,
                        spm_slicetime,
                        auto_spm_slicetime,
                        auto_nipy_slicetime)
from .slicetime_params import (STCParameters,
                               STCParametersInterface)

from .spatial import get_bounding_box
