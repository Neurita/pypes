# -*- coding: utf-8 -*-

from .config import PYPES_CFG as configuration
from .config import (get_config_setting,
                     check_mandatory_inputs,
                     node_settings,
                     setup_node,
                     check_atlas_file,
                     update_config)
from .environ import  spm_tpm_priors_path
from .files   import  (remove_ext,
                       get_extension,
                       get_affine,
                       get_data_dims,
                       get_vox_dims,
                       )
from .piping  import  (extend_trait_list,
                       fsl_merge,
                       joinstrings,
                       find_wf_node,
                       get_datasink,
                       get_input_node,
                       get_input_file_name,
                       )



