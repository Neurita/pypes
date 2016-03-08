# -*- coding: utf-8 -*-

from .environ import  spm_tpm_priors_path
from .files   import  (remove_ext,
                       get_extension,
                       get_affine,
                       get_data_dims,
                       get_vox_dims,
                       )

from .piping  import  (setup_node,
                       extend_trait_list,
                       fsl_merge,
                       joinstrings,
                       find_wf_node,
                       get_datasink,
                       get_input_node,
                       get_input_file_name,
                       )

from .spatial import get_bounding_box