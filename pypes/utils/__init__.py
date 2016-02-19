# -*- coding: utf-8 -*-

from .environ import  spm_tpm_priors_path
from .files   import  remove_ext, get_extension
from .piping  import  (extend_trait_list,
                       fsl_merge,
                       joinstrings,
                       find_wf_node,
                       get_datasink,
                       get_input_node,
                       get_input_file_name,
                       )