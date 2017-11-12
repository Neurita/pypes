# -*- coding: utf-8 -*-

from .environ import  spm_tpm_priors_path

from .files   import  (remove_ext,
                       rename,
                       get_extension,
                       get_affine,
                       get_data_dims,
                       get_vox_dims,
                       fetch_one_file,
                       extension_duplicates,)

from .piping  import  (extend_trait_list,
                       selectindex,
                       get_trait_value,
                       fsl_merge,
                       get_node,
                       joinstrings,
                       find_wf_node,
                       get_datasink,
                       get_input_node,
                       get_interface_node,
                       get_input_file_name,)

from .pandas import (add_table_headers,
                     write_tabbed_excel,)

