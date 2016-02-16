# -*- coding: utf-8 -*-
"""
A nipype.SelectFiles node based on hansel.Crumb
"""

from   copy import deepcopy

import nipype.pipeline.engine as pe
from   nipype.interfaces.io import SelectFiles
from   nipype.interfaces.utility import IdentityInterface

from .utils.piping import get_values_map_keys


def crumb_input(work_dir, data_crumb, crumb_arg_values, files_crumb_args, wf_name="input_files"):
    """ A workflow of IdentityInterface->SelectFiles set up using a hansel.Crumb

    Parameters
    ----------
    work_dir: str
        Path to the working directory of the workflow

    data_crumb: hansel.Crumb
        The crumb until the subject files.
        Example: Crumb('/home/hansel/data/{subject_id}/{session_id}/{modality}/{image_file})

    crumb_arg_values: Dict[str -> list]
        Example 1: {'session_id': ['session_0', 'session_1', 'session_2'],
                    'subject_id': ['hansel', 'gretel'],
                   }

        This will be input to an IdentityInterface node.

        **Note**: if any crumb argument is not being defined here, will use all the unique values using
        the Crumb `ls` function. Note that this will only work if the structure tree is the same for all
        subjects.

    files_crumb_args: Dict[str -> list of 2-tuple]
        Maps of crumb argument values to specify each file in the `data_crumb`.
        Example: {'anat': [('modality', 'anat'), ('image_file', 'anat_hc.nii.gz')],
                  'pet':  [('modality', 'pet'),  ('image_file', ''pet_fdg.nii.gz'')],
                 }

    wf_name: str
        Name of the workflow

    Nipype Outputs
    --------------
    select.{template_key}: path to existing file
        Will give you the path of the {file_name}s taken from the keys of `files_crumb_args`.

    infosrc.{field_name}: str
        Will give you the value of the field from `crumb_arg_values` keys.

    Returns
    -------
    wf: nipype Workflow

    """
    # Input workflow
    wf = pe.Workflow(name=wf_name, base_dir=work_dir)

    # check iterables in crumb
    param_args = deepcopy(crumb_arg_values) if crumb_arg_values else {}

    # create the lists of argument names
    crumb_args = list(data_crumb.all_args())
    arg_names  = list(param_args.keys()) if param_args else []
    arg_names.extend(get_values_map_keys(files_crumb_args))

    # check their size and expand them if needed#
    # TODO: this has to be fixed to list only folders that exist in the complete path, not only all possible combinations.
    n_crumbs   = len(crumb_args)
    n_names    = len(arg_names)
    if n_crumbs > n_names: # add the missing argument values to the crumb_arg_values dictionary
        rem_args = set(crumb_args) - set(arg_names)
        for arg in rem_args:
            param_args[arg] = sorted(list(set(data_crumb[arg])))

    elif n_crumbs < n_names:
        raise KeyError('Expected `crumb_arg_values` to have the argument names of {}, '
                       'but got these extra: {}.'.format(data_crumb, set(arg_names)-set(crumb_args)))

    # Infosource - a function free node to iterate over the list of subject names
    field_names = list(data_crumb.all_args())

    infosource = pe.Node(IdentityInterface(fields=field_names), name="infosrc")
    infosource.iterables = list(param_args.items())

    # SelectFiles
    # another option is to use the '*' in the selectfiles path for the argument that are not being defined in
    # `crumb_arg_values` but then, the crumb ignore_list and regexes will be ignored.
    # for arg_name in field_names:
    #     if arg_name not in crumb_arg_values:
    #         data_crumb.replace(**dict([(arg_name, '*')]))

    select = pe.Node(SelectFiles({fkey: data_crumb.replace(**dict(farg)).path
                                  for fkey, farg in files_crumb_args.items()},
                                 base_directory=data_crumb.split()[0]), name="select")

    # Connect, e.g., 'infosrc.subject_id' to 'select.subject_id'
    wf.connect([(infosource, select, [(field, field) for field in field_names])])

    return wf
