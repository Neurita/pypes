# -*- coding: utf-8 -*-
"""
Workflows to grab input file structures.
"""

import os.path as op
from   copy import deepcopy

import nipype.pipeline.engine as pe
from   nipype.interfaces.io import DataSink, SelectFiles
from   nipype.interfaces.utility import IdentityInterface

from .utils import extend_trait_list, joinstrings
from .utils.piping import get_values_map_keys


def build_crumb_workflow(wfname_attacher, data_crumb, in_out_kwargs, output_dir,
                         cache_dir='', crumb_replaces=None):
    """ Returns a workflow for the give `data_crumb` with the attached workflows
    given by `attach_functions`.

    Parameters
    ----------
    wfname_attacher: dict[Str] -> function
        Dictionary with name of the workflow and its corresponding
         attach function that will be in charge of attaching workflows
         to the main input/output workflow.

    data_crumb: hansel.Crumb
        The crumb until the subject files.
        Example: Crumb('/home/hansel/cobre/raw/{subject_id}/session_1/{modality}/{image_file})
        The last 2 crumb arguments of `data_crumb` must be '{modality}/{image}',
        which indicates each of the subject/session files.
        This argument will be replaced by the corresponding image name.

    in_out_kwargs: dict with keyword arguments
        This arguments are for the in_out_crumb_wf.
        Mainly 'files_crumb_args' which will declare the values each file
        type the crumb arguments in `data_crumb` must be replaced with.
        Example:
              {'files_crumb_args': {'anat':  [('modality', 'anat_1'),
                                              ('image',    'mprage.nii.gz')],
                                    'rest':  [('modality', 'rest_1'),
                                              ('image',    'rest.nii.gz')],
                                   },
              }

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path

    crumb_replaces: dict with keyword arguments
        Keyword arguments with values for the data_crumb crumb path.
    """
    if crumb_replaces is not None:
        data_crumb = data_crumb.replace(**crumb_replaces)

    if not data_crumb.exists():
        raise IOError("Expected an existing folder for `data_crumb`, got {}.".format(data_crumb))

    if not wfname_attacher:
        raise ValueError("Expected `wfname_attacher` to have at least one function, "
                         "got {}.".format(wfname_attacher))

    if not in_out_kwargs:
        raise ValueError("Expected `in_out_kwargs` to have at least the parameters for"
                         "`files_crumb_args`, got {}.".format(in_out_kwargs))

    # check some args
    if not cache_dir:
        cache_dir = op.join(op.dirname(output_dir), "wd")

    # generate the workflow
    main_wf = input_output_crumb(
                                 work_dir=cache_dir,
                                 data_crumb=data_crumb,
                                 output_dir=output_dir,
                                 crumb_arg_values=dict(**crumb_replaces),
                                 input_wf_name='input_files',
                                 **in_out_kwargs)

    for wf_name, attach_wf in wfname_attacher.items():
        main_wf = attach_wf(main_wf=main_wf, wf_name=wf_name)

    # move the crash files folder elsewhere
    main_wf.config["execution"]["crashdump_dir"] = op.join(main_wf.base_dir,
                                                           main_wf.name, "log")

    return main_wf


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
    param_args = deepcopy(crumb_arg_values)

    # create the lists of argument names
    crumb_args = list(data_crumb.keys())
    arg_names  = list(param_args.keys())
    arg_names.extend(get_values_map_keys(files_crumb_args))

    # check their size and expand them if needed
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
    field_names = list(data_crumb.keys())

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


def input_output_crumb(work_dir, data_crumb, output_dir, crumb_arg_values, files_crumb_args,
                       input_wf_name=None, wf_name="main_workflow"):
    """ Creates a workflow with the `subject_session_file` input nodes and an empty `datasink`.
    The 'datasink' must be connected afterwards in order to work.

    Parameters
    ----------
    work_dir: str
        Path to the workflow temporary folder

    data_crumb: hansel.Crumb
        The crumb until the subject files.
        Example: Crumb('/home/hansel/data/{subject_id}/{session_id}/{modality}/{image_file})

    output_dir: str
        Path to where the datasink will leave the results.

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
                  'pet':  [('modality', 'pet'), ('image_file', ''pet_fdg.nii.gz'')],
                 }

    input_wf_name: src
        Name of the root input-output workflow

    wf_name: str
        Name of the main workflow

    Returns
    -------
    wf: Workflow
    """
    # create the root workflow
    main_wf = pe.Workflow(name=wf_name, base_dir=work_dir)

    # datasink
    datasink = pe.Node(DataSink(parameterization=False,
                                base_directory=output_dir,),
                                name="datasink")

    # input workflow
    # (work_dir, data_crumb, crumb_arg_values, files_crumb_args, wf_name="input_files"):
    input_wf = crumb_input(work_dir=work_dir,
                           data_crumb=data_crumb,
                           crumb_arg_values=crumb_arg_values,
                           files_crumb_args=files_crumb_args,
                           wf_name=input_wf_name)

    # basic file name substitutions for the datasink
    file_args = get_values_map_keys(files_crumb_args)
    undef_args = [name for name in set(data_crumb.keys()) if name not in file_args]

    substitutions = [(name, "") for name in undef_args]
    substitutions.append(("__", "_"))

    datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions,
                                                      substitutions)

    # connect the input_wf to the datasink
    joinpath = pe.Node(joinstrings(len(undef_args)), name='joinpath')

    # Connect the infosrc node to the datasink
    input_joins = [('infosrc.{}'.format(name), 'arg{}'.format(arg_no+1))
                   for arg_no, name in enumerate(undef_args)]

    main_wf.connect([
                     (input_wf, joinpath, input_joins),

                     (joinpath, datasink, [("out", "container")]),
                    ])

    return main_wf


def get_input_file_name(input_node, fname_key):
    """ Return the name of the file given by the node key `fname_key` in `input_node`.

    Parameters
    ----------
    input_node: nipype Node
        a node with a file input interface (SelectFiles, for now).

    fname_key: str
        The key that is used to access the file path using `input_node`.
        Example: 'anat'

    Returns
    -------
    filename: str
        The base input file path from the input node.
    """
    if isinstance(input_node.interface, SelectFiles):
        try:
            fname = input_node.interface._templates[fname_key]
        except:
            raise AttributeError("Could not find a SelectFiles node called 'select' in main workflow.")
        else:
            return fname
    else:
        raise NotImplementedError('`get_input_file_name` has not been implemented for nodes'
                                  ' of type {}.'.format(type(input_node.interface)))
