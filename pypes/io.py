# -*- coding: utf-8 -*-
"""
Workflows to grab input file structures.
"""

import os.path as op

import nipype.pipeline.engine as pe
from   nipype.interfaces.io import DataSink, SelectFiles

from .crumb import crumb_input
from .utils import extend_trait_list, joinstrings
from .utils.piping import get_values_map_keys


def build_crumb_workflow(wfname_attacher, data_crumb, in_out_kwargs, output_dir,
                         cache_dir='', crumb_replaces=None, params=None):
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

    if not wfname_attacher or wfname_attacher is None:
        raise ValueError("Expected `wfname_attacher` to have at least one function, "
                         "got {}.".format(wfname_attacher))

    if not in_out_kwargs or in_out_kwargs is None:
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
                                 crumb_arg_values=dict(**crumb_replaces) if crumb_replaces is not None else None,
                                 input_wf_name='input_files',
                                 files_crumb_args=in_out_kwargs['files_crumb_args'])

    for wf_name, attach_wf in wfname_attacher.items():
        main_wf = attach_wf(main_wf=main_wf, wf_name=wf_name, params=params)

    # move the crash files folder elsewhere
    main_wf.config["execution"]["crashdump_dir"] = op.join(main_wf.base_dir,
                                                           main_wf.name, "log")

    return main_wf


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
    undef_args = [name for name in list(data_crumb.open_args()) if name not in file_args]

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
