# -*- coding: utf-8 -*-
"""
Workflows to grab input file structures.
"""
import os.path as op
import logging as log

import nipype.pipeline.engine as pe
from   nipype.interfaces.io import DataSink
from   nipype.interfaces.utility import IdentityInterface
from   hansel.utils import joint_value_map, valuesmap_to_dict

from .crumb  import DataCrumb
from .utils  import extend_trait_list, joinstrings
from .       import configuration


def build_crumb_workflow(wfname_attacher, data_crumb, in_out_kwargs, output_dir,
                         cache_dir='', wf_name="main_workflow"):
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
        At least one crumb arguments of `data_crumb` must be open,
        this argument will be replaced by the corresponding image name.

    in_out_kwargs: dict with keyword arguments
        This arguments are for the in_out_crumb_wf.
        Mainly 'files_crumb_args' which will declare the values each file
        type the crumb arguments in `data_crumb` must be replaced with.
        Example:
              {'anat':  [('modality', 'anat_1'),
                         ('image',    'mprage.nii.gz')],
              'rest':   [('modality', 'rest_1'),
                         ('image',    'rest.nii.gz')],
              }

    cache_dir: str
        The working directory of the workflow.

    output_dir: str
        The output folder path.

    wf_name: str
        Name of the main workflow.
    """
    if not data_crumb.exists():
        raise IOError("Expected an existing folder for `data_crumb`, got {}.".format(data_crumb))

    if not data_crumb.isabs():
        raise IOError("Expected an absolute Crumb path for `data_crumb`, got {}.".format(data_crumb))

    if not wfname_attacher or wfname_attacher is None:
        raise ValueError("Expected `wfname_attacher` to have at least one function, "
                         "got {}.".format(wfname_attacher))

    #if not in_out_kwargs or in_out_kwargs is None:
    #    raise ValueError("Expected `in_out_kwargs` to have at least the name for sets of parameters for "
    #                     " `data_crumb`, got {}.".format(in_out_kwargs))

    # check some args
    if not cache_dir:
        cache_dir = op.join(op.dirname(output_dir), "wd")

    # print the configuration parameters
    log.info('Using the following configuration parameters:')
    log.info(configuration)

    # generate the workflow
    main_wf = crumb_wf(work_dir=cache_dir,
                       data_crumb=data_crumb,
                       output_dir=output_dir,
                       file_templates=in_out_kwargs,
                       wf_name=wf_name)

    for wf_name, attach_wf in wfname_attacher.items():
        main_wf = attach_wf(main_wf=main_wf, wf_name=wf_name)

    # move the crash files folder elsewhere
    main_wf.config["execution"]["crashdump_dir"] = op.join(main_wf.base_dir,
                                                           main_wf.name, "log")

    log.info('Workflow created.')

    return main_wf


def crumb_wf(work_dir, data_crumb, output_dir, file_templates,
             wf_name="main_workflow"):
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

    file_templates: Dict[str -> list of 2-tuple]
        Maps of crumb argument values to specify each file in the `data_crumb`.
        Example: {'anat': [('modality', 'anat'), ('image_file', 'anat_hc.nii.gz')],
                  'pet':  [('modality', 'pet'),  ('image_file', 'pet_fdg.nii.gz')],
                 }

    wf_name: str
        Name of the main workflow

    Returns
    -------
    wf: Workflow
    """
    # create the root workflow
    wf = pe.Workflow(name=wf_name, base_dir=work_dir)

    # datasink
    datasink = pe.Node(DataSink(parameterization=False,
                                base_directory=output_dir,),
                       name="datasink")

    # input workflow
    # (work_dir, data_crumb, crumb_arg_values, files_crumb_args, wf_name="input_files"):
    select_files = pe.Node(DataCrumb(crumb=data_crumb,
                                     templates=file_templates,
                                     raise_on_empty=False),
                           name='selectfiles')

    # basic file name substitutions for the datasink
    undef_args = select_files.interface._infields
    substitutions = [(name, "") for name in undef_args]
    substitutions.append(("__", "_"))

    datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions,
                                                      substitutions)

    # Infosource - the information source that iterates over crumb values map from the filesystem
    infosource = pe.Node(interface=IdentityInterface(fields=undef_args), name="infosrc")
    infosource.iterables = list(valuesmap_to_dict(joint_value_map(data_crumb, undef_args)).items())
    infosource.synchronize = True

    # connect the input_wf to the datasink
    joinpath = pe.Node(joinstrings(len(undef_args)), name='joinpath')

    # Connect the infosrc node to the datasink
    input_joins = [(name, 'arg{}'.format(arg_no+1))
                   for arg_no, name in enumerate(undef_args)]

    wf.connect([
                (infosource,   select_files, [(field, field) for field in undef_args]),
                (select_files, joinpath,     input_joins),
                (joinpath,     datasink,     [("out", "container")]),
               ],
              )

    return wf


