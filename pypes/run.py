# -*- coding: utf-8 -*-
"""
Helper functions to build base workflow and run them
"""

import nipype.pipeline.engine   as pe
from   nipype.interfaces.io     import DataSink

from   .io    import subject_session_input, crumb_input_wf
from   .utils import extend_trait_list, joinstrings
from   .utils.piping import get_values_map_keys


def in_out_crumb_wf(work_dir, data_crumb, output_dir, crumb_arg_values, files_crumb_args,
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
    input_wf = crumb_input_wf(work_dir=work_dir,
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


def run_wf(wf, plugin='MultiProc', n_cpus=2, **plugin_kwargs):
    """ Execute `wf` with `plugin`.

    Parameters
    ----------
    wf: nipype Workflow

    plugin: str
        The pipeline execution plugin.
        See wf.run docstring for choices.

    n_cpus: int
        Number of CPUs to use with the 'MultiProc' plugin.

    plugin_kwargs: keyword argumens
        Keyword arguments for the plugin if using something different
        then 'MultiProc'.
    """
    # run the workflow according to `plugin`
    if plugin == "MultiProc" and n_cpus > 1:
        wf.run("MultiProc", plugin_args={"n_procs": n_cpus})
    elif not plugin or plugin is None or n_cpus <= 1:
        wf.run()
    else:
        wf.run(plugin=plugin, **plugin_kwargs)

