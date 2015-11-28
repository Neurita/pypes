"""
Helper functions to build base workflow and run them
"""

import nipype.pipeline.engine   as pe
from   nipype.interfaces.io     import DataSink

from   .input_files import subject_session_input
from   .utils       import extend_trait_list, joinpaths


def in_out_workflow(work_dir, data_dir, output_dir, session_names, file_names, subject_ids,
                    input_wf_name, wf_name="main_workflow"):
    """ Creates a workflow with the `subject_session_file` input nodes and an empty `datasink`.
    The 'datasink' must be connected in order to work.

    Parameters
    ----------
    work_dir: str
        Path to the workflow temporary folder

    data_dir: str
        Path to where the subject folders is.

    output_dir: str
        Path to where the datasink will leave the results.

    session_names: list of str
        Example: ['session_0']

    file_names: list of str
        Example: ['mprage.nii.gz', 'rest.nii.gz']
        Example: ['anat_1/mprage.nii.gz', 'rest_1/rest.nii.gz']

    subject_ids: list of str
        Use this if you want to limit the analysis to certain subject IDs.
        If `None` will pick the folders from os.listdir(data_dir).

    input_wf_name: src

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
    input_wf = subject_session_input(base_dir=data_dir,
                                     session_names=session_names,
                                     file_names=file_names,
                                     subject_ids=subject_ids,
                                     wf_name=input_wf_name)

    joinpath = pe.Node(joinpaths(), name='joinpath')

    # basic file name substitutions for the datasink
    substitutions = [("_subject_id", ""),
                     ("_session_id_", ""),
                    ]

    datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions,
                                                      substitutions)

    # Connect the infosrc node to the datasink
    main_wf.connect([
                      # datasink
                      (input_wf, joinpath, [("infosrc.subject_id", "arg1"),
                                            ("infosrc.session_id", "arg2")]),
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

