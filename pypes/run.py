"""

"""

import os.path as op

import nipype.pipeline.engine   as pe
from   nipype.interfaces.io     import DataSink

from   .input_files import subject_session_input
from   .anat        import attach_t1_preprocessing
from   .pet         import attach_pet_preprocessing
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


def run(wf_name="spm_t1_preproc", base_dir="", cache_dir="", output_dir="",
        year="", plugin="MultiProc", n_cpus=4):
    """

    Parameters
    ----------
    wf_name: str

    base_dir: str

    cache_dir: str

    output_dir: str

    year: str or int

    plugin: str

    n_cpus: int
    """
    if not year:
        data_dir = base_dir
    else:
        data_dir = op.join(base_dir, year)

    if not data_dir or not op.exists(data_dir):
        raise IOError("Expected an existing folder for `data_dir`, got {}.".format(data_dir))

    wfs = {"spm_t1_preproc": attach_t1_preprocessing,
           "spm_pet_preproc": attach_pet_preprocessing,
          }

    if wf_name not in wfs:
        raise ValueError("Expected `wf_name` to be in {}, got {}.".format(list(wfs.keys()),
                                                                          wf_name))

    # check some args
    if not output_dir:
        output_dir = op.join(op.dirname(data_dir), "out", year)

    if not cache_dir:
        cache_dir = op.join(op.dirname(data_dir), "wd", year)

    # generate the workflow
    main_wf = in_out_workflow(work_dir=cache_dir,
                              data_dir=data_dir,
                              output_dir=output_dir,
                              session_names=['session_0'],
                              file_names=['anat_hc.nii.gz', 'pet_fdg.nii.gz'],
                              subject_ids=None,
                              input_wf_name='input_files')

    wf = wfs[wf_name](main_wf=main_wf,
                      data_dir=data_dir,
                      work_dir=cache_dir,
                      output_dir=output_dir,)

    # move the crash files folder elsewhere
    wf.config["execution"]["crashdump_dir"] = op.join(wf.base_dir, wf.name, "log")

    # print the graph of the workflow in the working directory
    wf.write_graph("{}_colored_workflow".format(wf_name),  graph2use="colored")
    wf.write_graph("{}_exec_workflow".format(wf_name),     graph2use="exec")

    # run the workflow according to `plugin`
    if plugin == "MultiProc" and n_cpus > 1:
        wf.run("MultiProc", plugin_args={"n_procs": n_cpus})
    elif not plugin or plugin is None or n_cpus <= 1:
        wf.run()
    else:
        wf.run(plugin=plugin)

