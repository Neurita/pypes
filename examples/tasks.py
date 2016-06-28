"""
An example of how to build, plot and run an workflow.
Note the commented lines would set this as an `invoke` task function
 to allow this to be run from the command line.
"""
import os.path as op

from invoke import task
from hansel import Crumb

from pypes.run import run_wf
from pypes.plot import plot_workflow


def run_pype(workflow, plugin="MultiProc", n_cpus=4):
    try:
        # plot the graph in its working directory
        plot_workflow(workflow)

        # run it
        run_wf(workflow, plugin, n_cpus)
    except:
        raise
    finally:
        import pdb, sys
        print(sys.exc_info())
        pdb.post_mortem(sys.exc_info()[2])


@task
def michael_pype(wf_name="camino_dti_tract_and_pet", base_dir="", cache_dir="", output_dir="",
                 plugin="MultiProc", n_cpus=4):
    """ Run the

    ParametersA
    ----------
    wf_name: str

    base_dir: str

    cache_dir: str

    output_dir: str

    year: str or int

    plugin: str

    n_cpus: int
    """
    from pypes.datasets import clinical_crumb_workflow

    data_path = op.join(op.expanduser(base_dir), '{diagnosis}', '{subject_id}', '{session_id}', '{image}')
    data_crumb = Crumb(data_path, ignore_list=['.*'])

    atlas_file = HAMM_MNI

    wf = clinical_crumb_workflow(wf_name    = wf_name,
                                 data_crumb = data_crumb,
                                 cache_dir  = op.abspath(op.expanduser(cache_dir)) if cache_dir else '',
                                 output_dir = op.abspath(op.expanduser(output_dir)) if output_dir else '',
                                 params={'atlas_file': atlas_file},
                                 )
    run_pype(wf, plugin=plugin, n_cpus=n_cpus)


@task
def tum_pype(wf_name="spm_anat_preproc", base_dir="", cache_dir="", output_dir="", plugin="MultiProc", n_cpus=4):
    """ Run the

    ParametersA
    ----------
    wf_name: str

    base_dir: str

    cache_dir: str

    output_dir: str

    year: str or int

    plugin: str

    n_cpus: int
    """
    from pypes.datasets import clinical_crumb_workflow

    data_path = op.join(op.expanduser(base_dir), '{year}', '{subject_id}', '{session_id}', '{image}')
    data_crumb = Crumb(data_path, ignore_list=['.*'])

    atlas_file = HAMM_MNI

    wf = clinical_crumb_workflow(wf_name    = wf_name,
                                 data_crumb = data_crumb,
                                 cache_dir  = op.abspath(op.expanduser(cache_dir)) if cache_dir else '',
                                 output_dir = op.abspath(op.expanduser(output_dir)) if output_dir else '',
                                 params={'atlas_file': atlas_file}
                                 )
    run_pype(wf, plugin=plugin, n_cpus=n_cpus)


@task
def cobre_pype(wf_name="spm_anat_preproc", base_dir="", cache_dir="", output_dir="", plugin="MultiProc", n_cpus=4):
    """ Run the

    ParametersA
    ----------
    wf_name: str

    base_dir: str
        Base path to where the data is

    cache_dir: str

    output_dir: str

    year: str or int

    plugin: str

    n_cpus: int
    """
    from pypes.datasets import cobre_crumb_workflow

    data_path = op.join(op.expanduser(base_dir), '{subject_id}', 'session_1', '{modality}', '{image}')
    data_crumb = Crumb(data_path, ignore_list=['.*'])

    wf = cobre_crumb_workflow(wf_name    = wf_name,
                              data_crumb = data_crumb,
                              cache_dir  = op.abspath(op.expanduser(cache_dir)) if cache_dir else '',
                              output_dir = op.abspath(op.expanduser(output_dir)) if output_dir else '',
                             )
    run_pype(wf, plugin=plugin, n_cpus=n_cpus)
