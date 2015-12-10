"""
An example of how to build, plot and run an workflow.
Note the commented lines would set this as an `invoke` task function
 to allow this to be run from the command line.
"""

#from invoke import task

from pypes.run import run_wf
from pypes.plot import plot_workflow
from pypes.datasets import cobre_workflow


#@task
def run(wf_name="spm_anat_preproc", base_dir="", cache_dir="", output_dir="",
        plugin="MultiProc", n_cpus=4):
    """

    ParametersA
    ----------
    wf_name: str

    base_dir: str

    cache_dir: str

    output_dir: str

    plugin: str

    n_cpus: int
    """
    wf = cobre_workflow(wf_name=wf_name,
                        base_dir=base_dir,
                        cache_dir=cache_dir,
                        output_dir=output_dir,
                        )

    # plot the graph in its working directory
    plot_workflow(wf)

    try:
        # run it
        run_wf(wf, plugin, n_cpus)
    finally:
        import ipdb, sys
        print(sys.exc_info())
        ipdb.post_mortem(sys.exc_info()[2])
