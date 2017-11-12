# -*- coding: utf-8 -*-
"""
Helper functions to build base workflow and run them
"""
from neuro_pypes.plot import plot_workflow


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
    if plugin == "MultiProc" or n_cpus > 1:
        wf.run("MultiProc", plugin_args={"n_procs": n_cpus})
    elif not plugin or plugin is None or n_cpus <= 1:
        wf.run(plugin=None)
    else:
        wf.run(plugin=plugin, **plugin_kwargs)


def run_debug(workflow, plugin="MultiProc", n_cpus=4, **plugin_kwargs):
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

    try:
        # plot the graph in its working directory
        plot_workflow(workflow)

        # run it
        run_wf(workflow, plugin=plugin, n_cpus=n_cpus,
               **plugin_kwargs)
    except:
        import sys
        print(sys.exc_info())

        import pdb
        pdb.post_mortem(sys.exc_info()[2])

        raise
    else:
        print('Workflow successfully finished.')
