# -*- coding: utf-8 -*-
"""
Helper functions to plot nipype workflows
"""


def plot_workflow(wf, detailed=False):
    """ Plot the `wf` pipeline nodes in different ways in its working directory.

    Parameters
    ----------
    wf: nipype Workflow

    detailed: bool
        If True will also plot the detailed execution plot.
        This takes a very long time if the number of subjects if large.
    """
    # print the graph of the workflow in the working directory
    try:
        wf.write_graph("{}_colored_workflow".format(wf.name),
                       graph2use="colored")
    except RuntimeError as re:
        print('Error plotting colored workflow: {}.'.format(str(re)))
    else:
        if detailed:
            wf.write_graph("{}_exec_workflow".format(wf.name),
                           graph2use="exec")
