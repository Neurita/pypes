# -*- coding: utf-8 -*-
"""
Helper functions to plot nipype workflows
"""


def plot_workflow(wf):
    """ Plot `wf` nodes in 3 different ways in its working directory.

    Parameters
    ----------
    wf: nipype Workflow
    """
    # print the graph of the workflow in the working directory
    wf.write_graph("{}_colored_workflow".format(wf.name), graph2use="colored")
    wf.write_graph("{}_exec_workflow".format(wf.name),    graph2use="exec")