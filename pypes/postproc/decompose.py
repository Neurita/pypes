# -*- coding: utf-8 -*-
"""
Nipype interfaces to perform signal decomposition.
"""

import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io
from   nipype.interfaces.utility import IdentityInterface, Function

from   .plotting import plot_ica_results
from   ..nilearn.canica import CanICAInterface
from   ..nilearn.utils import concat_imgs
from   ..config  import setup_node, get_config_setting
from   ..utils   import (get_trait_value,
                         get_datasink,)


def attach_canica(main_wf, wf_name="canica", **kwargs):
    """ Attach a nilearn CanICA interface to `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    kwargs: dict[str]->str
        input_node: str
            Name of the input node from where to connect the source `input_connect`.

        input_connection: str
            Name of the connection to obtain the source files.

    Nipype Inputs for `main_wf`
    ---------------------------
    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    srcwf_name   = kwargs['input_node']
    srcconn_name = kwargs['input_connection']

    src_wf   = main_wf.get_node(srcwf_name)
    datasink = get_datasink(main_wf, name='datasink')

    base_outdir  = datasink.inputs.base_directory
    ica_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,),
                           name="{}_datasink".format(wf_name))

    # the list of the subjects files
    ica_subjs = pe.JoinNode(interface=IdentityInterface(fields=["ica_subjs"]),
                            joinsource="infosrc",
                            joinfield="ica_subjs",
                            name="ica_subjs")

    # warp each subject to the group template
    canica = setup_node(CanICAInterface(), name="{}_ica".format(wf_name),)

    # Connect the nodes
    main_wf.connect([
                     # file list input
                     (src_wf,     ica_subjs, [(srcconn_name, "ica_subjs")]),

                     # canica
                     (ica_subjs,  canica,    [("ica_subjs",  "in_files")]),

                     # canica output
                     (canica, ica_datasink,  [("components", "canica.@components")]),
                     (canica, ica_datasink,  [("loadings",   "canica.@loadings")]),
                     (canica, ica_datasink,  [("score",      "canica.@score")]),
                   ])
    return main_wf


def attach_concat_canica(main_wf, wf_name="canica", **kwargs):
    """ Attach a Concat and a nilearn CanICA interface to `main_wf`.

    The Concat node will merge all the files together in one 4D volume before delivering it to CanICA.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    kwargs: dict[str]->str
        input_node: str
            Name of the input node from where to connect the source `input_connect`.

        input_connection: str
            Name of the connection to obtain the source files.

    Nipype Inputs for `main_wf`
    ---------------------------
    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    srcwf_name   = kwargs['input_node']
    srcconn_name = kwargs['input_connection']

    src_wf   = main_wf.get_node(srcwf_name)
    datasink = get_datasink(main_wf, name='datasink')

    base_outdir  = datasink.inputs.base_directory
    ica_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,),
                           name="ica_datasink".format(wf_name))
    ica_datasink.inputs.container = 'ica_{}'.format(wf_name)

    # the list of the raw pet subjects
    ica_subjs = pe.JoinNode(interface=IdentityInterface(fields=["ica_subjs"]),
                            joinsource="infosrc",
                            joinfield="ica_subjs",
                            name="ica_subjs")

    # concat images
    concat = setup_node(Function(function=concat_imgs,
                                 input_names=["in_files"],
                                 output_names=["out_file"],),
                        name="concat")

    makelist = lambda x: [x]

    # warp each subject to the group template
    canica = setup_node(CanICAInterface(), name="{}_ica".format(wf_name),)
    algorithm = get_config_setting("{}_ica.algorithm".format(wf_name),
                                   default=get_config_setting('canica.algorithm',
                                   default=''))
    if algorithm:
        canica.inputs.algorithm = algorithm

    # Connect the nodes
    main_wf.connect([
                     # file list input
                     (src_wf, ica_subjs, [(srcconn_name, "ica_subjs")]),

                     # concat images
                     (ica_subjs, concat, [("ica_subjs", "in_files")]),

                     # canica
                     (concat, canica, [(("out_file", makelist), "in_files")]),

                     # canica output
                     (canica, ica_datasink, [("components", "@components"),
                                             ("loadings",   "@loadings"),
                                             ("score",      "@score"),
                                            ]),
                   ])

    # plot the ICA results?
    do_plot = get_config_setting('canica_extra.plot', default=True)
    if not do_plot:
        return main_wf

    # get the plot threshold from the ICA node or the config file (in that order).
    plot_thr = get_config_setting('canica_extra.plot_thr', default=0)
    plot_thr = get_trait_value(canica.inputs, 'threshold', default=plot_thr)

    # concat images
    plot_ica = setup_node(Function(function=plot_ica_results,
                                   input_names=["ica_result", "application", "mask_file", "zscore", "bg_img"],
                                   output_names=["all_icc_plot", "iccs_plot"],),
                          name="plot_ica")
    plot_ica.inputs.zscore      = plot_thr
    plot_ica.inputs.mask_file   = get_trait_value(canica.inputs, 'mask')
    plot_ica.inputs.application = 'nilearn'

    # Connect the plotting nodes
    main_wf.connect([
                     # canica
                     (canica,   plot_ica,     [("components",   "ica_result")]),

                     # canica output
                     (plot_ica, ica_datasink, [("all_icc_plot", "@all_icc_plot"),
                                               ("iccs_plot",    "@iccs_plot"),
                                              ]),
                     ])

    return main_wf