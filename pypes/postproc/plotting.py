# -*- coding: utf-8 -*-
"""
Helper functions to plot results.
"""


def plot_connectivity_matrix(connectivity_matrix, label_names):
    """ Plot the connectivity matrix. """
    # Display the correlation matrix
    import numpy as np
    from matplotlib import pyplot as plt

    plt.figure(figsize=(10, 10))

    # Mask out the major diagonal
    np.fill_diagonal(connectivity_matrix, 0)
    fig = plt.imshow(connectivity_matrix, interpolation="nearest", cmap="RdBu_r",
                     vmax=0.8, vmin=-0.8)
    plt.colorbar()
    # And display the labels
    _ = plt.xticks(range(len(label_names)), label_names, rotation=90)
    _ = plt.yticks(range(len(label_names)), label_names)

    return fig


def plot_all_components(components_img):
    """ Plot the components IC spatial maps in only one brain. """
    from nilearn.plotting import plot_prob_atlas
    from matplotlib import pyplot as plt

    fig = plt.figure()
    plot_prob_atlas(components_img, title='All ICA components', figure=fig)

    return fig


def plot_canica_components(components_img):
    """ Plot the components IC spatial maps in a grid. """
    from nilearn.image import iter_img
    from nilearn.plotting import plot_stat_map
    from matplotlib import pyplot as plt

    n_ics  = len(list(iter_img(components_img)))
    n_rows = int(n_ics/2)
    fig, axes = plt.subplots(n_rows, 2)
    fig.set_size_inches(6, 3*n_rows)

    for i, cur_img in enumerate(iter_img(components_img)):
        ax = plt.subplot(n_rows, 2, i+1)
        plot_stat_map(cur_img, display_mode="z", title="IC %d" % i,
                      cut_coords=1, colorbar=False, figure=fig, axes=ax)

    return fig