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
    plt.imshow(connectivity_matrix, interpolation="nearest", cmap="RdBu_r",
               vmax=0.8, vmin=-0.8)
    plt.colorbar()
    # And display the labels
    _ = plt.xticks(range(len(label_names)), label_names, rotation=90)
    _ = plt.yticks(range(len(label_names)), label_names)


def plot_all_components(components_img):
    """ Plot the components IC spatial maps in only one brain. """
    from nilearn.plotting import plot_prob_atlas, show

    p = plot_prob_atlas(components_img, title='All ICA components')
    show()
    p.close()


def plot_canica_components(components_img):
    """ Plot the components IC spatial maps in a grid. """
    from nilearn.image import iter_img
    from nilearn.plotting import plot_stat_map, show

    for i, cur_img in enumerate(iter_img(components_img)):
        plot_stat_map(cur_img, display_mode="z", title="IC %d" % i,
                      cut_coords=1, colorbar=False)

    show()
    p.close()