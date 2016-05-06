# -*- coding: utf-8 -*-
"""
Plotting utilities to use nilearn from nipype
"""


def plot_all_components(components_img, **kwargs):
    """ Plot the components IC spatial maps in only one brain. """
    from nilearn.plotting import plot_prob_atlas
    from matplotlib import pyplot as plt

    fig = plt.figure()
    p = plot_prob_atlas(components_img, title='All ICA components', figure=fig, **kwargs)
    p.close()

    return fig


def plot_canica_components(components_img, **kwargs):
    """ Plot the components IC spatial maps in a grid. """
    from nilearn.image import iter_img
    from nilearn.plotting import plot_stat_map
    from matplotlib import pyplot as plt

    n_ics  = len(list(iter_img(components_img)))
    n_rows = int(n_ics/2)
    fig, axes = plt.subplots(n_rows, 2)
    fig.set_size_inches(6, 3*n_rows)

    plts = []
    for i, cur_img in enumerate(iter_img(components_img)):
        ax = plt.subplot(n_rows, 2, i+1)
        p = plot_stat_map(cur_img, display_mode="z", title="IC %d" % i,
                          cut_coords=1, colorbar=False, figure=fig, axes=ax, **kwargs)
        plts.append(p)

    for p in plts:
        p.close()

    return fig