# -*- coding: utf-8 -*-
"""
Plotting utilities to use nilearn from nipype
"""


def plot_all_components(components_img, **kwargs):
    """ Plot the components IC spatial maps in only one brain. """
    from nilearn.plotting import plot_prob_atlas
    from matplotlib import pyplot as plt

    fig = plt.figure(facecolor='white')
    p   = plot_prob_atlas(components_img, figure=fig, draw_cross=False, **kwargs)
    p.close()

    return fig


def plot_canica_components(components_img, **kwargs):
    """ Plot the components IC spatial maps in a grid."""
    import math
    from nilearn.image import iter_img
    from nilearn.plotting import plot_stat_map
    from matplotlib import pyplot as plt
    from matplotlib import gridspec

    n_ics     = len(list(iter_img(components_img)))
    n_rows    = math.ceil(n_ics/2)
    fig = plt.figure(figsize=(6, 3*n_rows), facecolor='black')
    gs  = gridspec.GridSpec(n_rows, 2)

    plots = []
    for i, ic_img in enumerate(iter_img(components_img)):
        ax = plt.subplot(gs[i])
        p  = plot_stat_map(ic_img, display_mode="z", title="IC %d" % i,
                           cut_coords=1, colorbar=False, figure=fig, axes=ax, **kwargs)
        plots.append(p)

    for p in plots:
        p.close()

    return fig


def plot_multi_slices(img, cut_dir="z", n_cuts=20, n_cols=4, figsize=(10, 20),
                      title="", title_fontsize=32, plot_func=None, **kwargs):
    """ Create a plot of `n_cuts` of `img` organized distributed in `n_cols`.
    Parameters
    ----------
    img: niimg-like

    cut_dir: str
        Sectional direction; possible values are "x", "y" or "z".

    n_cuts: int
        Number of cuts in the plot.

    n_cols: int
        Maximum number of image columns in the plot.

    figsize: 2-tuple of int
        (w, h) size in inches of the figure.

    title: str
        The superior title of the figure.

    title_fontsize: int
        The size of the title font.

    plot_func: function
       Function to plot each slice.
       Default: nilearn.plotting.plot_stat_map

    kwargs: keyword arguments
        Input arguments for plot_func.

    Returns
    -------
    fig: matplotlib.figure
    """
    import math

    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    import nilearn.plotting as niplot
    import nilearn.image as niimg

    # there is another version without grouper, but it is less efficient.
    def grouper(iterable, n, fillvalue=None):
        "Collect data into fixed-length chunks or blocks"
        # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
        import itertools
        args = [iter(iterable)] * n
        return itertools.zip_longest(fillvalue=fillvalue, *args)

    if plot_func is None:
        plot_func = niplot.plot_stat_map

    _img = niimg.load_img(img)

    n_rows = math.ceil(n_cuts/n_cols)
    cuts   = niplot.find_cut_slices(_img, n_cuts=n_cuts, direction=cut_dir)

    fig = plt.figure(figsize=figsize, facecolor='black')
    gs  = gridspec.GridSpec(n_rows, 1)

    if title:
        fig.suptitle(title, fontsize=title_fontsize)

    kwargs.setdefault('threshold', None)
    plots = []
    for i, cut_chunks in enumerate(grouper(cuts, n_cols)):
        ax = plt.subplot(gs[i])

        cut_chunks = [cut for cut in cut_chunks if cut is not None]

        try:
            p = plot_func(_img,
                          display_mode=cut_dir,
                          cut_coords=cut_chunks,
                          colorbar=True,
                          figure=fig,
                          axes=ax,
                          **kwargs)
        except:
            raise
        else:
            plots.append(p)

    for p in plots:
        p.close()

    return fig