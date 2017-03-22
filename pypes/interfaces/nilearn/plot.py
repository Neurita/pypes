# -*- coding: utf-8 -*-
"""
Plotting utilities to use nilearn from nipype
"""
import logging


def plot_all_components(components_img, **kwargs):
    """ Plot the components IC spatial maps in only one brain. """
    from nilearn.plotting import plot_prob_atlas
    from matplotlib import pyplot as plt

    fig = plt.figure(facecolor='white')
    p   = plot_prob_atlas(components_img, figure=fig, draw_cross=False, **kwargs)
    p.close()

    return fig


def plot_ica_components(components_img, **kwargs):
    """ Plot the components IC spatial maps in a grid."""
    import math
    from nilearn.image import iter_img
    from nilearn.plotting import plot_stat_map
    from matplotlib import pyplot as plt
    from matplotlib import gridspec

    n_ics  = len(list(iter_img(components_img)))
    n_rows = math.ceil(n_ics/2)
    fig = plt.figure(figsize=(6, 3*n_rows), facecolor='black')
    gs  = gridspec.GridSpec(n_rows, 2)

    plots = []
    for i, ic_img in enumerate(iter_img(components_img)):
        ax = plt.subplot(gs[i])
        p  = plot_stat_map(ic_img, display_mode="z", title="IC {}".format(i+1),
                           cut_coords=1, colorbar=False, figure=fig, axes=ax, **kwargs)
        plots.append(p)

    for p in plots:
        p.close()

    return fig


def plot_multi_slices(img, cut_dir="z", n_cuts=20, n_cols=4, figsize=(2.5, 3),
                      title="", title_fontsize=32, plot_func=None, black_bg=True, **kwargs):
    """ Create a plot of `n_cuts` of `img` organized distributed in `n_cols`.
    Parameters
    ----------
    img: niimg-like

    cut_dir: str, optional
        Sectional direction; possible values are "x", "y" or "z".

    n_cuts: int, optional
        Number of cuts in the plot.

    n_cols: int, optional
        Maximum number of image columns in the plot.

    figsize: 2-tuple of int, optional
        (w, h) size in inches of the figure.

    title: str, optional
        The superior title of the figure.

    title_fontsize: int, optional
        The size of the title font.

    plot_func: function, optional
       Function to plot each slice.
       Default: nilearn.plotting.plot_stat_map

    black_bg: boolean, optional
        If True, the background of the image is set to be black.
        If you wish to save figures with a black background, you will need to pass
        "facecolor='k', edgecolor='k'" to matplotlib.pyplot.savefig.

    kwargs: keyword arguments, optional
        Input arguments for plot_func.

    Returns
    -------
    fig: matplotlib.figure
    """
    import math

    import matplotlib.pyplot as plt
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

    n_rows = 1
    if n_cuts > n_cols:
        n_rows = math.ceil(n_cuts/n_cols)

    spacing = kwargs.get('spacing', 'auto')
    cuts = niplot.find_cut_slices(_img, n_cuts=n_cuts,
                                  direction=cut_dir,
                                  spacing=spacing)

    # instantiate the figure
    if black_bg:
        facecolor = 'black'
        titlecolor = 'white'
    else:
        facecolor = 'white'
        titlecolor = 'black'

    figsize = figsize[0] * n_cols, figsize[1] * n_rows
    fig = plt.figure(figsize=figsize, facecolor=facecolor)
    gs  = gridspec.GridSpec(n_rows, 1)

    if title:
        fig.suptitle(title, fontsize=title_fontsize,
                     color=titlecolor)

    # for the superior plot
    put_colorbar = True

    # make the plots
    for i, cut_chunks in enumerate(grouper(cuts, n_cols)):
        ax = plt.subplot(gs[i])

        cut_chunks = [cut for cut in cut_chunks if cut is not None]

        try:
            p = plot_func(_img,
                          display_mode=cut_dir,
                          cut_coords=cut_chunks,
                          colorbar=put_colorbar,
                          figure=fig,
                          axes=ax,
                          black_bg=black_bg,
                          **kwargs)
        except IndexError:
            logging.warning('Could not plot for coords {}.'.format(cut_chunks))
        finally:
            put_colorbar = False

    return fig


def plot_ortho_slices(img, n_cuts=4, n_cols=6, figsize=(2.5, 3),
                      title="", title_fontsize=32, plot_func=None,
                      black_bg=True, **kwargs):
    """
    Parameters
    ----------
    img: niimg-like

    n_cuts: int, optional
        Number of cuts for each dimension (X, Y, and Z) in the plot.

    n_cols: int, optional
        Maximum number of image columns in the plot.

    figsize: 2-tuple of int, optional
        (w, h) size in inches of one slice image.

    title: str, optional
        The superior title of the figure.

    title_fontsize: int, optional
        The size of the title font.

    plot_func: function, optional
       Function to plot each slice.
       Default: nilearn.plotting.plot_stat_map

    black_bg: boolean, optional
        If True, the background of the image is set to be black.
        If you wish to save figures with a black background, you will need to pass
        "facecolor='k', edgecolor='k'" to matplotlib.pyplot.savefig.

    kwargs: keyword arguments, optional
        Input arguments for plot_func.

    Returns
    -------
    fig: matplotlib.figure
    """
    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    import nilearn.plotting as niplot
    import nilearn.image as niimg

    if plot_func is None:
        plot_func = niplot.plot_stat_map

    # load the image file
    _img = niimg.load_img(img)

    directions = ('x', 'y', 'z')

    # calculate the shape of the figure
    total_cuts = n_cuts * 3
    if total_cuts > n_cols:
        n_rows = 3
        n_cols = 1
        colorbard_idx = 0
    else:
        n_rows = 1
        n_cols = 3
        colorbard_idx = len(directions)-1

    # calculate the cut coordinates for each direction
    cuts = []
    spacing = kwargs.get('spacing', 'auto')
    for cut_dir in directions:
        dir_cuts = niplot.find_cut_slices(_img, n_cuts=n_cuts,
                                          direction=cut_dir,
                                          spacing=spacing)
        cuts.append((cut_dir, dir_cuts))

    # instantiate the figure
    if black_bg:
        facecolor = 'black'
        titlecolor = 'white'
    else:
        facecolor = 'white'
        titlecolor = 'black'

    figsize = figsize[0] * n_cols * n_cuts, figsize[1] * n_rows
    fig = plt.figure(figsize=figsize, facecolor=facecolor)
    gs  = gridspec.GridSpec(n_rows, n_cols)

    # put the title, if any
    if title:
        fig.suptitle(title, fontsize=title_fontsize,
                     color=titlecolor)

    # plot on the figure
    for i, (cut_dir, cut_coords) in enumerate(cuts):
        ax = plt.subplot(gs[i])

        put_colorbar = True if i == colorbard_idx else False
        try:
            p = plot_func(_img,
                          display_mode=cut_dir,
                          cut_coords=cut_coords,
                          colorbar=put_colorbar,
                          figure=fig,
                          axes=ax,
                          black_bg=black_bg,
                          **kwargs)
        except IndexError:
            logging.warning('Could not plot for coords {}.'.format(cut_coords))

    return fig


def plot_stat_overlay(stat_img, contour_img, bg_img, **kwargs):
    """Plot over bg_img a stat_img and the countour."""
    import nilearn.plotting as niplot

    if bg_img is not None:
        kwargs['bg_img'] = bg_img

    display = niplot.plot_stat_map(stat_img, **kwargs)
    display.add_contours(contour_img, filled=True, alpha=0.6, levels=[0.5], colors='g')
    return display


def plot_overlays(stat_imgs, contour_imgs, bg_img=None,
                  figsize=(2.5, 3), title_fontsize=32, title='', **kwargs):
    """Plots each contour_imgs as an overlay of its corresponding `stat_imgs`.
    `contour_imgs` and `stat_imgs` must have the same length."""

    from matplotlib import pyplot as plt
    from matplotlib import gridspec

    _stat_imgs = list(stat_imgs)
    _cnts_imgs = list(contour_imgs)
    n_stats    = len(_stat_imgs)
    n_conts    = len(_cnts_imgs)
    if n_stats != n_conts:
        raise AttributeError('The length of `stat_imgs` and `contour_imgs` are '
                             'different, got {} and {}.'.format(n_stats, n_conts))

    n_rows = n_conts
    n_cols = 3 # because I am doing the ortho plot

    figsize = figsize[0] * n_cols, figsize[1] * n_rows
    fig = plt.figure(figsize=figsize, facecolor='black')
    gs  = gridspec.GridSpec(n_rows, 1)

    # put the title, if any
    if title:
        fig.suptitle(title, fontsize=title_fontsize, color='white')

    for i, (simg, cimg) in enumerate(zip(_stat_imgs, _cnts_imgs)):
        ax = plt.subplot(gs[i])
        p = plot_stat_overlay(simg, cimg, bg_img=bg_img, figure=fig, axes=ax, **kwargs)

    return fig
