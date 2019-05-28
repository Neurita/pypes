#!python
import os
import pathlib
from functools import partial

import click
import pandas as pd

import hansel
from neuro_pypes.cli.plot_helpers import create_imglist_html
from neuro_pypes.cli.utils import (
    CONTEXT_SETTINGS,
    CrumbPath,
    ExistingDirPath,
    ExistingFilePath,
    UnexistingFilePath,
    check_not_none,
    _get_plot_file_pairs,
    Spreadsheet
)


# declare the CLI group
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '--bg',
    '--background',
    type=CrumbPath(),
    callback=check_not_none,
    required=True,
    help='The hansel.Crumb path to the background images.'
)
@click.option(
    '--fg',
    '--foreground',
    type=CrumbPath(),
    required=False,
    help='The hansel.Crumb path to the foreground images.'
)
@click.option(
    '-s',
    '--style',
    type=click.Choice(['roi', 'anat', 'quant', 'outline']),
    required=False,
    default='anat',
    help='The style of the plot:'
         'roi: discrete colored ROI overlay,\n'
         'anat: gray-scale plotting,\n'
         'quant: statistical overlay,\n'
         'outline: outlined foreground.'
)
@click.option(
    '-n',
    '--ncuts',
    type=int,
    default=4,
    help='Number of cuts.'
)
@click.option(
    '-c',
    '--ncols',
    type=int,
    default=4,
    help='Number of columns.'
)
@click.option(
    '-d',
    '--cut_dir',
    type=click.Choice(['x', 'y', 'z', 'o']),
    default='o',
    help='Direction of the slice cut. "o" is for the 3 orthogonal cuts.'
)
@click.option(
    '-a',
    '--alpha',
    type=float,
    default=0.5,
    help='Transparency alpha value for the foreground.'
)
@click.option(
    '-o',
    '--out_file',
    type=UnexistingFilePath,
    default='slices.html',
    help='The output file path, it will also create a folder with '
         'the same name without the extension to store the plot images.'
)
def plot(
    background: hansel.Crumb,
    foreground: hansel.Crumb,
    style: str,
    ncuts: int,
    ncols: int,
    cut_dir: str,
    alpha: float,
    out_file: click.Path
):
    """Plot slices of the `bg` (background) images, overlaid with the `fg` (foreground)
    images if provided.

    Examples: \n
    nitap plot --bg "/data/hansel/cobre/{sid}/{session}/anat.nii.gz" --fg "/data/hansel/cobre/{sid}/{session}/mni_in_anat_space.nii.gz"\n
    nitap plot --bg "/data/hansel/cobre/{sid}/{session}/anat.nii.gz"\n
    """
    import nilearn.plotting as niplot

    from neuro_pypes.interfaces.nilearn import plot_multi_slices, plot_ortho_slices

    plotting_styles = {
        'roi': niplot.plot_roi,
        'anat': niplot.plot_anat,
        'quant': niplot.plot_stat_map,
        'outline': niplot.plot_stat_map,
    }
    plotting_style = plotting_styles[style]

    plotting_funcs = {
        'x': partial(plot_multi_slices, cut_dir='x'),
        'y': partial(plot_multi_slices, cut_dir='y'),
        'z': partial(plot_multi_slices, cut_dir='z'),
        'o': plot_ortho_slices
    }
    plotting_func = partial(plotting_funcs[cut_dir], plot_func=plotting_style, n_cuts=ncuts, n_cols=ncols, alpha=alpha)

    index_file = pathlib.Path(out_file)
    output_plotdir = pathlib.Path(index_file.parent) / index_file.stem
    output_plotdir.mkdir(exist_ok=True)

    plot_files = []

    file_pairs = _get_plot_file_pairs(background, foreground)

    for bg_file, fg_file in file_pairs:
        out_plot_file = bg_file.replace(os.path.sep, '_')
        if fg_file is not None:
            out_plot_file += '-' + fg_file.replace(os.path.sep, '_')
            fig = plotting_func(str(fg_file), bg_img=str(bg_file))
        else:
            fig = plotting_func(str(bg_file))

        plot_file = output_plotdir / (out_plot_file + '.png')
        fig.savefig(str(plot_file), facecolor='k', edgecolor='k', bbox_inches='tight')
        plot_files.append(plot_file)
        click.echo('Created file://{}.'.format(plot_file))

    create_imglist_html(plot_files, output_filepath=index_file)
    click.echo('Created index file in file://{}.'.format(index_file))


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '-i',
    '--input',
    type=CrumbPath(),
    callback=check_not_none,
    required=True,
    help='The hansel.Crumb path to the background images.'
)
@click.option(
    '-e',
    '--extra',
    type=Spreadsheet(),
    required=False,
    help='A spreadsheet (CSV file) with extra data for the subjects in the input. \n'
         'The first row, with column names, must be match at least one argument.'
)
@click.option(
    '-o',
    '--out_file',
    type=UnexistingFilePath,
    default='motion_stats.xls',
    help='The output Excel spreadsheet with the subjects motion stats.'
)
def motion(input: hansel.Crumb, extra: pd.DataFrame, out_file: hansel.Crumb):
    """ Create in `out_path` an Excel spreadsheet with some of the motion statistics obtained from the
    `statistics_files` output of the nipype.RapidArt found in the hansel.Crumb `motion_file_cr`.

    Examples: \n
    nitap motion -i "/data/hansel/cobre/{sid}/{session}/rest/artifact_stats/motion_stats.json" -o motion.xls\n
    nitap motion -i "/data/nuk/out/{group}/{sid}/session_0/rest/artifact_stats/motion_stats.json" -o motion.xls\n
    """
    from neuro_pypes.fmri.utils import motion_stats_sheet

    crumb_args = list(input.open_args())
    df = motion_stats_sheet(input, crumb_args)

    if extra:
        extra_columns = set(extra.columns.values)
        matched_args = extra_columns.intersection(crumb_args)
        if not matched_args:
            click.fail('Found no matches in the spreadsheet file between: '
                       '"{}" and "{}".'.format(extra_columns, crumb_args))
        df.join(extra, on=matched_args)

    df.to_excel(out_file)
    print('Successfully wrote the motions spreadsheet in "{}".'.format(out_file))


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '-i',
    '--ica_results_dir',
    type=ExistingDirPath,
    callback=check_not_none,
    required=True,
    help='Path to the ICA analysis results folder.'
)
@click.option(
    '-a',
    '--app',
    type=click.Choice(('nilearn', 'sbm', 'gift', 'gift-group')),
    default='gift',
    required=False,
    help='An identification of the tool used for the ICA analysis.'
)
@click.option(
    '-l',
    '--labels',
    type=ExistingFilePath,
    required=False,
    help='A CSV file with two columns: "subject_id" and "group".'
    'The subject_ids must be in the paths contained in the Subject.mat'
    'file and match the `subjid_pat` argument.'
)
@click.option(
    '-m',
    '--mask',
    type=ExistingFilePath,
    required=False,
    help='Path to a mask file to filter out non-brain areas from the IC spatial map.'
)
@click.option(
    '-z',
    '--zscore',
    type=float,
    default=2.,
    required=False,
    help='Upper bound threshold for the IC spatial maps.'
)
@click.option(
    '-s',
    '--subjid_pat',
    type=str,
    default=r'(?P<patid>[a-z]{2}_[0-9]{6})',
    required=False,
    help='A search regex pattern that returns one group element that contains the subject id.'
         'This will be used to search for subject_id in the file paths contained in the Subjects.mat file.'
)
@click.option(
    '-o',
    '--out_file',
    type=UnexistingFilePath,
    required=False,
    default='loading_coefficients.xls',
    help='The output Excel spreadsheet with the loading coefficients.'
)
def loadings(ica_results_dir, app, labels, mask, background, zscore, subjid_pat, out_file):
    """
    Save the Excel loadings files in the `ica_out_dir`.
    The `out_file` has the loadings as is, with the subjects IDs and group.
    The other file has a `weighted_` preffix which has the loading signs changed according to
    the average correlation value of the "main" region of each of the IC spatial maps.

    Examples: \n
    nitap loadings -i "/data/hansel/cobre/rest_ica/experiment1" -o subject_loadings.xls\n
    nitap loadings -i "/data/hansel/cobre/rest_ica/experiment2" -o subject_loadings.xls\n
    """
    from pypes.networks.plotting import _pick_plotter

    check_blob_plot = 'check_sign_blobs.png'
    weighted_loadings_out_file = 'weighted_{}'.format(out_file)

    plotter_class = _pick_plotter(app)
    click.echo('Using {} to read the results data.'.format(plotter_class.__name__))

    plotter = plotter_class(ica_results_dir)
    plotter.fit(mask_file=mask, mode='+-', zscore=zscore)

    sdf = plotter.simple_loadings_df(group_labels_file=labels, subjid_pat=subjid_pat)
    sdf.to_excel(out_file)
    click.echo('Created {} with loading coefficients.'.format(os.path.abspath(out_file)))

    pdf = plotter.weighted_loadings_df(group_labels_file=labels, subjid_pat=subjid_pat)
    pdf.to_excel(weighted_loadings_out_file)
    click.echo('Created {} with weighted loading coefficients.'.format(os.path.abspath(weighted_loadings_out_file)))

    check_blob_plot = os.path.join(os.path.dirname(out_file), check_blob_plot)
    plotter.plot_icmaps_and_blobs(check_blob_plot, bg_img=background)
    click.echo(
        'Created {} with coloured blob plots to verify weighted loading coefficients.'.format(
            os.path.abspath(check_blob_plot)
        )
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '-i',
    '--ica_results_dir',
    type=ExistingDirPath,
    callback=check_not_none,
    required=True,
    help='Path to the ICA analysis results folder.'
)
@click.option(
    '-a',
    '--app',
    type=click.Choice(('nilearn', 'sbm', 'gift', 'gift-group')),
    default='gift',
    required=False,
    help='An identification of the tool used for the ICA analysis.'
)
@click.option(
    '-m',
    '--mask',
    type=ExistingFilePath,
    required=False,
    help='Path to a mask file to filter out non-brain areas from the IC spatial map.'
)
@click.option(
    '-b',
    '--background',
    type=ExistingFilePath,
    required=False,
    help='A background image for the blob-check report, to verify that the blobs'
         'taken into account for the loadings signs are correct.'
)
@click.option(
    '-z',
    '--zscore',
    type=float,
    default=2.,
    required=False,
    help='Upper bound threshold for the IC spatial maps.'
)
@click.option(
    '-o',
    '--out_dir',
    type=ExistingDirPath,
    required=False,
    help='The output path for the plot files. If None will use `ica_results_dir`.'
)
def spatial_maps(ica_results_dir, app, mask, background,  zscore, out_dir):
    """
    Plot the spatial maps in `ica_results_dir`.

    Examples: \n
    nitap spatial_maps -i "/data/hansel/cobre/rest_ica/experiment1" -a canica -z 2 -b /data/hansel/std_brains/templates/mni.nii\n
    nitap spatial_maps -i "/data/hansel/cobre/rest_ica/experiment2" -a gift -z 1\n
    """
    from neuro_pypes.ica.plotting import _pick_plotter

    plotter_class = _pick_plotter(app)
    click.echo('Using {} to read the results data.'.format(plotter_class.__name__))

    plotter = plotter_class(ica_results_dir)
    plotter.fit(mask_file=mask, mode='+-', zscore=zscore)

    if not out_dir:
        out_dir = ica_results_dir

    plotter.plot_icmaps(bg_img=background)
    click.echo('Created IC spatial map plots in {}.'.format(out_dir))


# @task
# def dcm2nii(input_crumb_path, output_dir, regex='fnmatch', ncpus=3):
#     """ Convert all DICOM files within `input_crumb_path` into NifTI in `output_folder`.
#
#     Will copy only the NifTI files reoriented by MRICron's dcm2nii command.
#     Will rename the NifTI files that are matched with recognized modalities to the short
#     modality name from config.ACQ_PATTERNS.
#
#     Parameters
#     ----------
#     input_dir: str
#         A crumb path str indicating the whole path until the DICOM files.
#         Example: '/home/hansel/data/{group}/{subj_id}/{session_id}/{acquisition}/{dcm_file}
#
#         The crumb argument just before the last one will be used as folder container reference
#         for the DICOM series.
#
#     output_dir: str
#         The root folder path where to save the tree of nifti files.
#         Example: '/home/hansel/nifti'
#         This function will create the same tree as the crumbs in input_crumb_path, hence
#         for the example above the output would have the following structure:
#         '/home/hansel/nifti/{group}/{subj_id}/{session_id}/{nifti_file}'
#
#         Where {nifti_file} will take the name from the {acquisition} or from the
#         patterns in ACQ_PATTERNS in `config.py` file.
#
#     regex: str
#         The regular expression syntax you may want to set in the Crumbs.
#         See hansel.Crumb documentation for this.
#
#     ncpus: int
#         If PARALLEL is set to True in `config.py` this says the number of
#         processes that will be launched for dcm2nii in parallel.
#     """
#     from boyle.dicom.convert import convert_dcm2nii
#
#     input_dir  = op.expanduser(input_crumb_path)
#     output_dir = op.expanduser(output_dir)
#
#     if not op.exists(output_dir):
#         log.info('Creating output folder {}.'.format(output_dir))
#         os.makedirs(output_dir)
#     else:
#         log.info('Output folder {} already exists, this will overwrite/merge '
#                  'whatever is inside.'.format(output_dir))
#
#     input_dir  = Crumb(input_dir, regex=regex, ignore_list=['.*'])
#
#     if not input_dir.has_crumbs():
#         raise ValueError('I am almost sure that this cannot work if you do not '
#                          'use crumb arguments in the input path, got {}.'.format(input_dir))
#
#     acq_folder_arg, last_in_arg = tuple(input_dir.all_args())[-2:]
#     out_arg_names = ['{' + arg + '}' for arg in tuple(input_dir.all_args())[:-1]]
#     output_dir    = Crumb(op.join(output_dir, *out_arg_names), regex=regex, ignore_list=['.*'])
#
#     src_dst = []
#     acquisitions = input_dir.ls(acq_folder_arg, make_crumbs=True)
#     for acq in acquisitions:
#         out_args = acq.arg_values.copy()
#         acq_out  = output_dir.replace(**out_args)
#
#         out_dir  = op.dirname (acq_out.path)
#         out_file = op.basename(acq_out.path) + '.nii.gz'
#         os.makedirs(out_dir, exist_ok=True)
#
#         src_dst.append((acq.split()[0], out_dir, out_file))
#
#     if PARALLEL and ncpus > 1:
#          import multiprocessing as mp
#          pool = mp.Pool(processes=ncpus)
#          results = [pool.apply_async(convert_dcm2nii, args=(dr, ss, dst)) for dr, ss, dst in src_dst]
#          _ = [p.get() for p in results]
#     else:
#          _ = [convert_dcm2nii(path, sess, dst) for path, sess, dst in src_dst]
#
#
# @task
# def decompress_dicoms(input_dir):
#     """ Decompress all *.dcm files recursively found in DICOM_DIR.
#     This uses 'gdcmconv --raw'.
#     It works when 'dcm2nii' shows the `Unsupported Transfer Syntax` error. This error is
#     usually caused by lack of JPEG2000 support in dcm2nii compilation.
#
#     Read more:
#     http://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage#Transfer_Syntaxes_and_Compressed_Images
#
#     Parameters
#     ----------
#     input_dir: str
#         Folder path
#
#     Notes
#     -----
#     The *.dcm files in `input_folder` will be overwritten.
#     """
#     import subprocess
#
#     dcmfiles = sorted(recursive_glob(input_dir, '*.dcm'))
#     for dcm in dcmfiles:
#         cmd = 'gdcmconv --raw -i "{0}" -o "{0}"'.format(dcm)
#         log.debug('Calling {}.'.format(cmd))
#         try:
#             subprocess.check_call(cmd, shell=True)
#         except:
#             pass

# TODO:
# -----------------------
# CORTICAL THICKNESS: http://localhost:8888/lab/tree/Cortical%20thickness%20measures%20statistics.ipynb
# ----------------------
