#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PyInvoke tasks examples for usage of workflows from Pypes.
"""

from __future__ import (absolute_import,
                        division,
                        print_function,
                        unicode_literals)

import os
import os.path as op
import logging

import pandas as pd
from hansel import Crumb
from invoke import task
from boyle.files.search  import recursive_glob

from pypes.run import run_debug, run_wf


log = logging.getLogger()

# to get std_brains with some atlases and templates:
# git clone https://github.com/Neurita/std_brains.git
STDB_DIR = '/home/hansel/data/std_brains'
HAMM_DIR = op.join(STDB_DIR, 'atlases', 'hammers')
HAMM_MNI = op.join(HAMM_DIR, 'Hammers_mith_atlas_n30r83_SPM5.nii.gz')

SPM_CANONICAL_BRAIN_2MM = op.join(STDB_DIR, 'templates', 'spm_canonical', 'single_subj_T1_brain.nii.gz')


def verbose_switch(verbose=False):
    if verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.getLogger().setLevel(log_level)


@task
def decompress_dicoms(ctx, input_dir):
    """ Decompress all *.dcm files recursively found in DICOM_DIR.
    This uses 'gdcmconv --raw'.
    It works when 'dcm2nii' shows the `Unsupported Transfer Syntax` error. This error is
    usually caused by lack of JPEG2000 support in dcm2nii compilation.

    Read more:
    http://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage#Transfer_Syntaxes_and_Compressed_Images

    Parameters
    ----------
    input_dir: str
        Folder path

    Notes
    -----
    The *.dcm files in `input_folder` will be overwritten.
    """
    import subprocess

    dcmfiles = sorted(recursive_glob(input_dir, '*.dcm'))
    for dcm in dcmfiles:
        cmd = 'gdcmconv --raw -i "{0}" -o "{0}"'.format(dcm)
        log.debug('Calling {}.'.format(cmd))
        try:
            subprocess.check_call(cmd, shell=True)
        except:
            pass


@task
def dcm2nii(ctx, input_crumb_path, output_dir, regex='fnmatch', ncpus=3):
    """ Convert all DICOM files within `input_crumb_path` into NifTI in `output_folder`.

    Will copy only the NifTI files reoriented by MRICron's dcm2nii command.
    Will rename the NifTI files that are matched with recognized modalities to the short
    modality name from config.ACQ_PATTERNS.

    Parameters
    ----------
    input_dir: str
        A crumb path str indicating the whole path until the DICOM files.
        Example: '/home/hansel/data/{group}/{subj_id}/{session_id}/{acquisition}/{dcm_file}

        The crumb argument just before the last one will be used as folder container reference
        for the DICOM series.

    output_dir: str
        The root folder path where to save the tree of nifti files.
        Example: '/home/hansel/nifti'
        This function will create the same tree as the crumbs in input_crumb_path, hence
        for the example above the output would have the following structure:
        '/home/hansel/nifti/{group}/{subj_id}/{session_id}/{nifti_file}'

        Where {nifti_file} will take the name from the {acquisition} or from the
        patterns in ACQ_PATTERNS in `config.py` file.

    regex: str
        The regular expression syntax you may want to set in the Crumbs.
        See hansel.Crumb documentation for this.

    ncpus: int
        this says the number of processes that will be launched for dcm2nii in parallel.
    """
    from boyle.dicom.convert import convert_dcm2nii

    input_dir  = op.expanduser(input_crumb_path)
    output_dir = op.expanduser(output_dir)

    if not op.exists(output_dir):
        log.info('Creating output folder {}.'.format(output_dir))
        os.makedirs(output_dir)
    else:
        log.info('Output folder {} already exists, this will overwrite/merge '
                 'whatever is inside.'.format(output_dir))

    input_dir  = Crumb(input_dir, regex=regex, ignore_list=['.*'])

    if not input_dir.has_crumbs():
        raise ValueError('I am almost sure that this cannot work if you do not '
                         'use crumb arguments in the input path, got {}.'.format(input_dir))

    acq_folder_arg, last_in_arg = tuple(input_dir.all_args())[-2:]
    out_arg_names = ['{' + arg + '}' for arg in tuple(input_dir.all_args())[:-1]]
    output_dir    = Crumb(op.join(output_dir, *out_arg_names), regex=regex, ignore_list=['.*'])

    src_dst = []
    acquisitions = input_dir.ls(acq_folder_arg, make_crumbs=True)
    for acq in acquisitions:
        out_args = acq.arg_values.copy()
        acq_out  = output_dir.replace(**out_args)

        out_dir  = op.dirname (acq_out.path)
        out_file = op.basename(acq_out.path) + '.nii.gz'
        os.makedirs(out_dir, exist_ok=True)

        src_dst.append((acq.split()[0], out_dir, out_file))

    if ncpus > 1:
         import multiprocessing as mp
         pool = mp.Pool(processes=ncpus)
         results = [pool.apply_async(convert_dcm2nii, args=(dr, ss, dst)) for dr, ss, dst in src_dst]
         _ = [p.get() for p in results]
    else:
         _ = [convert_dcm2nii(path, sess, dst) for path, sess, dst in src_dst]


@task
def clinical_pype(ctx, wf_name="spm_anat_preproc", base_dir="",
                  cache_dir="", output_dir="", settings_file='',
                  plugin="MultiProc", n_cpus=4):
    """ Run the basic pipeline.

    Parameters
    ----------
    wf_name: str

    base_dir: str

    cache_dir: str

    output_dir: str

    year: str or int

    plugin: str

    n_cpus: int
    """
    from pypes.datasets import clinical_crumb_workflow

    data_path = op.join(op.expanduser(base_dir), '{year}', '{subject_id}', '{session_id}', '{image}')
    data_crumb = Crumb(data_path, ignore_list=['.*'])

    atlas_file = HAMM_MNI

    wf = clinical_crumb_workflow(wf_name     = wf_name,
                                 data_crumb  = data_crumb,
                                 cache_dir   = op.abspath(op.expanduser(cache_dir)) if cache_dir else '',
                                 output_dir  = op.abspath(op.expanduser(output_dir)) if output_dir else '',
                                 config_file = settings_file,
                                 params={'atlas_file': atlas_file},
                                 )

    if n_cpus > 1:
        run_wf(wf, plugin=plugin, n_cpus=n_cpus)
    else:
        run_wf(wf, plugin=None)


@task
def run_canica(ctx, input_crumb, output_dir, cache_dir="", mask_file="", algorithm='canica', comps=30, smooth_fwhm=8,
               wf_name="", settings_file=""):
    """ Perform ICA (CanICA or DictLearning) on the files given by `input_crumb`.
    Parameters
    ----------
    input_crumb: str
        Crumb path that will give a list of the input files for ICA.
        The last open argument and its pattern of the `input_crumb` will be used as a reference for the input image
        file for the ICA. So, put a crumb argument with fixed expression in the basename of the path, e.g.:

        `/home/hansel/cobre/{sid}/session_0/{img:rest.nii.gz}`.

    mask_file: str
        Path to a mask file to select the image regions that
        This file must have the same dimensions as all the files listed from `input_crumb`.

    algorithm: str
        Name of the ICA algorithme.
        Choices: 'canica', 'dictlearning'

    comps: int
        Number of components to extract from the ICA.

    Outputs
    -------
    The results will be stored in `output_dir`.
    """
    from functools import partial

    from pypes.config import update_config
    from pypes.io     import build_crumb_workflow
    from pypes.ica    import attach_concat_canica


    # set the configuration parameters
    if settings_file:
        update_config(settings_file)

    # expanduser in inputs paths:
    cache_dir  = op.expanduser(cache_dir)
    output_dir = op.expanduser(output_dir)
    if not cache_dir:
        cache_dir = op.join(output_dir, '.pypes_cache')

    # base folder depending if using MR-PET pipeline or PET-only
    data_crumb = Crumb(input_crumb, ignore_list=['.*'])

    # more configs
    if not wf_name:
        wf_name = algorithm

    if comps:
        update_config({wf_name + '_ica.n_components'  : comps})

    update_config({wf_name + '_ica.algorithm':      algorithm})
    update_config({wf_name + '_ica.mask':           mask_file})
    update_config({wf_name + '_ica.smoothing_fwhm': smooth_fwhm})
    update_config({wf_name + '_ica.do_cca':         True})
    update_config({wf_name + '_ica.standardize':    True})

    update_config({wf_name + '_ica.n_init':         20})
    update_config({wf_name + '_ica.n_jobs':         -1})
    update_config({'plot_ica.bg_img':  SPM_CANONICAL_BRAIN_2MM})

    # the input folder and files
    files_crumb_args = {}
    _, arg_name = data_crumb._last_open_arg()
    files_crumb_args['input_img'] = [(arg_name, data_crumb.patterns.get(arg_name, ""))]

    kwargs = dict()
    kwargs['input_connection'] = 'input_img'
    kwargs['input_node']       = 'selectfiles'

    # build the workflow
    wf = build_crumb_workflow({wf_name: partial(attach_concat_canica, **kwargs)},
                              data_crumb=data_crumb,
                              in_out_kwargs=files_crumb_args,
                              output_dir=output_dir,
                              cache_dir=cache_dir,)

    wf.remove_nodes([wf.get_node('datasink')])

    run_wf(wf)


@task
def run_gift(ctx, input_dir, output_dir, mask_file, zscore_plot=2):
    """ Perform MIALAB's GIFT InfoMax.

    Uses the gift_batch_template.m file to create the GIFT batch input.

    Parameters
    ----------
    input_dir: str

    output_dir: str

    mask_file: str

    zscore_plot: float

    Examples
    --------
    $ inv run_gift -i /home/hansel/data/thomas/ica_in -o /home/hansel/data/thomas/ica_out

    Outputs
    -------
    The results will be stored in output_dir/{algorithm}_{preproc_type}_{group}
    ...
    """
    import os
    import io
    import os.path as op
    import subprocess

    from jinja2 import Template
    from pypes.ica import plot_ica_results

    tmp_file = 'gift_batch_template.m'
    tmp_str = Template(io.open(tmp_file).read())
    tmp_str = tmp_str.render(input_dir=input_dir,
                             output_dir=output_dir,
                             mask_file=mask_file)
    batch_file = op.abspath('gift_filled_template.m')

    io.open(batch_file, 'w').write(tmp_str)

    cmd = 'matlab -nodesktop -nosplash -r "icatb_batch_file_run(\'{}\'); exit();"'.format(batch_file)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    os.remove(batch_file)

    bg_img = op.expanduser(SPM_CANONICAL_BRAIN_2MM)
    return plot_ica_results(output_dir, application='gift',
                            mask_file=mask_file, zscore=zscore_plot, bg_img=bg_img)


@task
def run_sbm(ctx, input_dir, output_dir, mask_file, zscore_plot=2):
    """ Perform MIALAB's SBM InfoMax.
    Uses the sbm_batch_template.m file to create the SBM batch input.

    Parameters
    ----------
    input_dir: str

    output_dir: str

    mask_file: str

    zscore_plot: float

    Examples
    --------
    $ inv run_sbm -i /home/hansel/data/thomas/ica_in -o /home/hansel/data/thomas/ica_out

    Outputs
    -------
    The results will be stored in output_dir/{algorithm}_{preproc_type}_{group}
    ...
    """
    import os
    import io
    import os.path as op
    import subprocess

    from jinja2 import Template
    from pypes.ica import plot_ica_results

    input_glob = op.join(input_dir, '*.nii')

    tmp_file = 'sbm_batch_template.m'
    tmp_str = Template(io.open(tmp_file).read())
    tmp_str = tmp_str.render(input_glob=input_glob,
                             output_dir=output_dir,
                             out_prefix='sbm_',
                             mask_file=mask_file)

    batch_file = op.abspath('sbm_filled_template.m')
    io.open(batch_file, 'w').write(tmp_str)

    cmd = 'matlab -nodesktop -nosplash -r "icatb_batch_file_run(\'{}\'); exit();"'.format(batch_file)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    os.remove(batch_file)

    bg_img = op.expanduser(SPM_CANONICAL_BRAIN_2MM)
    return plot_ica_results(output_dir, application='sbm', mask_file=mask_file,
                            zscore=zscore_plot, bg_img=bg_img)


@task
def cobre_pype(ctx, wf_name="spm_anat_rest_preproc", base_dir="", cache_dir="", output_dir="", settings_file="",
               plugin=None, n_cpus=4):
    """ Run the

    ParametersA
    ----------
    wf_name: str

    base_dir: str
        Base path to where the data is

    cache_dir: str

    output_dir: str

    year: str or int

    plugin: str

    n_cpus: int
    """
    from pypes.datasets import cobre_crumb_workflow

    data_path = op.join(op.expanduser(base_dir), '{subject_id}', 'session_1', '{modality}', '{image}')
    data_crumb = Crumb(data_path, ignore_list=['.*'])

    wf = cobre_crumb_workflow(wf_name     = wf_name,
                              data_crumb  = data_crumb,
                              cache_dir   = op.abspath(op.expanduser(cache_dir)) if cache_dir else '',
                              output_dir  = op.abspath(op.expanduser(output_dir)) if output_dir else '',
                              config_file = settings_file,
                              params={'atlas_file': HAMM_MNI},
                             )

    run_wf(wf, plugin=plugin, n_cpus=n_cpus)


@task(autoprint=True)
def plot_ica_results(ctx, ica_result, application, mask_file='', mode='+-', zscore=0, bg_img=None):
    """ Use nilearn through pypes to plot results from CanICA, DictLearning, MIALAB GIFT or SBM,
    given the ICA result folder path.

    Parameters
    ----------
    ica_result: str
        Path to the ICA output folder or the ICA components volume file.

    application: str
        Choicese: ('canica', 'sbm', 'gift', 'gift-group')

    mask_file: str
        Path to the brain mask file to be used for thresholding.

    mode: str
        Choices: '+' for positive threshold,
                 '+-' for positive and negative threshold and
                 '-' for negative threshold.

    zscore: int
        Value of the Z-score thresholding.

    bg_img: str
        Path to a background image.
        If empty will use the SPM canonical brain image at 2mm.
    """
    from pypes.ica import plot_ica_results

    if bg_img is None:
        bg_img = op.expanduser(SPM_CANONICAL_BRAIN_2MM)

    return plot_ica_results(ica_result,
                            application=application,
                            mask_file=mask_file,
                            zscore=float(zscore),
                            mode=mode,
                            bg_img=bg_img)


@task
def motion_stats_sheet(ctx, motion_file_cr, crumb_fields, out_path):
    """ Create in `out_path` an Excel spreadsheet with some of the motion statistics obtained from the
    `statistics_files` output of the nipype.RapidArt found in the hansel.Crumb `motion_file_cr`.

    Parameters
    ----------
    motion_file_cr: str

    crumb_fields: list of str

    out_path: str

    Examples
    --------
    >>> inv motion_stats_sheet \
    >>> --motion-file-cr "/home/hansel/data/out/{group}/{patient_id}/{session}/rest/artifact_stats/motion_stats.json" \
    >>> --crumb-fields "['group', 'patient_id', 'session']" \
    >>> --out-path "/home/hansel/data/motion_stats.xls"
    """
    import json
    from collections import OrderedDict

    from hansel import Crumb

    def get_motion_record(mtn_file_cr, crumb_fields):
        """ Return an OrderedDict of the information found in the `mtn_file_cr` and also
        `crumb_fields` Crumb argument values."""
        stats = json.load(open(str(mtn_file_cr)))

        outliers = stats[1]
        motion_norm = stats[3]['motion_norm']

        #outliers_hdr = list(outliers.keys())
        motion_hdr   = ['{}_motion_norm'.format(k) for k in motion_norm.keys()]

        mtn_record = OrderedDict()
        for fn in crumb_fields:
            mtn_record[fn] = mtn_file_cr[fn][0]

        mtn_record.update(outliers)

        for hdr, fn in zip(motion_hdr, motion_norm):
            mtn_record[hdr] = motion_norm[fn]

        return mtn_record

    # process the input
    motion_file_cr = Crumb(motion_file_cr)
    crumb_fields   = [crf.strip() for crf in crumb_fields[1:-1].replace("'", "").split(',')]

    # create the motion records
    motionstats = [get_motion_record(stats_file, crumb_fields) for stats_file in motion_file_cr.ls()]

    # create a pandas Dataframe out of it
    df = pd.DataFrame.from_records(motionstats, columns=motionstats[0].keys())

    # save it into an excel file
    df.to_excel(out_path)


@task
def ica_sbm_loadings_sheet(ctx, ica_out_dir, labels_file="", mask="", bg_img=None, zscore=2.,
                           subjid_pat=r'(?P<patid>[a-z]{2}_[0-9]{6})'):
    """
    Save the Excel loadings files in the `ica_out_dir`.
    One file is `subject_loadings.xls` which has the loadings as is, with the subjects IDs and group.
    The other file is `subject_group_loadings.xls` which has the loading signs changed according to
    the average correlation value of the "main" region of each of the IC spatial maps.

    Parameters
    ----------
    ica_out_dir: str
        Path to the SBM ICA analysis output folder.

    labels_file: str
        A CSV file with two columns: "subject_id" and "group".
        The subject_ids must be in the paths contained in the Subject.mat
        file and match the `subjid_pat` argument.

    mask: str
        Path to a mask file to select only brain area from the IC spatial map.

    bg_img: str
        A background image for the blob plots check report, to verify that the blobs
        taken into account for the loadings signs are correct.

    zscore: float
        Value to threshold the IC spatial maps to obtain the IC spatial map "main" region.

    subjid_pat: regext str
        A search regex pattern that returns one group element that
        contains the subject id.
        This will be used to search for subject_id in the file paths
        contained in the Subjects.mat file.
    """
    from pypes.ica.plotting import SBMICAResultsPlotter

    rawloadings_filename   = 'subject_loadings.xls'
    grouploadings_filename = 'subject_weighted_loadings.xls'
    check_blob_plot        = 'check_sign_blobs.png'

    plotter = SBMICAResultsPlotter(ica_out_dir)
    plotter.fit(mask_file=mask, mode='+-', zscore=zscore)

    # generate and save the simple loadings sheet
    sdf = plotter.simple_loadings_df(group_labels_file=labels_file, subjid_pat=subjid_pat)
    sdf.to_excel(op.join(ica_out_dir, rawloadings_filename))

    # generate and save the group-processed loadings sheet
    pdf = plotter.weighted_loadings_df(group_labels_file=labels_file, subjid_pat=subjid_pat)
    pdf.to_excel(op.join(ica_out_dir, grouploadings_filename))

    # plot blobs over IC maps for checking
    check_blob_plot = op.join(ica_out_dir, check_blob_plot)
    plotter.plot_icmaps_and_blobs(check_blob_plot, bg_img=bg_img)
