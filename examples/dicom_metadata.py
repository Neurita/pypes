#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Processing tasks for clinical data extracted from DICOM files.
"""
from __future__ import (absolute_import,
                        division,
                        print_function,
                        unicode_literals)

import os
import os.path as op
import logging
from collections import OrderedDict

import pandas as pd
from hansel import Crumb
from invoke import task


log = logging.getLogger()


def verbose_switch(verbose=False):
    if verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.getLogger().setLevel(log_level)


@task(autoprint=True)
def subj_data_from_dicoms(ctx, crumb_path, arg_name, verbose=False):
    """ Print a list of folder_name -> NUK id. The NUK ID is calculated from
    the first DICOM file found in the end of the `dicom_path`.

    Parameters
    ----------
    crumb_path: str
        Path with Crumbs to the DICOM files, e.g.,
        /home/hansel/data/{subj_id}/{session}/{acq}/{dcm_file}

    arg_name: str
        Name of the argument in `dicom_path` of the subj_id

    Returns
    -------
    subj_data: dict of subj records
        A dict with records of the information extracted from the DICOM files
        as well as the calculated NUK Pseudonym.
    """
    if verbose:
        verbose_switch(verbose)

    crumb = Crumb(op.expanduser(op.abspath(crumb_path)), ignore_list=['.*'])
    if not crumb.has_crumbs():
        raise ValueError('Expected a path with crumb arguments, e.g., '
                         '"/home/hansel/data/{group}/{sid}/{session}"')
    subj_nuks = []
    for path in crumb.ls(arg_name):
        log.info('Reading DICOMs in {}.'.format(path))
        subj_path = path.split()[0]
        subj = _read_dcm_until_valid(subj_path)
        if subj is None:
            log.info('Could not find a valid DICOM in {}.'.format(subj_path))
        else:
            subj_nuks.append(subj)

    return subj_nuks


def _read_dcm_until_valid(subj_path):
    """ Look for all DCM files within subj_path until read a valid
    data using _read_subjdata_from_dcm."""
    from boyle.dicom.utils import get_dicom_files

    dcm_iter  = get_dicom_files(subj_path)
    for dcm in dcm_iter:
        try:
            subj = _read_subjdata_from_dcm(subj_path, dcm)
        except Exception as exc:
            log.debug(str(exc))
            continue
        else:
            return subj


def _read_subjdata_from_dcm(subj_path, dcm_obj):
    """ Read Patient data fields in `dcm_obj` and return a dictionary."""
    log.debug('Reading {}.'.format(subj_path))
    subj = OrderedDict()
    subj['name']             = dcm_obj.PatientName.given_name
    subj['surname']          = dcm_obj.PatientName.family_name
    subj['birthdate']        = pd.Timestamp(dcm_obj.PatientBirthDate)
    subj['1st_session_date'] = pd.Timestamp(dcm_obj.AcquisitionDate) # + ' ' + dcm.AcquisitionTime)
    subj['dcm_folder']       = subj_path
    return subj


@task(autoprint=True)
def create_group_name_conversion(ctx, crumb_path, arg_name, outfile):
    """ Create a .csv file with the information extracted with `nuk_ids_from_dicoms`.

    Parameters
    ----------
    crumb_path: str
        Path with Crumbs to the DICOM files, e.g.,
        /home/hansel/data/{subj_id}/{session}/{acq}/{dcm_file}

    arg_name: str
        Name of the argument in `dicom_path` of the subj_id

    outfile: str

    Returns
    -------
    subj_df: pandas DataFrame
        A DataFrame with the information extracted from the DICOM files
        as well as the calculated NUK Pseudonym.
    """
    subj_data = subj_data_from_dicoms(ctx, crumb_path, arg_name)

    subj_df = pd.DataFrame.from_records(subj_data)
    subj_df.to_csv(outfile, sep=',', index=False)

    log.info('File created in {}.'.format(op.abspath(outfile)))

    return subj_df


@task
def rename_to_nuk_id_from_dicoms(ctx, crumb, arg, verbose_only=False):
    """ Print a list of folder_name -> NUK id. The NUK ID is calculated from
    the first DICOM file found in the end of the `dicom_path`.

    Parameters
    ----------
    crumb: str
        Path with Crumbs to the DICOM files, e.g.,
        /home/hansel/data/{subj_id}/{session}/{acq}/{dcm_file}

    arg: str
        Name of the argument in `dicom_path` of the subject identification.

    verbose_only: bool

    Returns
    -------
    subj_data: dict of subj records
        A dict with records of the information extracted from the DICOM files
        as well as the calculated NUK Pseudonym.
    """
    subj_data = subj_data_from_dicoms(ctx, crumb, arg)

    id_path = crumb.split(arg)[0] + arg + '}'

    rename_to_nukid(crumb_path=id_path,
                    arg_name=arg,
                    subj_data=subj_data,
                    verbose_only=verbose_only)

    return subj_data


@task(autoprint=True)
def rename_to_nuk_id_from_csv(ctx, crumb_path, arg_name, data_file, verbose_only=False):
    """ Print a list of folder_name -> NUK id. The NUK ID is calculated from
    the first DICOM file found in the end of the `dicom_path`.

    Parameters
    ----------
    crumb_path: str
        Path with Crumbs to the subject folders files, e.g.,
        /home/hansel/data/{subj_id}

    arg_name: str
        Name of the argument in `crumb_path` of the subject identification.

    data_file: str
        Path to a csv file with the necessary data to run rename_to_nukid
        for each row.

    verbose_only: bool

    Returns
    -------
    src_dsts: list of 2-tuples of str
    """
    subj_data = pd.read_csv(data_file)

    return rename_to_nukid(crumb_path=crumb_path,
                           arg_name=arg_name,
                           subj_data=subj_data,
                           verbose_only=verbose_only)


def rename_to_nukid(crumb_path, arg_name, subj_data, verbose_only=False):
    """ Rename the folders at the `arg_name` level using `subj_data` records.
    Will rename from subj_data['DCM Folder'] to subj_data['NUK Pseudonym'].

    Parameters
    ----------
    crumb_path: str
        Path with Crumbs to the DICOM files, e.g.,
        /home/hansel/data/{subj_id}

    arg_name: str
        Name of the argument in `dicom_path` of the subject identification.
        These names should be the same as the ones in 'DCM Folder' value.

    outfile: str

    verbose_only: bool

    Returns
    -------
    src_dsts: list of 2-tuples of str
    """

    def rename_many(src_dsts, verbose_only=False):
        """ For each 2-tuple in src_dsts of file/folder paths will rename the first element to the second.

        Parameters
        ----------
        src_dsts : list of 2-tuple

        verbose_only: bool
            Will not perform the operation will only print them.
        """
        for (src, dst) in src_dsts:
            if not op.exists(src):
                raise IOError('Could not find source file {}.'.format(src))

            if op.exists(dst):
                if src == dst:
                    continue
                else:
                    raise IOError('Destination path {} already exists.'.format(dst))

            log.info('mv {} -> {}'.format(src, dst))
            if not verbose_only:
                os.rename(src, dst)

    if isinstance(subj_data, pd.DataFrame):
        subjs = subj_data.to_records()
    else:
        subjs = subj_data

    if not Crumb.has_crumbs(Crumb(crumb_path)):
        raise ValueError('Expected a path with crumb arguments, e.g., '
                         '"/home/hansel/data/{group}/{sid}/{session}"')

    crumb = Crumb(op.expanduser(op.abspath(crumb_path)), ignore_list=['.*'])
    src_dsts = []
    for subj in subjs:
        src_crs = crumb.replace(**{arg_name: subj['DCM Folder']}).unfold()

        for src_cr in src_crs:
            dst_args = src_cr.arg_values.copy()
            dst_args[arg_name] = subj['nukid']

            dst_cr = crumb.replace(**dst_args)
            if not op.exists(src_cr.path):
                raise IOError('Could not find folder {} for subject {}.'.format(src_cr, subj))

            if Crumb.has_crumbs(dst_cr.path):
                raise KeyError('The destination path should be fully specified, got {}.'.format(dst_cr))

            src_dsts.append((src_cr.path, dst_cr.path))

    rename_many(src_dsts=src_dsts, verbose_only=verbose_only)

    return src_dsts