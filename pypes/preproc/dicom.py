# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
 Helper functions and nipype interface for
 reading DICOM files, specially from Siemens acquisitions
"""


def dcm_ascii_hdr(dcm_file):
    """ Return the CSA ASCII header from a Siemens DICOM file.

    References
    ----------
    http://www.nmr.mgh.harvard.edu/~greve/dicom-unpack
    https://groups.google.com/forum/#!msg/pydicom/-CmY-yK8Y7A/XaR4cUpizvgJ

    Note
    ----
    This only works for Siemens DICOM files.
    """
    def get_csa_header(dcm_file):
        import pydicom
        from nibabel.nicom import csareader as csar
        dcm_data = pydicom.read_file(dcm_file)
        return csar.get_csa_header(dcm_data, 'series')

    csa = get_csa_header(dcm_file)
    return csa['tags']['MrPhoenixProtocol']['items'][0]


def split_dcm_ahdr(ahdr):
    """ Return the protocol meta data and the ascconv part separately
    in list of lines.

    Parameters
    ----------
    ahdr: str

    Returns
    -------
    meta: list of str

    ascconv: list of str

    References
    ----------
    https://groups.google.com/forum/#!msg/pydicom/-CmY-yK8Y7A/XaR4cUpizvgJ

    Note
    ----
    This only works for Siemens DICOM files.
    """
    parts = ahdr.split('### ASCCONV BEGIN ###')

    meta   = parts[0].split('\n')
    ascconv = parts[1].split('### ASCCONV END ###')[0].split('\n')

    return meta, ascconv