# -*- coding: utf-8 -*-
"""
Resting state networks post-hoc analysis tools.
"""
from collections import OrderedDict
import os.path as op

import numpy as np
import nilearn.plotting as niplot
import nilearn.masking as nimask
import nilearn.image as niimg


class RSN_templates:
    """ A class to parse and return useful values from the
    RSN templates given by GIFT.

    Parameters
    ----------
    img_file: str
        Path to a img-like RSN templates file.

    txt_file: str
        Path to a text file with the labels for `img_file`.
        The expected structure of each line of
        this file is  <name of the RSN>, <list of indices>.
        For example:

        Basal Ganglia, 21
        Auditory, 17
        Sensorimotor, 7, 23, 24, 38, 56, 29
        Visual, 46, 64, 67, 48, 39, 59
        Default-Mode, 50, 53, 25, 68
        Attentional, 34, 60, 52, 72, 71, 55
        Frontal, 42, 20, 47, 49


    start_from_one: bool
        If True it means that the `txt_file` volume indices start from 1.
        From 0 otherwise.
    """
    def __init__(self, img_file, txt_file, start_from_one=True):
        self._img_file       = img_file
        self._txt_file       = txt_file
        self._start_from_one = start_from_one

        self.network_names = self._network_names()
        self.img           = niimg.load_img(self._img_file)
        self._self_check()

    def iter_rsns(self, mask=None, affine=None):
        """ Yield name and image of each RSN volume from the text file."""
        for idx, name in self.network_names.items():
            yield name, niimg.index_img(self.img, idx)

    def _network_names(self):
        names_idx = self._updt_networks_blob_idx()
        return {idx: name for name, idxs in names_idx.items() for idx in idxs}

    def _self_check(self):
        if len(self.network_names) == self.img.shape[3]:
            return

        print('Then number of volumes in the image is different from the number of labels in the text file.\n'
              'I am going to use only the ones in the text file.')

    def _updt_networks_blob_idx(self):
        """ Read the text file and return a dict[str->List[int]] with network names and blob indices."""
        lines = [l.strip() for l in open(self._txt_file).readlines()]

        netblobs = OrderedDict()
        for l in lines:
            pcs = l.split(',')
            netname = pcs[0]
            blobs   = [int(idx) for idx in pcs[1:]]
            if self._start_from_one:
                blobs = [i-1 for i in blobs]

            netblobs[netname] = blobs
        return netblobs

    def plot_all(self):
        names = rsns.network_names
        for idx, rsn in enumerate(niimg.iter_img(rsns.img)):
            disp = niplot.plot_roi(rsn, title=names.get(idx, None))

    # def correlations_with(self, ics):
    #     """
    #     """
    #     from scipy.stats import pearsonr
    #
    #     correlations = []
    #
    #     for rsn_idx, (name, rsn) in enumerate(allens.iter_rsns()):
    #         for ic_idx, ic in enumerate(niimg.iter_img(ics)):
    #             rsn_transf = niimg.resample_to_img(rsn, ic, copy=True)
    #
    #             rsn_masked = nimask.apply_mask(rsn_transf, mask_file)
    #             ic_masked  = nimask.apply_mask(ic,         mask_file)
    #
    #             r, p = pearsonr(rsn_masked, ic_masked)
    #
    #             correlations.append((rsn_idx, ic_idx, r))
    #
    #     from collections import defaultdict
    #     corr = defaultdict(dict)
    #
    #     for idx, (rsn_idx, ic_idx, r) in enumerate(correlations):
    #         corr[ic_idx+1][rsn_idx] = r
    #
    #     pd.DataFrame.from_dict(corr).to_excel('/home/alexandre/data/thomas/ica_out/8mm/fmri_no-grptemplate_30ICs/allen_rsn_correlations.xls')
    #     [print(name) for idx, name in allens.network_names.items()]

    def _joined_rsn_blobs(self, indices):
        """ Return a NiftiImage containing a merge of the RSN blobs indexed by `indices`."""
        rsn_img = niimg.load_img(self._img_file)
        oimg = niimg.index_img(rsn_img, indices[0])
        for idx in indices[1:]:
            oimg = niimg.math_img('oimg + next', oimg=oimg, next=niimg.index_img(rsn_img, idx))

        return niimg.math_img('oimg.astype(bool)', oimg=oimg)
