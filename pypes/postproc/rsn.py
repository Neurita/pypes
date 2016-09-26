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


from collections import OrderedDict, defaultdict
from boyle.more_collections import DefaultOrderedDict

from scipy.stats import pearsonr
import nilearn.image as niimg
import pandas as pd


class RestingNetworksTemplates:
    """ A class to parse and return useful values from RSN templates.

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
    
    def iter_rsns(self):
        """ Yield name and image of each RSN volume from the text file."""
        for idx, name in self.network_names.items():
            yield name, niimg.index_img(self.img, idx)
    
    def _network_names(self):
        names_idx = self._updt_networks_blob_idx()
        return OrderedDict([(idx, name) for name, idxs in names_idx.items() for idx in idxs])
      
    def _self_check(self):
        if len(self.network_names) == self.img.shape[3]:
            return

        print('The number of volumes in the image is different from the number of labels in the text file.\n'
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
 
    def _joined_rsn_blobs(self, indices):
        """ Return a NiftiImage containing a merge of the RSN blobs indexed by `indices`."""
        rsn_img = niimg.load_img(self._img_file)
        oimg = niimg.index_img(rsn_img, indices[0])
        for idx in indices[1:]:
            oimg = niimg.math_img('oimg + next', oimg=oimg, next=niimg.index_img(rsn_img, idx))
      
        return niimg.math_img('oimg.astype(bool)', oimg=oimg)

    def correlations(self, ic_imgs, mask_file):
        """ Correlation values of each RSN to each IC map in `ic_imgs` masked by `mask_file`.
        
        Parameters
        ----------
        ic_imgs: list of niimg-like
        
        mask_file: niimg-like
        
        Returns
        -------
        corrs_df: pandas.DataFrame
        """
        correlations = DefaultOrderedDict(dict)

        n_ics = ic_imgs.shape[-1]
        mask_trnsf = niimg.resample_to_img(mask_file, niimg.index_img(ic_imgs, 0), interpolation='nearest',
                                           copy=True)

        for rsn_idx, (name, rsn) in enumerate(self.iter_rsns()):
            rsn_transf = niimg.resample_to_img(rsn, niimg.index_img(ic_imgs, 0), copy=True)
            rsn_masked = nimask.apply_mask(rsn_transf, mask_trnsf)

            for ic_idx, ic in enumerate(niimg.iter_img(ic_imgs)):
                ic_masked  = nimask.apply_mask(ic, mask_trnsf)
                r, p = pearsonr(rsn_masked, ic_masked)
                correlations[ic_idx+1][rsn_idx] = r

        corrs_df = pd.DataFrame.from_dict(correlations)

        corrs_df['IC']  = allens.network_names.keys()
        corrs_df['IC']  = corrs_df['IC'] + 1

        corrs_df['RSN'] = allens.network_names.values()      
        
        return corrs_df
    
    def goodness_of_fit(self, ic_imgs, mask_file, rsn_thr=4.0):
        """ Goodness-of-fit values described as in Zhou et al., 2010, Brain.
        
        Parameters
        ----------
        ic_imgs: list of niimg-like
        
        mask_file: niimg-like
        
        rsn_thr: float
        
        Returns
        -------
        gof_df: pandas.DataFrame

        Notes
        -----
        These ICN templates were thresholded at a z-score 4.0 to be visually comparable to the
        consistent ICNs published by Damoiseaux et al. (2006). A minor modification of previous
        goodness-of-fit methods (Seeley et al., 2007b, 2009) was included here for template 
        matching, with goodness-of-fit scores calculated by multiplying 
        (i) the average z-score difference between voxels falling within the template and
        voxels falling outside the template; and 
        (ii) the difference in the percentage of positive z-score voxels inside and outside the template.
        
        This goodness-of-fit algorithm proves less vulnerable to inter-subject variability in shape, 
        size, location and strength of each ICN.
        """
        gofs = DefaultOrderedDict(dict)

        # threshold the RSN templates
        if rsn_thr > 0:
            thr_rsns = [(name, spatial_map(rsn, thr=rsn_thr, mode='+-')) for name, rsn in self.iter_rsns()]       
        else:
            thr_rsns = self.iter_rsns()
        
        # for each RSN template
        for rsn_idx, (name, rsn) in enumerate(thr_rsns):          
           
            # for each IC map
            for ic_idx, ic in enumerate(niimg.iter_img(ic_imgs)):

                subj_ic = filter_icc(ic, thr=rsn_thr, mask=mask_file)

                # prepare the RSN masks
                rsn_brain_mask = niimg.resample_to_img(mask_file, rsn, interpolation='nearest')
                rsn_in  = niimg.math_img('np.abs(img) > 0', img=rsn)
                rsn_out = niimg.math_img('img == 0',        img=rsn)

                rsn_out = niimg.math_img('mask * img', mask=rsn_brain_mask, img=rsn_out)

                rsn_in  = niimg.resample_to_img(rsn_in,  subj_ic, interpolation='nearest')
                rsn_out = niimg.resample_to_img(rsn_out, subj_ic, interpolation='nearest')

                # apply the mask
                zscore_in  = niimg.math_img('mask * img', mask=rsn_in,  img=subj_ic).get_data()
                zscore_out = niimg.math_img('mask * img', mask=rsn_out, img=subj_ic).get_data()

                #gof_term1
                # calculate the the average z-score difference between voxels falling 
                # within the template and voxels falling outside the template
                gof_term1  = zscore_in.mean() - zscore_out.mean()

                #gof_term2
                # the difference in the percentage of positive z-score voxels inside and outside the template.
                n_pos_zscore_in  = np.sum(zscore_in  > 0)
                n_pos_zscore_out = np.sum(zscore_out > 0)
                n_pos_zscore_tot = n_pos_zscore_in + n_pos_zscore_out
                
                if n_pos_zscore_tot != 0:
                    n_pos_zscore_pcnt = 100 / n_pos_zscore_tot
                    gof_term2 = (n_pos_zscore_in - n_pos_zscore_out) * n_pos_zscore_pcnt
                else:
                    gof_term2 = 0
                
                # global gof
                gof = gof_term1 * gof_term2

                # add the result
                gofs[ic_idx+1][rsn_idx] = gof
        
        return gofs
