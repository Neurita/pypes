# -*- coding: utf-8 -*-
"""
Nipype interfaces to canica and dictlearning in nilearn.decomposition
"""
import os.path as op

import numpy as np

from nipype.interfaces.base import (BaseInterface,
                                    BaseInterfaceInputSpec,
                                    TraitedSpec,
                                    InputMultiPath,
                                    OutputMultiPath,
                                    traits)
from nilearn.decomposition import CanICA, DictLearning

from ...utils import get_trait_value


class CanICAInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(traits.File(desc="NifTI image file(s) from where to extract the data. \n"
                                               "If more than one, all should be spatially normalized.",
                                          exists=True, mandatory=True))
    mask = traits.File(desc="Mask to be used on data. If an instance of masker is passed, then its mask will be used.\n"
                            " If no mask is given, it will be computed automatically by a MultiNiftiMasker \n"
                            "with default parameters.",
                       exists=True)
    algorithm = traits.Enum(['canica', 'dictlearning'], desc="Desired nilearn ICA method.")
    confounds = traits.File(desc="CSV file path."
                                 "This parameter is passed to nilearn.signal.clean. "
                                 "Please see the related documentation for details.",
                            exists=True)
    n_components = traits.Int(desc="Number of components to extract.", default_value=20, usedefault=True)
    n_init = traits.Int(desc="CanICA only: The number of times the fastICA algorithm is restarted.", default_value=10, usedefault=True)
    n_epochs = traits.Int(desc="DictLearning only: Number of epochs the algorithm should run on the data.", default_value=1, usedefault=True)
    alpha = traits.Float(desc="DictLearning only: Sparsity controlling parameter.", default_value=1.0, usedefault=True)
    do_cca = traits.Bool(desc="CanICA only: Indicate if a Canonical Correlation Analysis must be run after the PCA.",
                         default_value=True, usedefault=True)
    memory = traits.Str(desc="Used to cache the masking process. By default, no caching is done. "
                        "If a string is given, it is the path to the caching directory.",
                        default_value='.', usedefault=True,)
    memory_level = traits.Int(desc="Rough estimator of the amount of memory used by caching. "
                              "Higher value means more memory for caching.", default_value=1, usedefault=True)
    threshold = traits.Either(traits.Bool, traits.Float,
                              desc="CanICA only: If None, no thresholding is applied.\n"
                                   "If ‘auto’, then we apply a thresholding that will keep the n_voxels, \n"
                                   "more intense voxels across all the maps, n_voxels being the number of voxels \n"
                                   "in a brain volume. A float value indicates the ratio of voxels to keep "
                                   "(2. means that the maps will together have 2 x n_voxels non-zero voxels ).",)
    random_state = traits.Int(desc="Pseudo number generator state used for random sampling.",)
    smoothing_fwhm = traits.Float(desc="If smoothing_fwhm is defined, it gives the full-width half maximum in "
                                       "millimeters of the spatial smoothing to apply to the signal.",)
    standardize = traits.Bool(desc="If standardize is True, the time-series are centered and normed: "
                                   "their mean is put to 0 and their variance to 1 in the time dimension.",
                              default_value=True, usedefault=True)
    n_jobs = traits.Int(desc="The number of CPUs to use to do the computation. -1 means 'all CPUs', "
                             "-2 'all CPUs but one', and so on.",
                        default_value=1, usedefault=True)


class CanICAOutputSpec(TraitedSpec):
    components = traits.File(desc="A nifti file with the reconstructed volume for each loading.")
    score      = traits.File(desc="Numpy txt file that holds the score for each subjects."
                                  "Score is two dimensional if per_component is True."
                                  "First dimension is squeezed if the number of subjects is one")
    loadings   = OutputMultiPath(traits.File(desc="For each subject, each sample, loadings for each "
                                                  "decomposition components shape: "
                                                  "number of subjects * (number of scans, number of regions))"))


class CanICAInterface(BaseInterface):
    """ Nipype Interface to NiLearn methods to perform Canonical Independent Component Analysis.

    For more information look at: nilearn.decomposition.CanICA
    """
    input_spec = CanICAInputSpec
    output_spec = CanICAOutputSpec

    def _run_interface(self, runtime):
        algorithm         = get_trait_value(self.inputs, 'algorithm',      default='canica')
        mask              = get_trait_value(self.inputs, 'mask')
        n_components      = get_trait_value(self.inputs, 'n_components')
        do_cca            = get_trait_value(self.inputs, 'do_cca')
        smoothing_fwhm    = get_trait_value(self.inputs, 'smoothing_fwhm', default=None)
        standardize       = get_trait_value(self.inputs, 'standardize',    default=None)
        threshold         = get_trait_value(self.inputs, 'threshold',      default=None)
        random_state      = get_trait_value(self.inputs, 'random_state',   default=None)
        n_init            = get_trait_value(self.inputs, 'n_init')
        n_jobs            = get_trait_value(self.inputs, 'n_jobs')
        n_epochs          = get_trait_value(self.inputs, 'n_epochs')
        alpha             = get_trait_value(self.inputs, 'alpha')
        memory            = get_trait_value(self.inputs, 'memory')
        memory_level      = get_trait_value(self.inputs, 'memory_level')
        confounds         = get_trait_value(self.inputs, 'confounds')

        # init the estimator
        if algorithm == 'canica':
            self._estimator = CanICA(mask=mask,
                                     n_components=n_components,
                                     threshold=threshold,
                                     random_state=random_state,
                                     standardize=standardize,
                                     smoothing_fwhm=smoothing_fwhm,
                                     do_cca=do_cca,
                                     verbose=1,
                                     n_init=n_init,
                                     memory=memory,
                                     memory_level=memory_level,
                                     n_jobs=n_jobs,
                                     )

        elif algorithm == 'dictlearning':
            self._estimator = DictLearning(mask=mask,
                                           n_components=n_components,
                                           random_state=random_state,
                                           standardize=standardize,
                                           smoothing_fwhm=smoothing_fwhm,
                                           verbose=1,
                                           n_epochs=n_epochs,
                                           alpha=alpha,
                                           memory=memory,
                                           memory_level=memory_level,
                                           n_jobs=n_jobs,
                                           )

        # set output file names
        self._estimator_name = algorithm
        self._confounds = confounds

        self._reconstructed_img_file = '{}_resting_state.nii.gz'.format(self._estimator_name)
        self._score_file   = '{}_score.txt'.format(self._estimator_name)
        self._loading_file = '{}_{}_loading.txt'

        # fit and transform
        self._estimator.fit(self.inputs.in_files, confounds=self._confounds)
        self._score    = self._estimator.score    (self.inputs.in_files, confounds=self._confounds)
        self._loadings = self._estimator.transform(self.inputs.in_files, confounds=self._confounds)

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()

        masker = self._estimator.masker_
        # Drop output maps to a Nifti file
        components_img = masker.inverse_transform(self._estimator.components_)
        components_img.to_filename(self._reconstructed_img_file)

        # save the score array
        if isinstance(self._score, float):
            with open(self._score_file, 'w') as f:
                f.write(str("%.10f" % self._score))
        else:
            np.savetxt(self._score_file, self._score, fmt='%.10f')

        # save the loadings files
        self._loading_files = []
        for idx, loadings in enumerate(self._loadings):
            loading_file = self._loading_file.format(self._estimator_name, idx)
            np.savetxt(loading_file, loadings, fmt='%.10f')
            self._loading_files.append(loading_file)

        outputs['components'] = op.abspath(self._reconstructed_img_file)
        outputs['score']      = op.abspath(self._score_file)
        outputs['loadings']   = [op.abspath(lf) for lf in self._loading_files]
        return outputs