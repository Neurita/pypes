# -*- coding: utf-8 -*-
"""
Nipype interfaces to calculate connectivity measures using nilearn.
"""
import os.path as op

import numpy as np
import nilearn.connectome
from   nilearn.input_data     import NiftiMapsMasker, NiftiLabelsMasker
#from   boyle.nifti.roi        import get_roilist_from_atlas
from   nipype.interfaces.base import (BaseInterface,
                                      TraitedSpec,
                                      InputMultiPath,
                                      BaseInterfaceInputSpec,
                                      traits,)

from ...utils import get_trait_value


class ConnectivityCorrelationInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(traits.File(desc="NifTI image file(s) from where to extract the data. \n"
                                               "If more than one (3D volumes), all should be spatially normalized.",
                                          exists=True, mandatory=True))
    atlas_file = traits.File(desc="Atlas image file defining the connectivity ROIs.\n"
                                  "Must be spatially normalized to in_files.",
                             exists=True, mandatory=True)
    atlas_type = traits.Enum("probabilistic", "labels",
                             desc="The type of atlas.",
                             default="labels")

    # masker options
    smoothing_fwhm = traits.Float(desc="If smoothing_fwhm is defined, it gives the full-width half maximum in "
                                       "millimeters of the spatial smoothing to apply to the signal.",)
    standardize = traits.Bool(desc="If standardize is True, the time-series are centered and normed: "
                                   "their mean is put to 0 and their variance to 1 in the time dimension.",
                              default_value=False)
    resampling_target = traits.Enum("mask", "maps", "data", "labels", "",
                                    desc="Gives which image gives the final shape/size. "
                                         "This depends on the `atlas_type`. "
                                         "For 'probabilistic' you must use 'mask', 'maps' or None; while for"
                                         "'labels' you must use 'data', 'labels' or None."
                                         "Have a look on nilearn docs for more information.")

    # connectome options
    kind = traits.Enum ("correlation", "partial correlation", "tangent", "covariance", "precision",
                        desc="The connectivity matrix kind.", default='covariance')


class ConnectivityCorrelationOutputSpec(TraitedSpec):
    connectivity = traits.File(desc="Numpy text file with the connectivity matrix.")
    timeseries   = traits.File(desc="Numpy text file with the time-series or 4th dimension data matrix extracted "
                                    "from the atlas ROIs.")


class ConnectivityCorrelationInterface(BaseInterface):
    """ Nipype Interface to NiLearn methods to calculate connectivity matrices using 4D data and a
    spatially normalized ROI atlas.

    For more information look at: nilearn.connectome.ConnectivityMeasure
    """
    input_spec = ConnectivityCorrelationInputSpec
    output_spec = ConnectivityCorrelationOutputSpec

    def _run_interface(self, runtime):
        atlas_type        = get_trait_value(self.inputs, 'atlas_type')
        conn_kind         = get_trait_value(self.inputs, 'kind',)
        #rois_list         = get_trait_value(self.inputs, 'rois_list',         default=None)
        smoothing_fwhm    = get_trait_value(self.inputs, 'smoothing_fwhm',    default=None)
        standardize       = get_trait_value(self.inputs, 'standardize',       default=None)
        resampling_target = get_trait_value(self.inputs, 'resampling_target', default=None)

        self._time_series_file = op.abspath('conn_timeseries.txt')
        self._conn_mat_file    = op.abspath('connectivity.txt')

        ## TODO: add parameter to choose the ROI labels to be used.
        # if rois_list is None:
        #     rois_vals = get_roilist_from_atlas(self.inputs.atlas_file)
        # else:
        #     try:
        #         rois_vals = np.loadtxt(rois_list, dtype=int)
        #     except:
        #         raise IOError('Error reading ROIs list file {}.'.format(rois_list))

        in_files = self.inputs.in_files
        if len(in_files) == 1:
            in_files = in_files[0]

        if atlas_type == 'probabilistic':
            AtlasMasker = NiftiMapsMasker
        elif atlas_type == 'labels':
            AtlasMasker = NiftiLabelsMasker

        masker = AtlasMasker(self.inputs.atlas_file,
                             standardize=standardize,
                             smoothing_fwhm=smoothing_fwhm,
                             resampling_target=resampling_target,
                             memory='nilearn_cache',
                             verbose=5)

        self._time_series = masker.fit_transform(in_files)

        conn_measure   = nilearn.connectome.ConnectivityMeasure(kind=conn_kind)
        self._conn_mat = conn_measure.fit_transform([self._time_series])

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()

        np.savetxt(self._time_series_file, self._time_series,        fmt='%.10f')
        np.savetxt(self._conn_mat_file,    self._conn_mat.squeeze(), fmt='%.10f')

        outputs['timeseries'  ] = self._time_series_file
        outputs['connectivity'] = self._conn_mat_file
        return outputs

