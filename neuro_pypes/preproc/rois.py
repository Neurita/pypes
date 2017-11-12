# -*- coding: utf-8 -*-
"""
Nipype processing nodes for ROI images
"""
import os.path as op


from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec
from nipype.utils.filemanip import split_filename

import nibabel as nib
from   boyle.nifti.roi import drain_rois

from ..utils.files import niftiimg_out


@niftiimg_out
def drain_rois_img(img):
    return drain_rois(img)


class DrainROIsInputSpec(BaseInterfaceInputSpec):
    img = File(exists=True, desc='volume with ROIS to be drained', mandatory=True)


class DrainROIsOutputSpec(TraitedSpec):
    out = File(exists=True, desc="the volume with the drained ROIs")


class DrainROIs(BaseInterface):
    input_spec  = DrainROIsInputSpec
    output_spec = DrainROIsOutputSpec

    def _run_interface(self, runtime):
        fname = self.inputs.img
        nu_img = drain_rois_img(fname)

        _, base, _ = split_filename(fname)
        nib.save(nu_img, base + '_drained.nii.gz')

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        fname = self.inputs.img
        _, base, _ = split_filename(fname)
        outputs["out"] = op.abspath(base + '_drained.nii.gz')
        return outputs
