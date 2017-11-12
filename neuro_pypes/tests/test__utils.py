# -*- coding: utf-8 -*-
from neuro_pypes._utils import format_pair_list


def test_format_pair_list():

    anat_fbasename = 'anat_hc'

    regexp_subst = [
                     (r"/{anat}_.*corrected_seg8.mat$", "/{anat}_to_mni_affine.mat"),
                     (r"/m{anat}.*_corrected.nii$",     "/{anat}_biascorrected.nii"),
                     (r"/w{anat}.*_biascorrected.nii$", "/{anat}_mni.nii"),
                     (r"/y_{anat}.*nii$",               "/{anat}_to_mni_field.nii"),
                     (r"/iy_{anat}.*nii$",              "/{anat}_to_mni_inv_field.nii"),
                     (r"/mwc1{anat}.*nii$",             "/{anat}_gm_mod_w2tpm.nii"),
                     (r"/mwc2{anat}.*nii$",             "/{anat}_wm_mod_w2tpm.nii"),
                     (r"/mwc3{anat}.*nii$",             "/{anat}_csf_mod_w2tpm.nii"),
                     (r"/mwc4{anat}.*nii$",             "/{anat}_nobrain_mod_w2tpm.nii"),
                     (r"/c1{anat}.*nii$",               "/{anat}_gm.nii"),
                     (r"/c2{anat}.*nii$",               "/{anat}_wm.nii"),
                     (r"/c3{anat}.*nii$",               "/{anat}_csf.nii"),
                     (r"/c4{anat}.*nii$",               "/{anat}_nobrain.nii"),
                     (r"/c5{anat}.*nii$",               "/{anat}_nobrain_mask.nii"),
                   ]

    result = format_pair_list(regexp_subst, anat=anat_fbasename)

    assert(result == [
                      (r"/anat_hc_.*corrected_seg8.mat$", "/anat_hc_to_mni_affine.mat"),
                      (r"/manat_hc.*_corrected.nii$",     "/anat_hc_biascorrected.nii"),
                      (r"/wanat_hc.*_biascorrected.nii$", "/anat_hc_mni.nii"),
                      (r"/y_anat_hc.*nii$",               "/anat_hc_to_mni_field.nii"),
                      (r"/iy_anat_hc.*nii$",              "/anat_hc_to_mni_inv_field.nii"),
                      (r"/mwc1anat_hc.*nii$",             "/anat_hc_gm_mod_w2tpm.nii"),
                      (r"/mwc2anat_hc.*nii$",             "/anat_hc_wm_mod_w2tpm.nii"),
                      (r"/mwc3anat_hc.*nii$",             "/anat_hc_csf_mod_w2tpm.nii"),
                      (r"/mwc4anat_hc.*nii$",             "/anat_hc_nobrain_mod_w2tpm.nii"),
                      (r"/c1anat_hc.*nii$",               "/anat_hc_gm.nii"),
                      (r"/c2anat_hc.*nii$",               "/anat_hc_wm.nii"),
                      (r"/c3anat_hc.*nii$",               "/anat_hc_csf.nii"),
                      (r"/c4anat_hc.*nii$",               "/anat_hc_nobrain.nii"),
                      (r"/c5anat_hc.*nii$",               "/anat_hc_nobrain_mask.nii"),
                     ])