#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os

from hansel import Crumb
from hansel.operations import joint_value_map, valuesmap_to_dict
import nipype.pipeline.engine as pe
from nipype.algorithms.misc import Gunzip
from nipype.interfaces import spm, fsl
from nipype.interfaces.utility import IdentityInterface, Function, Select
from nipype.interfaces.io import DataSink
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.interfaces.base import traits

from neuro_pypes.crumb import DataCrumb
from neuro_pypes.preproc.slicetime_params import STCParametersInterface
from neuro_pypes.interfaces.nilearn import math_img
from neuro_pypes.preproc import get_bounding_box
from neuro_pypes._utils import flatten_list
from neuro_pypes.utils import (
    remove_ext,
    joinstrings,
    selectindex,
    extend_trait_list
)


wf_name = 'spm_rest_preprocessing'

#work_dir = os.path.expanduser(f'~/data/neuro_pypes/{wf_name}/')
work_dir = os.path.expanduser(f'/data/neuro_pypes/{wf_name}/')

#input_dir = os.path.expanduser('~/projects/neuro/multimodal_test_data/raw')
input_dir = os.path.expanduser('/data/raw')

output_dir = os.path.join(work_dir, 'out')
cache_dir = os.path.join(work_dir, 'wd')

data_path = os.path.join(os.path.expanduser(input_dir), '{subject_id}', '{session}', '{image}')
data_crumb = Crumb(data_path, ignore_list=['.*'])
crumb_modalities = {
    'anat': [('image', 'anat_hc.nii.gz')],
    'rest': [('image', 'rest.nii.gz')]
}

anat_voxel_sizes = [1, 1, 1]

fmri_smoothing_kernel_fwhm = 8


wf = pe.Workflow(name=wf_name, base_dir=work_dir)

# ------------------------------------------------------------------------------------------------
# DATA INPUT AND SINK
# ------------------------------------------------------------------------------------------------
datasource = pe.Node(
    DataCrumb(crumb=data_crumb, templates=crumb_modalities, raise_on_empty=False),
    name='selectfiles'
)

datasink = pe.Node(
    DataSink(parameterization=False, base_directory=output_dir, ),
    name="datasink"
)

# basic file name substitutions for the datasink
undef_args = datasource.interface._infields
substitutions = [(name, "") for name in undef_args]
substitutions.append(("__", "_"))

# datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions, substitutions)

# Infosource - the information source that iterates over crumb values map from the filesystem
infosource = pe.Node(interface=IdentityInterface(fields=undef_args), name="infosrc")
infosource.iterables = list(valuesmap_to_dict(joint_value_map(data_crumb, undef_args)).items())
infosource.synchronize = True

# connect the input_wf to the datasink
joinpath = pe.Node(joinstrings(len(undef_args)), name='joinpath')

# Connect the infosrc node to the datasink
input_joins = [(name, 'arg{}'.format(arg_no + 1)) for arg_no, name in enumerate(undef_args)]

wf.connect([
    (infosource, datasource, [(field, field) for field in undef_args]),
    (datasource, joinpath, input_joins),
    (joinpath, datasink, [("out", "container")]),
])


# ------------------------------------------------------------------------------------------------
# ANAT
# ------------------------------------------------------------------------------------------------

# T1 preprocessing nodes

# ANTs N4 Bias field correction
# n4 = N4BiasFieldCorrection()
# n4.inputs.dimension = 3
# n4.inputs.bspline_fitting_distance = 300
# n4.inputs.shrink_factor = 3
# n4.inputs.n_iterations = [50, 50, 30, 20]
# n4.inputs.convergence_threshold = 1e-6
# n4.inputs.save_bias = True
# n4.inputs.input_image = traits.Undefined
# biascor = pe.Node(n4, name="bias_correction")

gunzip_anat = pe.Node(Gunzip(), name="gunzip_anat")

# SPM New Segment
spm_info = spm.Info()
priors_path = os.path.join(spm_info.path(), 'tpm', 'TPM.nii')
segment = spm.NewSegment()
tissue1 = ((priors_path, 1), 1, (True,  True),   (True,  True))
tissue2 = ((priors_path, 2), 1, (True,  True),   (True,  True))
tissue3 = ((priors_path, 3), 2, (True,  True),   (True,  True))
tissue4 = ((priors_path, 4), 3, (True,  True),   (True,  True))
tissue5 = ((priors_path, 5), 4, (True,  False),  (False, False))
tissue6 = ((priors_path, 6), 2, (False, False),  (False, False))
segment.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
segment.inputs.channel_info = (0.0001, 60, (True, True))
segment.inputs.write_deformation_fields = [True, True]
segment.inputs.channel_files = traits.Undefined
segment = pe.Node(segment, name="new_segment")

# Apply deformations
normalize_anat = spm.Normalize12(jobtype='write')
normalize_anat.inputs.write_voxel_sizes = anat_voxel_sizes
normalize_anat.inputs.deformation_file = traits.Undefined
normalize_anat.inputs.image_to_align = traits.Undefined
normalize_anat.inputs.write_bounding_box = traits.Undefined
warp_anat = pe.Node(normalize_anat, name="warp_anat")

tpm_bbox = pe.Node(
    Function(function=get_bounding_box, input_names=["in_file"], output_names=["bbox"]),
    name="tpm_bbox"
)
tpm_bbox.inputs.in_file = priors_path

# calculate brain mask from tissue maps
tissues = pe.Node(
    IdentityInterface(fields=["gm", "wm", "csf"], mandatory_inputs=True),
    name="tissues"
)
brain_mask = pe.Node(
    Function(
        function=math_img,
        input_names=["formula", "out_file", "gm", "wm", "csf"],
        output_names=["out_file"],
        imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
        name='brain_mask'
)
brain_mask.inputs.out_file = "tissues_brain_mask.nii.gz"
brain_mask.inputs.formula  = "np.abs(gm + wm + csf) > 0"

# Connect the nodes
wf.connect([
    # input to biasfieldcorrection
#     (datasource, biascor, [("anat", "input_image")]),

    # new segment
#     (biascor,      gunzip_anat, [("output_image", "in_file")]),
    (datasource, gunzip_anat, [("anat", "in_file")]),
    (gunzip_anat,  segment, [("out_file", "channel_files")]),

    # Normalize12
    (segment,   warp_anat,  [("forward_deformation_field", "deformation_file")]),
    (segment,   warp_anat,  [("bias_corrected_images",     "apply_to_files")]),
    (tpm_bbox,  warp_anat,  [("bbox",                      "write_bounding_box")]),

    # brain mask from tissues
    (segment, tissues,[
        (("native_class_images", selectindex, 0), "gm"),
        (("native_class_images", selectindex, 1), "wm"),
        (("native_class_images", selectindex, 2), "csf"),
    ]),

    (tissues, brain_mask, [("gm", "gm"), ("wm", "wm"), ("csf", "csf"),]),

    # output
    (warp_anat, datasink,  [("normalized_files",           "anat.@mni")]),
    (segment,   datasink,  [("modulated_class_images",     "anat.tissues.warped"),
                            ("native_class_images",        "anat.tissues.native"),
                            ("transformation_mat",         "anat.transform.@linear"),
                            ("forward_deformation_field",  "anat.transform.@forward"),
                            ("inverse_deformation_field",  "anat.transform.@inverse"),
                            ("bias_corrected_images",      "anat.@biascor")]),
    (brain_mask, datasink, [("out_file",                   "anat.@brain_mask")]),
])



def _sum_one_to_each(slice_order): # SPM starts count from 1
    return [i+1 for i in slice_order]

def _sum_one(num):
    return num + 1

def _pick_first(sequence):
    return sequence[0]


from nipype.interfaces.nipy.preprocess import Trim, ComputeMask

# ------------------------------------------------------------------------------------------------
# FMRI Clean
# ------------------------------------------------------------------------------------------------

# rs-fMRI preprocessing nodes
trim = pe.Node(Trim(), name="trim")

# slice-timing correction
params = pe.Node(STCParametersInterface(), name='stc_params')
params.inputs.time_repetition = 2
params.inputs.slice_mode = 'alt_inc'

gunzip = pe.Node(Gunzip(), name="gunzip")

stc = spm.SliceTiming()
stc.inputs.in_files = traits.Undefined
stc.inputs.out_prefix = 'stc'
slice_timing = pe.Node(stc, name='slice_timing')

wf.connect([
    # trim
    (datasource, trim, [("rest", "in_file")]),

    # slice time correction
    (trim, params, [("out_file", "in_files")]),

    # processing nodes
    (params, gunzip, [(("in_files", _pick_first), "in_file")]),
    (params, slice_timing, [
        (("slice_order", _sum_one_to_each), "slice_order"),
        (("ref_slice",   _sum_one), "ref_slice"),
        ("num_slices", "num_slices"),
        ("time_acquisition", "time_acquisition"),
        ("time_repetition", "time_repetition"),
    ]),

    (gunzip, slice_timing, [("out_file", "in_files")]),

])


# ------------------------------------------------------------------------------------------------
# FMRI Warp, Align, Filtering, Smoothing
# ------------------------------------------------------------------------------------------------
from nipype.interfaces.nipy import SpaceTimeRealigner
from nipype.algorithms.confounds import TSNR
from nipype.algorithms.rapidart import ArtifactDetect

from neuro_pypes.fmri.nuisance import rest_noise_filter_wf
from neuro_pypes.interfaces.nilearn import mean_img, smooth_img


realign = pe.Node(SpaceTimeRealigner(), name='realign')

# average
average = pe.Node(
    Function(
        function=mean_img,
        input_names=["in_file"],
        output_names=["out_file"],
        imports=['from neuro_pypes.interfaces.nilearn import ni2file']
    ),
    name='average_epi'
)

mean_gunzip = pe.Node(Gunzip(), name="mean_gunzip")

# co-registration nodes
coreg = spm.Coregister()
coreg.inputs.cost_function = "mi"
coreg.inputs.jobtype = 'estwrite'

coregister = pe.Node(coreg, name="coregister_fmri")
brain_sel = pe.Node(Select(index=[0, 1, 2]), name="brain_sel")

# brain mask made with EPI
epi_mask = pe.Node(ComputeMask(), name='epi_mask')

# brain mask made with the merge of the tissue segmentations
tissue_mask = pe.Node(fsl.MultiImageMaths(), name='tissue_mask')
tissue_mask.inputs.op_string = "-add %s -add %s -abs -kernel gauss 4 -dilM -ero -kernel gauss 1 -dilM -bin"
tissue_mask.inputs.out_file = "tissue_brain_mask.nii.gz"

# select tissues
gm_select = pe.Node(Select(index=[0]), name="gm_sel")
wmcsf_select = pe.Node(Select(index=[1, 2]), name="wmcsf_sel")

# noise filter
wm_select = pe.Node(Select(index=[1]), name="wm_sel")
csf_select = pe.Node(Select(index=[2]), name="csf_sel")


# anat to fMRI registration inputs
wf.connect([
#     (biascorr, coregister), [("output_image", "source")],
    (datasource, coregister, [("anat", "source")]),
    (segment, brain_sel, [("native_class_images", "inlist")]),
])


wf.connect([
    # motion correction
    (slice_timing, realign, [("timecorrected_files", "in_file")]),

    # coregistration target
    (realign, average, [("out_file", "in_file")]),
    (average, mean_gunzip, [("out_file", "in_file")]),
    (mean_gunzip, coregister, [("out_file", "target")]),

    # epi brain mask
    (average, epi_mask, [("out_file", "mean_volume")]),

    # coregistration
    (brain_sel, coregister, [(("out", flatten_list), "apply_to_files")]),

    # tissue brain mask
    (coregister, gm_select, [("coregistered_files", "inlist")]),
    (coregister, wmcsf_select, [("coregistered_files", "inlist")]),
    (gm_select, tissue_mask, [(("out", flatten_list), "in_file")]),
    (wmcsf_select, tissue_mask, [(("out", flatten_list), "operand_files")]),

    # nuisance correction
    (coregister, wm_select, [("coregistered_files", "inlist",)]),
    (coregister, csf_select, [("coregistered_files", "inlist",)]),
])


# ------------------------------------------------------------------------------------------------
# FMRI Noise removal
# ------------------------------------------------------------------------------------------------
from neuro_pypes.preproc import motion_regressors, extract_noise_components, create_regressors
from neuro_pypes.utils import selectindex, rename

# CompCor rsfMRI filters (at least compcor_csf should be True).
filters = {
    'compcor_csf': True,
    'compcor_wm': False,
    'gsr': False
}

# Compute TSNR on realigned data regressing polynomial up to order 2
tsnr = pe.Node(TSNR(regress_poly=2), name='tsnr')

# Use :class:`nipype.algorithms.rapidart` to determine which of the
# images in the functional series are outliers based on deviations in
# intensity or movement.
art = pe.Node(ArtifactDetect(), name="rapidart_artifacts")
# # Threshold to use to detect motion-related outliers when composite motion is being used
art.inputs.use_differences = [True, False]
art.inputs.use_norm = True
art.inputs.zintensity_threshold = 2
art.inputs.use_norm = True
art.inputs.norm_threshold = 1
art.inputs.mask_type = 'file'
art.inputs.parameter_source = 'NiPy'

# Compute motion regressors
motion_regs = pe.Node(
    Function(
        input_names=['motion_params', 'order', 'derivatives'],
        output_names=['out_files'],
        function=motion_regressors
    ),
    name='motion_regressors'
)
# motion regressors upto given order and derivative
# motion + d(motion)/dt + d2(motion)/dt2 (linear + quadratic)
motion_regs.inputs.order = 0
motion_regs.inputs.derivatives = 1

# Create a filter to remove motion and art confounds
motart_pars = pe.Node(
    Function(
        input_names=['motion_params', 'comp_norm', 'outliers', 'detrend_poly'],
        output_names=['out_files'],
        function=create_regressors
    ),
    name='motart_parameters'
)
# # number of polynomials to add to detrend
motart_pars.inputs.detrend_poly = 2

motion_filter = pe.Node(
    fsl.GLM(
        out_f_name='F_mcart.nii.gz',
        out_pf_name='pF_mcart.nii.gz',
        demean=True
    ),
    name='motion_filter'
)

# Noise confound regressors
compcor_pars = pe.Node(
    Function(
        input_names=['realigned_file', 'mask_file', 'num_components', 'extra_regressors'],
        output_names=['components_file'],
        function=extract_noise_components
    ),
    name='compcor_pars'
)
# Number of principal components to calculate when running CompCor. 5 or 6 is recommended.
compcor_pars.inputs.num_components = 6

compcor_filter = pe.Node(
    fsl.GLM(out_f_name='F.nii.gz', out_pf_name='pF.nii.gz', demean=True),
    name='compcor_filter'
)

# Global signal regression
gsr_pars = pe.Node(
    Function(
        input_names=['realigned_file', 'mask_file', 'num_components', 'extra_regressors'],
        output_names=['components_file'],
        function=extract_noise_components
    ),
    name='gsr_pars'
)
# Number of principal components to calculate when running Global Signal Regression. 1 is recommended.
gsr_pars.inputs.num_components: 1

gsr_filter = pe.Node(
    fsl.GLM(out_f_name='F_gsr.nii.gz', out_pf_name='pF_gsr.nii.gz', demean=True),
    name='gsr_filter'
)

wf.connect([
    # tsnr
    (realign, tsnr, [
        ("out_file", "in_file"),
    ]),

    # artifact detection
    (tissue_mask, art, [("out_file", "mask_file")]),
    (realign, art, [
        ("out_file", "realigned_files"),
        ("par_file", "realignment_parameters")
    ]),

    # calculte motion regressors
    (realign, motion_regs, [
        ("par_file", "motion_params")
    ]),

    # create motion and confound regressors parameters file
    (art, motart_pars, [
        ("norm_files", "comp_norm"),
        ("outlier_files", "outliers"),
    ]),
    (motion_regs, motart_pars, [
        ("out_files", "motion_params")
    ]),

    # motion filtering
    (realign, motion_filter, [
        ("out_file", "in_file"),
        (("out_file", rename, "_filtermotart"), "out_res_name"),
    ]),
    (motart_pars, motion_filter, [
        (("out_files", selectindex, 0), "design")
    ]),
])

wf.connect([
    # output
    (tsnr, datasink, [("tsnr_file", "rest.@tsnr")]),

    (motart_pars, datasink, [("out_files", "rest.@motion_regressors")]),
    (motion_filter, datasink, [("out_res", "rest.@motion_corrected")]),
    (art, datasink, [
        ("displacement_files", "rest.artifact_stats.@displacement"),
        ("intensity_files", "rest.artifact_stats.@intensity"),
        ("norm_files", "rest.artifact_stats.@norm"),
        ("outlier_files", "rest.artifact_stats.@outliers"),
        ("plot_files", "rest.artifact_stats.@plots"),
        ("statistic_files", "rest.artifact_stats.@stats"),
    ]),
])


last_filter = motion_filter

# compcor filter
if filters['compcor_csf'] or filters['compcor_wm']:
    wf.connect([
        # calculate compcor regressor and parameters file
        (motart_pars, compcor_pars, [(("out_files", selectindex, 0), "extra_regressors"), ]),
        (motion_filter, compcor_pars, [("out_res", "realigned_file"), ]),

        # the compcor filter
        (motion_filter, compcor_filter, [("out_res", "in_file"),
                                         (("out_res", rename, "_cleaned"), "out_res_name"),
                                         ]),
        (compcor_pars, compcor_filter, [("components_file", "design")]),
        (tissue_mask, compcor_filter, [("out_file", "mask")]),

        # output
        (compcor_pars, datasink, [("components_file", "rest.@compcor_regressors")]),
    ])
    last_filter = compcor_filter

# global signal regression
if filters['gsr']:
    wf.connect([
        # calculate gsr regressors parameters file
        (last_filter, gsr_pars, [("out_res", "realigned_file")]),
        (tissue_mask, gsr_pars, [("out_file", "mask_file")]),

        # the output file name
        (tissue_mask, gsr_filter, [("out_file", "mask")]),
        (last_filter, gsr_filter, [
            ("out_res", "in_file"),
            (("out_res", rename, "_gsr"), "out_res_name"),
        ]),
        (gsr_pars, gsr_filter, [("components_file", "design")]),

        # output
        (gsr_pars, datasink, [("components_file", "rest.@gsr_regressors")]),
    ])
    last_filter = gsr_filter

# connect the final nuisance correction output node
wf.connect([(last_filter, datasink, [("out_res", "rest.@nuis_corrected")]), ])

if filters['compcor_csf'] and filters['compcor_wm']:
    mask_merge = setup_node(Merge(2), name="mask_merge")
    wf.connect([
        ## the mask for the compcor filter
        (wm_select, mask_merge, [(("out", flatten_list), "in1")]),
        (csf_select, mask_merge, [(("out", flatten_list), "in2")]),
        (mask_merge, compcor_pars, [("out", "mask_file")]),
    ])

elif filters['compcor_csf']:
    wf.connect([
        ## the mask for the compcor filter
        (csf_select, compcor_pars, [(("out", flatten_list), "mask_file")]),
    ])

elif filters['compcor_wm']:
    wf.connect([
        ## the mask for the compcor filter
        (wm_select, compcor_pars, [(("out", flatten_list), "mask_file")]),
    ])


# In[ ]:


from neuro_pypes.fmri.filter import bandpass_filter
from neuro_pypes.interfaces.nilearn import smooth_img

# bandpass filtering
bandpass = pe.Node(
    Function(
        input_names=['files', 'lowpass_freq', 'highpass_freq', 'tr'],
        output_names=['out_files'],
        function=bandpass_filter
    ),
    name='bandpass'
)
bandpass.inputs.lowpass_freq = 0.1
bandpass.inputs.highpass_freq = 0.01

# smooth
smooth = pe.Node(
    Function(
        function=smooth_img,
        input_names=["in_file", "fwhm"],
        output_names=["out_file"],
        imports=['from neuro_pypes.interfaces.nilearn import ni2file']
    ),
    name="smooth"
)
smooth.inputs.fwhm = fmri_smoothing_kernel_fwhm
smooth.inputs.out_file = "smooth_std_{}.nii.gz".format(wf_name)


wf.connect([
    # temporal filtering
    (last_filter, bandpass, [("out_res", "files")]),

    # (realign,     bandpass,    [("out_file", "files")]),
    (params, bandpass, [("time_repetition", "tr")]),
    (bandpass, smooth, [("out_files", "in_file")]),

    # output
    (epi_mask, datasink, [("brain_mask", "rest.@epi_brain_mask")]),
    (tissue_mask, datasink, [("out_file", "rest.@tissues_brain_mask")]),
    (realign, datasink, [
        ("out_file", "rest.@realigned"),
        ("par_file", "rest.@motion_params"),
    ]),
    (coregister, datasink, [
        ("coregistered_files", "rest.@tissues"),
        ("coregistered_source", "rest.@anat"),
    ]),
    (average, datasink, [("out_file", "rest.@avg_epi")]),
    (bandpass, datasink, [("out_files", "rest.@time_filtered")]),
    (smooth, datasink, [("out_file", "rest.@smooth")]),
])


if __name__ == '__main__':
    n_cpus = 1

    if n_cpus > 1:
        wf.run(plugin=plugin, plugin_args={"n_procs": n_cpus})
    else:
        wf.run(plugin=None)
