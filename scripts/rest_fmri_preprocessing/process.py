import os

import nipype.pipeline.engine as pe
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import IdentityInterface

from neuro_pypes.crumb import DataCrumb

# from nipype.algorithms.misc import Gunzip
# from nipype.interfaces import fsl
# from nipype.interfaces.nipy.preprocess import Trim, ComputeMask
# from nipype.interfaces.utility import Function, Select, IdentityInterface

# from neuro_pypes._utils import format_pair_list, flatten_list
# from neuro_pypes.config import setup_node, get_config_setting
# from neuro_pypes.fmri.filter import bandpass_filter
# from neuro_pypes.fmri.nuisance import rest_noise_filter_wf
# from neuro_pypes.interfaces.nilearn import mean_img, smooth_img
# from neuro_pypes.preproc import (
#     auto_spm_slicetime,
#     nipy_motion_correction,
#     spm_coregister
# )
# from neuro_pypes.utils import (remove_ext,
#                                extend_trait_list,
#                                get_input_node,
#                                get_interface_node,
#                                get_datasink,
#                                get_input_file_name,
#                                extension_duplicates)

# ------------------------------------------------------------------------------------------------
# GLOBAL VARIABLES
# ------------------------------------------------------------------------------------------------

wf_name = 'rest_fmri_preprocess'

# STDB_DIR = os.path.expanduser('~/projects/neuro/std_brains')
# SPM_DIR = os.path.expanduser('~/Software/spm_mcr')
# BASE_DIR = os.path.expanduser('~/Data/nuk/petnet')
# data_dir = os.path.join(BASE_DIR, 'raw')
# cache_dir = os.path.join(BASE_DIR, 'wd')
# output_dir = os.path.join(BASE_DIR, 'out')
# plugin = None
# n_cpus = 5

# HAMM_DIR = os.path.join(STDB_DIR, 'atlases', 'hammers')
# HAMM_MNI = os.path.join(HAMM_DIR, 'Hammers_mith_atlas_n30r83_SPM5.nii.gz')
# HAMM_LABELS = os.path.join(HAMM_DIR, 'labels.txt')

# SPM_CANONICAL_BRAIN_2MM = os.path.join(STDB_DIR, 'templates', 'spm_canonical', 'single_subj_T1_brain.nii.gz')

# # template files
# PET_MNI  = os.path.join(STDB_DIR, 'templates', 'spm_canonical', 'pet.nii')
# MNI_MASK = os.path.join(STDB_DIR, 'templates', 'avg152T1_brain.nii.gz')

# # data input/output os.path.dirname(__file__) means alongside with this .py file
# #settings_file = os.path.join(os.path.dirname(__file__), 'pypes_config.yml')
# settings_file = ('/home/iripp/projects/alex/nuk_experiments/MRPET_15/Preproc_30_60pi_pet_recon_NEW/pypes_config.yml')
# #nipype_cfg_file = os.path.join(os.path.dirname(__file__), 'nipype.cfg')
# nipype_cfg_file = ('/home/iripp/projects/alex/nuk_experiments/MRPET_15/Preproc_30_60pi_pet_recon_NEW/nipype.cfg')



# mrpet15_preproc_wf_2 = dict([
#     ("spm_anat_preproc", attach_spm_anat_preprocessing),
#     ("spm_pet_preproc", attach_spm_pet_preprocessing),
#     ("spm_mrpet_preproc", attach_spm_mrpet_preprocessing),
#     ("spm_pet_grouptemplate", attach_spm_pet_grouptemplate),
# ])

# data_path = os.path.join(os.path.expanduser(data_dir), '{session}', '{subject_id}', '{scan}', '{image}')
# data_crumb = Crumb(data_path, ignore_list=['.*'])


# crumb_modalities = {
#     'anat': [('scan', 'T1'), ('image', 'Head_MPRAGE_highContrast.nii.gz')],
#     'pet': [('scan', 'PET_recon_first_30min'), ('image', 'pet_recon.nii.gz')],
# }


# ------------------------------------------------------------------------------------------------
# DATA INPUT AND SINK
# ------------------------------------------------------------------------------------------------

wf = pe.Workflow(name=wf_name, base_dir=work_dir)

# datasink
datasink = pe.Node(
    DataSink(parameterization=False, base_directory=output_dir, ),
    name="datasink"
)

# input workflow
# (work_dir, data_crumb, crumb_arg_values, files_crumb_args, wf_name="input_files"):
select_files = pe.Node(
    DataCrumb(crumb=data_crumb, templates=file_templates, raise_on_empty=False),
    name='selectfiles'
)

# basic file name substitutions for the datasink
undef_args = select_files.interface._infields
substitutions = [(name, "") for name in undef_args]
substitutions.append(("__", "_"))

datasink.inputs.substitutions = extend_trait_list(datasink.inputs.substitutions,
                                                    substitutions)

# Infosource - the information source that iterates over crumb values map from the filesystem
infosource = pe.Node(interface=IdentityInterface(fields=undef_args), name="infosrc")
infosource.iterables = list(valuesmap_to_dict(joint_value_map(data_crumb, undef_args)).items())
infosource.synchronize = True

# connect the input_wf to the datasink
joinpath = pe.Node(joinstrings(len(undef_args)), name='joinpath')

# Connect the infosrc node to the datasink
input_joins = [(name, 'arg{}'.format(arg_no + 1))
                for arg_no, name in enumerate(undef_args)]

wf.connect([
    (infosource, select_files, [(field, field) for field in undef_args]),
    (select_files, joinpath, input_joins),
    (joinpath, datasink, [("out", "container")]),
],
)

# ------------------------------------------------------------------------------------------------
# ANAT
# ------------------------------------------------------------------------------------------------

    # input node
    anat_input = pe.Node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
                         name="anat_input")

    # atlas registration
    if do_atlas and not isdefined(anat_input.inputs.atlas_file):
        anat_input.inputs.set(atlas_file=atlas_file)

    # T1 preprocessing nodes
    biascor     = setup_node(biasfield_correct(),      name="bias_correction")
    gunzip_anat = setup_node(Gunzip(),                 name="gunzip_anat")
    segment     = setup_node(spm_segment(),            name="new_segment")
    warp_anat   = setup_node(spm_apply_deformations(), name="warp_anat")

    tpm_bbox = setup_node(Function(function=get_bounding_box,
                                   input_names=["in_file"],
                                   output_names=["bbox"]),
                          name="tpm_bbox")
    tpm_bbox.inputs.in_file = spm_tpm_priors_path()

    # calculate brain mask from tissue maps
    tissues = setup_node(IdentityInterface(fields=["gm", "wm", "csf"], mandatory_inputs=True),
                         name="tissues")

    brain_mask = setup_node(Function(function=math_img,
                                     input_names=["formula", "out_file", "gm", "wm", "csf"],
                                     output_names=["out_file"],
                                     imports=['from neuro_pypes.interfaces.nilearn import ni2file']),
                            name='brain_mask')
    brain_mask.inputs.out_file = "tissues_brain_mask.nii.gz"
    brain_mask.inputs.formula  = "np.abs(gm + wm + csf) > 0"

    # output node
    anat_output = pe.Node(IdentityInterface(fields=out_fields), name="anat_output")

    # Connect the nodes
    wf.connect([
                # input to biasfieldcorrection
                (anat_input,   biascor    , [("in_file",      "input_image")]),

                # new segment
                (biascor,      gunzip_anat, [("output_image", "in_file")]),
                (gunzip_anat,  segment,     [("out_file",     "channel_files")]),

                # Normalize12
                (segment,   warp_anat,  [("forward_deformation_field", "deformation_file")]),
                (segment,   warp_anat,  [("bias_corrected_images",     "apply_to_files")]),
                (tpm_bbox,  warp_anat,  [("bbox",                      "write_bounding_box")]),

                # brain mask from tissues
                (segment, tissues,  [(("native_class_images", selectindex, 0), "gm"),
                                     (("native_class_images", selectindex, 1), "wm"),
                                     (("native_class_images", selectindex, 2), "csf"),
                                    ]),

                (tissues, brain_mask, [("gm", "gm"), ("wm", "wm"), ("csf", "csf"),]),

                # output
                (warp_anat, anat_output, [("normalized_files",           "anat_mni")]),
                (segment,   anat_output, [("modulated_class_images",     "tissues_warped"),
                                          ("native_class_images",        "tissues_native"),
                                          ("transformation_mat",         "affine_transform"),
                                          ("forward_deformation_field",  "warp_forward"),
                                          ("inverse_deformation_field",  "warp_inverse"),
                                          ("bias_corrected_images",      "anat_biascorr")]),
                (brain_mask, anat_output, [("out_file",                  "brain_mask")]),
              ])

    # atlas warping nodes
    if do_atlas:
        gunzip_atlas = pe.Node(Gunzip(), name="gunzip_atlas")
        warp_atlas   = setup_node(spm_apply_deformations(), name="warp_atlas")
        anat_bbox    = setup_node(Function(function=get_bounding_box,
                                           input_names=["in_file"],
                                           output_names=["bbox"]),
                                  name="anat_bbox")

        # set the warping interpolation to nearest neighbour.
        warp_atlas.inputs.write_interp = 0

        # connect the atlas registration nodes
        wf.connect([
                    (anat_input,    gunzip_atlas, [("atlas_file",                 "in_file")]),
                    (gunzip_anat,   anat_bbox,    [("out_file",                   "in_file")]),
                    (gunzip_atlas,  warp_atlas,   [("out_file",                   "apply_to_files")]),
                    (segment,       warp_atlas,   [("inverse_deformation_field",  "deformation_file")]),
                    (anat_bbox,     warp_atlas,   [("bbox",                       "write_bounding_box")]),
                    (warp_atlas,    anat_output,  [("normalized_files",           "atlas_anat")]),
                  ])

# Create the workflow object
wf = pe.Workflow(name=wf_name)

wf.connect([(in_files, anat_wf,  [("anat", "anat_input.in_file")]),
            (anat_wf,  datasink, [
                ("anat_output.anat_mni",         "anat.@mni"),
                ("anat_output.tissues_warped",   "anat.tissues.warped"),
                ("anat_output.tissues_native",   "anat.tissues.native"),
                ("anat_output.affine_transform", "anat.transform.@linear"),
                ("anat_output.warp_forward",     "anat.transform.@forward"),
                ("anat_output.warp_inverse",     "anat.transform.@inverse"),
                ("anat_output.anat_biascorr",    "anat.@biascor"),
                ("anat_output.brain_mask",       "anat.@brain_mask"),
            ]),
        ])

# check optional outputs
if do_atlas:
    wf.connect([(anat_wf, datasink, [("anat_output.atlas_anat", "anat.@atlas")]),])

do_cortical_thickness = get_config_setting('anat_preproc.do_cortical_thickness', False)
if do_cortical_thickness:
    wf.connect([(anat_wf, datasink, [("anat_output.cortical_thickness",  "anat.@cortical_thickness"),
                                            ("anat_output.warped_white_matter", "anat.@warped_white_matter"),
                                            ]),
                    ])

# ------------------------------------------------------------------------------------------------
# FMRI
# ------------------------------------------------------------------------------------------------

# # specify input and output fields
# in_fields = [
#     "in_file",
#     "anat",
#     "atlas_anat",
#     "coreg_target",
#     "tissues",
#     "lowpass_freq",
#     "highpass_freq",
# ]

# out_fields = [
#     "motion_corrected",
#     "motion_params",
#     "tissues",
#     "anat",
#     "avg_epi",
#     "time_filtered",
#     "smooth",
#     "tsnr_file",
#     "epi_brain_mask",
#     "tissues_brain_mask",
#     "motion_regressors",
#     "compcor_regressors",
#     "gsr_regressors",
#     "nuis_corrected",
#     "art_displacement_files",
#     "art_intensity_files",
#     "art_norm_files",
#     "art_outlier_files",
#     "art_plot_files",
#     "art_statistic_files",
# ]

# # input identities
# rest_input = setup_node(IdentityInterface(fields=in_fields, mandatory_inputs=True),
#                         name="rest_input")

# # rs-fMRI preprocessing nodes
# trim = setup_node(Trim(), name="trim")

# stc_wf = auto_spm_slicetime()
# realign = setup_node(nipy_motion_correction(), name='realign')

# # average
# average = setup_node(
#     Function(
#         function=mean_img,
#         input_names=["in_file"],
#         output_names=["out_file"],
#         imports=['from neuro_pypes.interfaces.nilearn import ni2file']
#     ),
#     name='average_epi'
# )

# mean_gunzip = setup_node(Gunzip(), name="mean_gunzip")

# # co-registration nodes
# coreg = setup_node(spm_coregister(cost_function="mi"), name="coreg_fmri")
# brain_sel = setup_node(Select(index=[0, 1, 2]), name="brain_sel")

# # brain mask made with EPI
# epi_mask = setup_node(ComputeMask(), name='epi_mask')

# # brain mask made with the merge of the tissue segmentations
# tissue_mask = setup_node(fsl.MultiImageMaths(), name='tissue_mask')
# tissue_mask.inputs.op_string = "-add %s -add %s -abs -kernel gauss 4 -dilM -ero -kernel gauss 1 -dilM -bin"
# tissue_mask.inputs.out_file = "tissue_brain_mask.nii.gz"

# # select tissues
# gm_select = setup_node(Select(index=[0]), name="gm_sel")
# wmcsf_select = setup_node(Select(index=[1, 2]), name="wmcsf_sel")

# # noise filter
# noise_wf = rest_noise_filter_wf()
# wm_select = setup_node(Select(index=[1]), name="wm_sel")
# csf_select = setup_node(Select(index=[2]), name="csf_sel")

# # bandpass filtering
# bandpass = setup_node(
#     Function(
#         input_names=['files', 'lowpass_freq', 'highpass_freq', 'tr'],
#         output_names=['out_files'],
#         function=bandpass_filter
#     ),
#     name='bandpass'
# )

# # smooth
# smooth = setup_node(
#     Function(
#         function=smooth_img,
#         input_names=["in_file", "fwhm"],
#         output_names=["out_file"],
#         imports=['from neuro_pypes.interfaces.nilearn import ni2file']
#     ),
#     name="smooth"
# )
# smooth.inputs.fwhm = get_config_setting('fmri_smooth.fwhm', default=8)
# smooth.inputs.out_file = "smooth_std_{}.nii.gz".format(wf_name)

# # output identities
# rest_output = setup_node(IdentityInterface(fields=out_fields), name="rest_output")

# # Connect the nodes
# wf.connect([
#     # trim
#     (rest_input, trim, [("in_file", "in_file")]),

#     # slice time correction
#     (trim, stc_wf, [("out_file", "stc_input.in_file")]),

#     # motion correction
#     (stc_wf, realign, [("stc_output.timecorrected_files", "in_file")]),

#     # coregistration target
#     (realign, average, [("out_file", "in_file")]),
#     (average, mean_gunzip, [("out_file", "in_file")]),
#     (mean_gunzip, coreg, [("out_file", "target")]),

#     # epi brain mask
#     (average, epi_mask, [("out_file", "mean_volume")]),

#     # coregistration
#     (rest_input, coreg, [("anat", "source")]),
#     (rest_input, brain_sel, [("tissues", "inlist")]),
#     (brain_sel, coreg, [(("out", flatten_list), "apply_to_files")]),

#     # tissue brain mask
#     (coreg, gm_select, [("coregistered_files", "inlist")]),
#     (coreg, wmcsf_select, [("coregistered_files", "inlist")]),
#     (gm_select, tissue_mask, [(("out", flatten_list), "in_file")]),
#     (wmcsf_select, tissue_mask, [(("out", flatten_list), "operand_files")]),

#     # nuisance correction
#     (coreg, wm_select, [("coregistered_files", "inlist",)]),
#     (coreg, csf_select, [("coregistered_files", "inlist",)]),
#     (realign, noise_wf, [("out_file", "rest_noise_input.in_file",)]),
#     (tissue_mask, noise_wf, [("out_file", "rest_noise_input.brain_mask")]),
#     (wm_select, noise_wf, [(("out", flatten_list), "rest_noise_input.wm_mask")]),
#     (csf_select, noise_wf, [(("out", flatten_list), "rest_noise_input.csf_mask")]),

#     (realign, noise_wf, [("par_file", "rest_noise_input.motion_params",)]),

#     # temporal filtering
#     (noise_wf, bandpass, [("rest_noise_output.nuis_corrected", "files")]),
#     # (realign,     bandpass,    [("out_file", "files")]),
#     (stc_wf, bandpass, [("stc_output.time_repetition", "tr")]),
#     (rest_input, bandpass, [
#         ("lowpass_freq", "lowpass_freq"),
#         ("highpass_freq", "highpass_freq"),
#     ]),
#     (bandpass, smooth, [("out_files", "in_file")]),

#     # output
#     (epi_mask, rest_output, [("brain_mask", "epi_brain_mask")]),
#     (tissue_mask, rest_output, [("out_file", "tissues_brain_mask")]),
#     (realign, rest_output, [
#         ("out_file", "motion_corrected"),
#         ("par_file", "motion_params"),
#     ]),
#     (coreg, rest_output, [
#         ("coregistered_files", "tissues"),
#         ("coregistered_source", "anat"),
#     ]),
#     (noise_wf, rest_output, [
#         ("rest_noise_output.motion_regressors", "motion_regressors"),
#         ("rest_noise_output.compcor_regressors", "compcor_regressors"),
#         ("rest_noise_output.gsr_regressors", "gsr_regressors"),
#         ("rest_noise_output.nuis_corrected", "nuis_corrected"),
#         ("rest_noise_output.tsnr_file", "tsnr_file"),
#         ("rest_noise_output.art_displacement_files", "art_displacement_files"),
#         ("rest_noise_output.art_intensity_files", "art_intensity_files"),
#         ("rest_noise_output.art_norm_files", "art_norm_files"),
#         ("rest_noise_output.art_outlier_files", "art_outlier_files"),
#         ("rest_noise_output.art_plot_files", "art_plot_files"),
#         ("rest_noise_output.art_statistic_files", "art_statistic_files"),
#     ]),
#     (average, rest_output, [("out_file", "avg_epi")]),
#     (bandpass, rest_output, [("out_files", "time_filtered")]),
#     (smooth, rest_output, [("out_file", "smooth")]),
# ])
