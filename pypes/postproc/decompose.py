# -*- coding: utf-8 -*-
"""
Nipype interfaces to perform signal decomposition.
"""
import os.path as op

import numpy as np
from   nilearn.decomposition     import CanICA
import nipype.pipeline.engine    as pe
from   nipype.interfaces         import io
from   nipype.interfaces.utility import IdentityInterface, Function
from   nipype.interfaces.base    import (BaseInterface,
                                         TraitedSpec,
                                         InputMultiPath,
                                         OutputMultiPath,
                                         BaseInterfaceInputSpec,
                                         traits,)

from   ..nilearn import concat_imgs
from   ..config import setup_node
from   ..utils  import (get_trait_value,
                        get_datasink,)


class CanICAInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(traits.File(desc="NifTI image file(s) from where to extract the data. \n"
                                               "If more than one, all should be spatially normalized.",
                                          exists=True, mandatory=True))
    mask = traits.File(desc="Mask to be used on data. If an instance of masker is passed, then its mask will be used.\n"
                            " If no mask is given, it will be computed automatically by a MultiNiftiMasker \n"
                            "with default parameters.",
                       exists=True)
    confounds = traits.File(desc="CSV file path."
                                 "This parameter is passed to nilearn.signal.clean. "
                                 "Please see the related documentation for details.",
                            exists=True)
    n_components = traits.Int(desc="Number of components to extract.",
                              default_value=20)
    do_cca = traits.Bool(desc="Indicate if a Canonical Correlation Analysis must be run after the PCA.",
                         default_value=True,)
    threshold = traits.Either(traits.Bool, traits.Float,
                              desc="If None, no thresholding is applied.\n"
                                   "If ‘auto’, then we apply a thresholding that will keep the n_voxels, \n"
                                   "more intense voxels across all the maps, n_voxels being the number of voxels \n"
                                   "in a brain volume. A float value indicates the ratio of voxels to keep "
                                   "(2. means that the maps will together have 2 x n_voxels non-zero voxels ).",)
    random_state = traits.Int(desc="Pseudo number generator state used for random sampling.",)
    smoothing_fwhm = traits.Float(desc="If smoothing_fwhm is defined, it gives the full-width half maximum in "
                                       "millimeters of the spatial smoothing to apply to the signal.",
                                 )
    standardize = traits.Bool(desc="If standardize is True, the time-series are centered and normed: "
                                   "their mean is put to 0 and their variance to 1 in the time dimension.",
                              default_value=True)


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

        n_components      = get_trait_value(self.inputs, 'n_components')
        do_cca            = get_trait_value(self.inputs, 'do_cca')
        smoothing_fwhm    = get_trait_value(self.inputs, 'smoothing_fwhm', default=None)
        standardize       = get_trait_value(self.inputs, 'standardize',    default=None)
        threshold         = get_trait_value(self.inputs, 'threshold',      default=None)
        random_state      = get_trait_value(self.inputs, 'random_state',   default=None)
        self._confounds   = get_trait_value(self.inputs, 'confounds',)

        # init the estimator
        self._estimator_name = 'canica'
        self._estimator = CanICA(n_components=n_components,
                                 threshold=threshold,
                                 random_state=random_state,
                                 standardize=standardize,
                                 smoothing_fwhm=smoothing_fwhm,
                                 do_cca=do_cca,
                                 verbose=1)

        # set output file names
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


def attach_canica(main_wf, wf_name="canica", **kwargs):
    """ Attach a nilearn CanICA interface to `main_wf`.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    kwargs: dict[str]->str
        input_node: str
            Name of the input node from where to connect the source `input_connect`.

        input_connection: str
            Name of the connection to obtain the source files.

    Nipype Inputs for `main_wf`
    ---------------------------
    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    srcwf_name   = kwargs['input_node']
    srcconn_name = kwargs['input_connection']

    src_wf   = main_wf.get_node(srcwf_name)
    datasink = get_datasink(main_wf, name='datasink')

    base_outdir  = datasink.inputs.base_directory
    ica_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,),
                           name="{}_datasink".format(wf_name))

    # the list of the raw pet subjects
    ica_subjs = pe.JoinNode(interface=IdentityInterface(fields=["ica_subjs"]),
                            joinsource="infosrc",
                            joinfield="ica_subjs",
                            name="ica_subjs")

    # warp each subject to the group template
    canica = setup_node(CanICAInterface(), name="{}_ica".format(wf_name),)

    # Connect the nodes
    main_wf.connect([
                     # file list input
                     (src_wf,       ica_subjs, [(srcconn_name,   "ica_subjs")]),

                     # canica
                     (ica_subjs,  canica,      [("ica_subjs",    "in_files")]),

                     # canica output
                     (canica, ica_datasink, [("components", "canica.@components")]),
                     (canica, ica_datasink, [("loadings",   "canica.@loadings")]),
                     (canica, ica_datasink, [("score",      "canica.@score")]),
                   ])
    return main_wf


def attach_concat_canica(main_wf, wf_name="canica", **kwargs):
    """ Attach a Concat and a nilearn CanICA interface to `main_wf`.

    The Concat node will merge all the files together in one 4D volume before delivering it to CanICA.

    Parameters
    ----------
    main_wf: nipype Workflow

    wf_name: str
        Name of the preprocessing workflow

    kwargs: dict[str]->str
        input_node: str
            Name of the input node from where to connect the source `input_connect`.

        input_connection: str
            Name of the connection to obtain the source files.

    Nipype Inputs for `main_wf`
    ---------------------------
    datasink: nipype Node

    Returns
    -------
    main_wf: nipype Workflow
    """
    # Dependency workflows
    srcwf_name   = kwargs['input_node']
    srcconn_name = kwargs['input_connection']

    src_wf   = main_wf.get_node(srcwf_name)
    datasink = get_datasink(main_wf, name='datasink')

    base_outdir  = datasink.inputs.base_directory
    ica_datasink = pe.Node(io.DataSink(parameterization=False,
                                       base_directory=base_outdir,),
                           name="ica_datasink".format(wf_name))
    ica_datasink.inputs.container = 'ica_{}'.format(wf_name)

    # the list of the raw pet subjects
    ica_subjs = pe.JoinNode(interface=IdentityInterface(fields=["ica_subjs"]),
                            joinsource="infosrc",
                            joinfield="ica_subjs",
                            name="ica_subjs")

    # concat images
    concat = setup_node(Function(function=concat_imgs,
                                 input_names=["in_files", "out_file"],
                                 output_names=["out_file"],),
                        name="concat")
    concat.inputs.out_file = 'concat_img.nii.gz'

    makelist = lambda x: [x]

    # warp each subject to the group template
    canica = setup_node(CanICAInterface(), name="{}_ica".format(wf_name),)

    # Connect the nodes
    main_wf.connect([
                     # file list input
                     (src_wf, ica_subjs, [(srcconn_name, "ica_subjs")]),

                     # concat images
                     (ica_subjs, concat, [("ica_subjs", "in_files")]),

                     # canica
                     (concat, canica, [(("out_file", makelist), "in_files")]),

                     # canica output
                     (canica, ica_datasink, [("components", "@components"),
                                             ("loadings",   "@loadings"),
                                             ("score",      "@score"),
                                            ]),
                   ])
    return main_wf