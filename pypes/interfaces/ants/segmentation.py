# -*- coding: utf-8 -*-
"""
Nipype interface for ANTs.

The implementation of this interface is motivated by this study:
http://dx.doi.org/10.1016/j.nicl.2016.05.017
"""

from nipype.interfaces.base import TraitedSpec, File, traits, isdefined
from nipype.interfaces.ants.base import ANTSCommand, ANTSCommandInputSpec
from nipype.utils.filemanip import split_filename
from nipype.external.due import BibTeX


class KellyKapowskiInputSpec(ANTSCommandInputSpec):
    dimension = traits.Enum(3, 2, argstr='--image-dimensionality %d', usedefault=True,
                            desc='image dimension (2 or 3)')

    segmentation_image = File(exists=True, argstr='--segmentation-image "%s"', mandatory=True,
                              desc="A segmentation image must be supplied labeling the gray and white matters.\n"
                                   "Default values = 2 and 3, respectively.",)

    gray_matter_label = traits.Int(2, usedefault=True,
                                   desc="The label value for the gray matter label in the segmentation_image.")

    white_matter_label = traits.Int(3, usedefault=True,
                                    desc="The label value for the white matter label in the segmentation_image.")

    gray_matter_prob_image = File(exists=True, argstr='--gray-matter-probability-image "%s"',
                                  desc="In addition to the segmentation image, a gray matter probability image can be\n"
                                       "used. If no such image is supplied, one is created using the segmentation image\n"
                                       "and a variance of 1.0 mm.")

    white_matter_prob_image = File(exists=True, argstr='--white-matter-probability-image "%s"',
                                   desc="In addition to the segmentation image, a white matter probability image can be\n"
                                       "used. If no such image is supplied, one is created using the segmentation image\n"
                                       "and a variance of 1.0 mm.")

    convergence = traits.Str(default="[50,0.001,10]", argstr='--convergence "%s"', usedefault=True,
                             desc="Convergence is determined by fitting a line to the normalized energy profile of\n"
                                  "the last N iterations (where N is specified by the window size) and determining\n"
                                  "the slope which is then compared with the convergence threshold.",)

    thickness_prior_estimate = traits.Float(10, usedefault=True, argstr="--thickness-prior-estimate %f",
                                            desc="Provides a prior constraint on the final thickness measurement in mm.")

    thickness_prior_image = File(exists=True, argstr='--thickness-prior-image "%s"',
                                 desc="An image containing spatially varying prior thickness values.")

    gradient_step = traits.Float(0.025, usedefault=True, argstr="--gradient-step %f",
                                 desc="Gradient step size for the optimization.")

    smoothing_variance = traits.Float(1.0, argstr="--smoothing-variance %f",
                                      desc="Defines the Gaussian smoothing of the hit and total images.")

    smoothing_velocity_field = traits.Float(1.5, argstr="--smoothing-velocity-field-parameter %f",
                                            desc="Defines the Gaussian smoothing of the velocity field (default = 1.5).\n"
                                            "If the b-spline smoothing option is chosen, then this defines the \n"
                                            "isotropic mesh spacing for the smoothing spline (default = 15).")

    use_bspline_smoothing = traits.Bool(argstr="--use-bspline-smoothing 1",
                                        desc="Sets the option for B-spline smoothing of the velocity field.")

    number_integration_points = traits.Int(10, argstr="--number-of-integration-points %d",
                                           desc="Number of compositions of the diffeomorphism per iteration.")

    max_invert_displacement_field_iters = traits.Int(20, argstr="--maximum-number-of-invert-displacement-field-iterations %d",
                                                     desc="Maximum number of iterations for estimating the invert \n"
                                                          "displacement field.")

    cortical_thickness = File(argstr='--output "%s"', keep_extension=True,
                              name_source=["segmentation_image"], name_template='%s_cortical_thickness',
                              desc='Filename for the cortical thickness.', hash_files=False)

    warped_white_matter = File(name_source=["segmentation_image"], keep_extension=True,
                               name_template='%s_warped_white_matter',
                               desc='Filename for the warped white matter file.', hash_files=False)


class KellyKapowskiOutputSpec(TraitedSpec):
    cortical_thickness = File(desc="A thickness map defined in the segmented gray matter.")
    warped_white_matter = File(desc="A warped white matter image.")


class KellyKapowski(ANTSCommand):
    """ Nipype Interface to ANTs' KellyKapowski, also known as DiReCT.

    DiReCT is a registration based estimate of cortical thickness. It was published
    in S. R. Das, B. B. Avants, M. Grossman, and J. C. Gee, Registration based
    cortical thickness measurement, Neuroimage 2009, 45:867--879.

    Examples
    --------
    >>> from nipype.interfaces.ants.segmentation import KellyKapowski
    >>> kk = KellyKapowski()
    >>> kk.inputs.dimension = 3
    >>> kk.inputs.segmentation_image = "segmentation0.nii.gz"
    >>> kk.inputs.convergence = "[45,0.0,10]"
    >>> kk.inputs.gradient_step = 0.025
    >>> kk.inputs.smoothing_variance = 1.0
    >>> kk.inputs.smoothing_velocity_field = 1.5
    >>> #kk.inputs.use_bspline_smoothing = False
    >>> kk.inputs.number_integration_points = 10
    >>> kk.inputs.thickness_prior_estimate = 10
    >>> kk.cmdline # doctest: +ALLOW_UNICODE
    u'KellyKapowski --convergence "[45,0.0,10]" \
--output "[segmentation0_cortical_thickness.nii.gz,segmentation0_warped_white_matter.nii.gz]" \
--image-dimensionality 3 --gradient-step 0.025000 --number-of-integration-points 10 \
--segmentation-image "[segmentation0.nii.gz,2,3]" --smoothing-variance 1.000000 \
--smoothing-velocity-field-parameter 1.500000 --thickness-prior-estimate 10.000000'

    """
    _cmd = "KellyKapowski"
    input_spec = KellyKapowskiInputSpec
    output_spec = KellyKapowskiOutputSpec

    references_ = [{'entry': BibTeX("@book{Das2009867,"
                                    "author={Sandhitsu R. Das and Brian B. Avants and Murray Grossman and James C. Gee},"
                                    "title={Registration based cortical thickness measurement.},"
                                    "journal={NeuroImage},"
                                    "volume={45},"
                                    "number={37},"
                                    "pages={867--879},"
                                    "year={2009},"
                                    "issn={1053-8119},"
                                    "url={http://www.sciencedirect.com/science/article/pii/S1053811908012780},"
                                    "doi={http://dx.doi.org/10.1016/j.neuroimage.2008.12.016}"
                                    "}"),
                    'description': 'The details on the implementation of DiReCT.',
                    'tags': ['implementation'],
                    }]

    def _parse_inputs(self, skip=None):
        if skip is None:
            skip = []
        skip += ['warped_white_matter', 'gray_matter_label', 'white_matter_label']
        return super(KellyKapowski, self)._parse_inputs(skip=skip)

    def _gen_filename(self, name):
        if name == 'cortical_thickness':
            output = self.inputs.cortical_thickness
            if not isdefined(output):
                _, name, ext = split_filename(self.inputs.segmentation_image)
                output = name + '_cortical_thickness' + ext
            return output

        if name == 'warped_white_matter':
            output = self.inputs.warped_white_matter
            if not isdefined(output):
                _, name, ext = split_filename(self.inputs.segmentation_image)
                output = name + '_warped_white_matter' + ext
            return output

        return None

    def _format_arg(self, opt, spec, val):
        if opt == "segmentation_image":
            newval = '[{0},{1},{2}]'.format(self.inputs.segmentation_image,
                                            self.inputs.gray_matter_label,
                                            self.inputs.white_matter_label)
            return spec.argstr % newval

        if opt == "cortical_thickness":
            ct = self._gen_filename("cortical_thickness")
            wm = self._gen_filename("warped_white_matter")
            newval = '[{},{}]'.format(ct, wm)
            return spec.argstr % newval

        return super(KellyKapowski, self)._format_arg(opt, spec, val)
