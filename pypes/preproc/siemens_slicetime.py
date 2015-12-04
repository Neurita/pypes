# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""SPM wrappers for slice timing correction of Siemens acquisitions data
"""
import os.path as op
import nibabel as nib

# Nipype imports
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, isdefined,
                                    traits, InputMultiPath, File)
from nipype.interfaces.spm.base import (SPMCommand, scans_for_fname,
                                        func_is_3d, Info,
                                        scans_for_fnames, SPMCommandInputSpec)
from nipype.utils.filemanip import (fname_presuffix, filename_to_list,
                                    list_to_filename, split_filename)

from .._utils import check_equal


class SiemensSliceTimingInputSpec(SPMCommandInputSpec):
    in_files = InputMultiPath(traits.Either(traits.List(File(exists=True)),
                                            File(exists=True)), field='scans',
                              desc='list of filenames to apply slice timing',
                              mandatory=True, copyfile=False)
    num_slices = traits.Int(field='nslices',
                            desc='number of slices in a volume')
    time_repetition = traits.Float(field='tr',
                                   desc=('time between volume acquisitions'
                                         '(start to start time)'),)
    time_acquisition = traits.Float(field='ta',
                                    desc=('time of volume acquisition. usually'
                                          'calculated as TR-(TR/num_slices)'),)
    slice_order = traits.List(field='so',
                              desc='1-based order in which slices are acquired',)
    slice_mode = traits.Enum('unknown', 'seq_inc', 'seq_dec', 'alt_inc', 'alt_dec', 'alt_inc2', 'alt_dec2',
                              field=None, usedefault=True,
                              desc='Description of the slice ordering method.',)
    ref_slice = traits.Int(field='refslice',
                           desc='1-based Number of the reference slice',)
    out_prefix = traits.String('a', field='prefix', usedefault=True,
                               desc='slicetimed output prefix')


class SiemensSliceTimingOutputSpec(TraitedSpec):
    timecorrected_files = OutputMultiPath(traits.Either(traits.List(File(exists=True)),
                                                        File(exists=True)),
                                          desc='slice time corrected files')



class SiemensSliceTiming(SPMCommand):
    """Use spm to perform slice timing correction.
    This subclass checks the input file to automatically
    set slice ordering parameters.

    Auto detection of slice order only if images
    are from Siemens and converted with dcm2nii from Nov 2013 or later.

    http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=19

    Examples
    --------

    >>> from nipype.interfaces.spm import SliceTiming
    >>> st = SliceTiming()
    >>> st.inputs.in_files = 'functional.nii'
    >>> st.inputs.num_slices = 32
    >>> st.inputs.time_repetition = 6.0
    >>> st.inputs.time_acquisition = 6. - 6./32.
    >>> st.inputs.slice_order = list(range(32,0,-1))
    >>> st.inputs.ref_slice = 1
    >>> st.run() # doctest: +SKIP

    """
    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if not isdefined(val):
            if opt == 'time_repetition':
                val = self.set_time_repetition()

            if opt == 'time_acquisition':
                val = self.set_time_acquisition()

            if opt == 'slice_order':
                val = self.set_slice_order()

            if opt == 'ref_slice':
                val = self.set_ref_slice()

            if opt == 'num_slices':
                val = self.set_num_slices()

        super(SiemensSliceTiming, self)._format_arg(opt, spec, val)


    def _check_in_files(self):
        if isinstance(self.inputs.in_files, str):
            in_files = [self.inputs.in_files]
        else:
            in_files = self.inputs.in_files

        for f in in_files:
            if not op.exists(f):
                raise IOError('Expected an existing file in `in_files`')

        return in_files

    @staticmethod
    def _get_n_slices(in_file):
        return nib.load(in_file).header.get_n_slices()

    @staticmethod
    def _get_time_repetition(in_file):
        return nib.load(in_file).header.get_dim_info()[2]

    @staticmethod
    def _get_time_acquisition(in_file):
        return TA = (TRsec/nslices)*(nslices-1);

    def set_slice_order(self):
        # if slice_order == 0 %attempt to autodetect slice order
        #     fid = fopen(fMRIname1);
        #     fseek(fid,122,'bof');
        #     slice_order = fread(fid,1,'uint8')
        #     fclose(fid);
        #     if (slice_order > kNIFTI_SLICE_UNKNOWN) && (slice_order <= kNIFTI_SLICE_ALT_DEC2)
        #         fprintf('Auto-detected slice order as %d\n',slice_order);
        #     else
        #         fprintf('%s error: unable to auto-detect slice order. Please manually specify slice order or use recent versions of dcm2nii.\n');
        #         return;
        #     end;
        # end
        # if (slice_order == kNIFTI_SLICE_ALT_INC2) || (slice_order == kNIFTI_SLICE_ALT_DEC2) %sequential
        #     isSiemens = true;
        # end;
        # if (slice_order == kNIFTI_SLICE_SEQ_INC) || (slice_order == kNIFTI_SLICE_SEQ_DEC) %sequential
        #     so = [1:1:nslices];
        # else % if sequential else Interleaved
        #     if (mod(nslices,2) == 0) && (isSiemens) %even number of slices, Siemens
        #         so =[2:2:nslices 1:2:nslices ];
        #     else
        #         so =[1:2:nslices 2:2:nslices];
        #     end
        # end
        # if (mod(slice_order,2) == 0) %isDescending
        #     so = (nslices+1)-so;
        # end; %isDescending
        # so

    def _time_acquisition(self):
        error_msg = 'The time acquisition calculated from all the `in_files` are not the same, got {}.'
        return self._check_all_equal(self._get_time_acquisition, error_msg)


    def _check_all_equal(self, func, error_msg=None):
        in_files = self._check_in_files()

        values = [func(f) for f in in_files]
        if not check_equal(values):
            if error_msg is None:
                error_msg = 'The values from {} are not the same, got {{}}.'.format(func.__name__)
            raise ValueError(error_msg.format(values))

        return values[0]

    def _n_slices(self):
        error_msg = 'The number of z slices for all `in_files` are not the same, got {}.'
        return self._check_all_equal(self._get_n_slices, error_msg)

    def set_time_repetition(self):
        if isdefined(self.inputs.time_repetition):
            return self.inputs.time_repetition

        error_msg = 'The TR calculated from all the `in_files` are not the same, got {}.'
        self.inputs.re self._check_all_equal(self._get_time_repetition, error_msg)



class SiemensSliceTiming(SPMCommand):
    """Use spm to perform slice timing correction.
    Some options are automated for Siemens acquisitions
    compared to spm.SliceTiming correction.

    Auto detection of slice order, i.e., slice_order == [0] and slice_mode == 'unknown'
    only works for images from Siemens and converted with dcm2nii from Nov 2013 or later.

    http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=19

    Nipype Inputs
    -------------
    All the inputs are the same as SliceTiming, although most are not mandatory anymore.

    slice_mode: str
        Choices:
            'unknown': auto detect if images are from Siemens and converted with dcm2nii from Nov 2013 or later #kNIFTI_SLICE_UNKNOWN
            'seq_inc': sequential ascending kNIFTI_SLICE_SEQ_INC = 1; %1,2,3,4
            'seq_dec': sequential descending kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
            'alt_inc': Siemens: interleaved ascending with odd number of slices, interleaved for other vendors kNIFTI_SLICE_ALT_INC = 3; %1,3,2,4
            'alt_dec': descending interleaved kNIFTI_SLICE_ALT_DEC = 4; %4,2,3,1
            'alt_inc2': Siemens interleaved ascending with even number of slices kNIFTI_SLICE_ALT_INC2 = 5; %2,4,1,3
            'alt_dec2': Siemens interleaved descending with even number of slices kNIFTI_SLICE_ALT_DEC2 = 6; %3,1,4,2

        Default: 'unknown'
        If left to default will read the TR from the nifti image header.

    Examples
    --------

    >>> from nipype.interfaces.spm import SliceTiming
    >>> st = SliceTiming()
    >>> st.inputs.in_files = 'functional.nii'
    >>> st.inputs.num_slices = 32
    >>> st.inputs.time_repetition = 6.0
    >>> st.inputs.time_acquisition = 6. - 6./32.
    >>> st.inputs.slice_order = list(range(32,0,-1))
    >>> st.inputs.ref_slice = 1
    >>> st.run() # doctest: +SKIP

    """

    input_spec = SiemensSliceTimingInputSpec
    output_spec = SiemensSliceTimingOutputSpec

    _jobtype = 'temporal'
    _jobname = 'st'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if opt == 'in_files':
            return scans_for_fnames(filename_to_list(val),
                                    keep4d=False,
                                    separate_sessions=True)
        return super(SiemensSliceTiming, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['timecorrected_files'] = []

        filelist = filename_to_list(self.inputs.in_files)
        for f in filelist:
            if isinstance(f, list):
                run = [fname_presuffix(in_f, prefix=self.inputs.out_prefix) for in_f in f]
            else:
                run = fname_presuffix(f, prefix=self.inputs.out_prefix)
            outputs['timecorrected_files'].append(run)
        return outputs
