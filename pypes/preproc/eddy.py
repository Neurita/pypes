"""
# Comments on the `eddy` tool from FSL FDT.

A description of the tool:
http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Eddy/UsersGuide

Our problem to run this tool instead of the good-old `eddy_correct` is the `--acqp` argument, an
acquisitions parameters file.

A detailed description of the --acpq input file is here:
http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_do_I_know_what_to_put_into_my_--acqp_file

In the following subsections I describe each of the fields needed to check and build the acquisitions parameters file.

## Phase Encoding Direction
Dicom header field: (0018,1312) InPlanePhaseEncodingDirection
The phase encoding direction is the OPPOSITE of frequency encoding direction:
- 'COL' = A/P (freqDir L/R),
- 'ROW' = L/R (freqDir A/P)


Nifti header field: "descrip.phaseDir:'+'" is for 'COL' in the DICOM InPlanePhaseEncodingDirection value.
So if you have only one phase encoding oriendation and a '+' in the previous header field,
the first 3 columns for the `--acqp` parameter file should be:
0 1 0

indicating that the scan was acquired with phase-encoding in the anterior-posterior direction.

For more info:
http://web.stanford.edu/group/vista/cgi-bin/wiki/index.php/DTI_Preprocessing_User_Manual#Frequency_vs_Phase_Encode_Direction
https://xwiki.nbirn.org:8443/xwiki/bin/view/Function-BIRN/PhaseEncodeDirectionIssues


## Effective Echo Spacing (aka dwell time)

Effective Echo Spacing (s) = 1/(BandwidthPerPixelPhaseEncode * MatrixSizePhase)

effective echo spacing = 1 / [(0019,1028) * (0051,100b component #1)] (from the archives)
https://www.jiscmail.ac.uk/cgi-bin/webadmin?A3=ind1303&L=FSL&E=quoted-printable&P=29358&B=--B_3444933351_15386849&T=text%2Fhtml;%20charset=ISO-8859-1&pending=

The dwell time is in the nifti header `descrip.dwell`, in seconds (or look at the field `time_units`).
http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.fsl.epi.html

More info:
http://lcni.uoregon.edu/kb-articles/kb-0003

## EPI factor

The EPI factor is not included in the nifti header.
You can read it using the Grassroots DICOM tool called `gdcmdump`, for example:
>>> gdcmdump -C IM-0126-0001.dcm | grep 'EPIFactor'
sFastImaging.lEPIFactor                  = 128

More info:
http://dicomlookup.com/default.htm


# The fourth element of the acquisitions parameter file

The fourth element in each row is the time (in seconds) between reading the center of the first echo and reading the
center of the last echo.
It is the "dwell time" multiplied by "number of PE steps - 1" and it is also the reciprocal of the PE bandwidth/pixel.

Total readout time (FSL) = (number of echoes - 1) * echo spacing

Total Readout Time (SPM) = 1/(BandwidthPerPixelPhaseEncode)
Since the Bandwidth Per Pixel Phase Encode is in Hz, this will give the readout time in seconds


# The `---index` argument

The index argument is a text file with a row of numbers. Each number
indicates what line (starting from 1) in the `acqp` file corresponds to
each volume in the DTI acquisition.

# What to do now with the `dcm2nii` files?

I see two options to calculate the `acqp` lines with these files.

1. We already have the `dwell` but we don't have the EPI factor.
We know that the standard in Siemens is 128 and we could stick to that.

2. Use the `slice_duration * 0.001` which is very near the calculated value.

# Summary

So, for example, if we had these acquisition parameters:

```
Phase enc. dir. P >> A
Echo spacing 0.75 [ms]
EPI factor 128
```

We should put in the `acqp` file this line:
0 1 0 0.095

"""