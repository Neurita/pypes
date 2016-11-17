# neurita/pypes

This is a Python module for brain PET and multimodal MRI processing.

It contains specific pipelines to pre-process this type of images.

First you will need to convert those DICOM files to NifTI files.
I recommend you to use [dcm2niix](https://github.com/rordenlab/dcm2niix) for this.
The you will probably want to clean and rename those files to have a homogeneous database.

I've tested these on ~140 subjects with very good results.
Without modifying parameters, 95% of the images were correctly registered and processed.
Very few subjects needed manual editings such as neck removal or AC-PC reorientation.

The spatial normalization here are only done with SPM12 for now.
