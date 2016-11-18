# neurita/pypes

This is a Python module for brain PET and multimodal MRI processing.

It contains specific pipelines to pre-process this type of images.

First you will need to convert those DICOM files to NifTI files.
I recommend you to use [dcm2niix](https://github.com/rordenlab/dcm2niix) for this.
The you will probably want to clean and rename those files to have a homogeneous database.

I've tested these on ~140 subjects with very good results.
Without modifying parameters, 95% of the images were correctly registered and processed.
Very few subjects needed manual editings such as neck removal or AC-PC reorientation.


## License

Licensed under the Apache License, Version 2.0.


## Author

The main author of this project is [Alexandre Savio](http://alexsavio.github.io/).
This work was done while he was employed by MD Dr. Igor Yakushev
in the Nuklear Medizin Department (NUK) of the Klinikum rechts der Isar of
the Technical University of Munich, Germany.

### Many thanks to other contributors:

- [Michael Schutte](https://github.com/schutte) for improvements and the DTI/Camino
pipeline while doing his Master Thesis in the NUK.
