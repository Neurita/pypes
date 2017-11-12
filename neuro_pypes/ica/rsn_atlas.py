# -*- coding: utf-8 -*-
"""
Resting state networks atlas container class
"""
from collections import OrderedDict

import nilearn.plotting as niplot
import nilearn.image as niimg


class RestingStateNetworks:
    """ A container class to parse and return useful values/images
    from RSN templates.

    Parameters
    ----------
    img_file: str
        Path to a img-like RSN templates file.

    txt_file: str
        Path to a text file with the labels for `img_file`.
        The expected structure of each line of
        this file is  <name of the RSN>, <list of indices>.
        For example:

        Basal Ganglia, 21
        Auditory, 17
        Sensorimotor, 7, 23, 24, 38, 56, 29
        Visual, 46, 64, 67, 48, 39, 59
        Default-Mode, 50, 53, 25, 68
        Attentional, 34, 60, 52, 72, 71, 55
        Frontal, 42, 20, 47, 49

    start_from_one: bool
        If True it means that the `txt_file` volume indices start from 1.
        From 0 otherwise. Be careful, the default is True!
    """
    def __init__(self, img_file, txt_file, start_from_one=True):
        self._img_file       = img_file
        self._txt_file       = txt_file
        self._start_from_one = start_from_one

        self.network_names = self._network_names()
        self._img          = niimg.load_img(self._img_file)
        self._self_check()

    def iter_networks(self):
        """Yield idx (what is in the text_file) and
        image of each RSN volume from the text file."""
        for idx, name in self.network_names.items():
            yield idx, self._get_img(idx)

    def iter_networks_names(self):
        """Yield idx (what is in the text_file) and
        image of each RSN volume from the text file."""
        for idx, name in self.network_names.items():
            yield idx, name

    def _network_names(self):
        """Return OrderedDict[int]->str, with the index and the name of each
        RSN."""
        names_idx = self._read_labels_file()
        return OrderedDict([(idx, name) for name, idxs in names_idx.items()
                                        for idx in idxs])

    def _get_img(self, network_index):
        """ Return one RSN given the index in the labels file."""
        img_idx = self._img_index(network_index)
        return niimg.index_img(self._img, img_idx)

    def _img_index(self, network_index):
        """Return the correspoding image index for the given network index."""
        if self._start_from_one:
            return network_index - 1

        return network_index

    def _self_check(self):
        """Simple content check."""
        n_labels = len(self.network_names)
        n_images = self._img.shape[-1]

        if n_labels == n_images:
            return

        if n_labels > n_images:
            raise ValueError('The number of labels is larger than the number '
                             'of images. Got {} and {}.'.format(n_labels, n_images))

        # print('The number of volumes in the image is different from the number '
        #       'of labels in the text file.\n I am going to use only the ones '
        #       ' in the text file.')

    def _read_labels_file(self):
        """ Read the text file and return a dict[str->List[int]] with network
        names and blob indices.
        """
        lines = [l.rstrip('\n') for l in open(self._txt_file).readlines()]

        netblobs = OrderedDict()
        for l in lines:
            pcs = l.split(',')
            netname = pcs[0]
            blobs   = [int(idx) for idx in pcs[1:]]

            netblobs[netname] = blobs
        return netblobs

    def plot_all(self):
        names = self.network_names
        for idx, rsn in enumerate(niimg.iter_img(self._img)):
            disp = niplot.plot_roi(rsn, title=names.get(idx, None))

    def join_networks(self, network_indices):
        """Return a NiftiImage containing a binarised version of the sum of
        the RSN images of each of the `network_indices`."""
        oimg = self._get_img(network_indices[0]).get_data()
        for idx in network_indices[1:]:
            oimg += self._get_img(idx).get_data()

        return niimg.new_img_like(self._get_img(network_indices[0]),
                                  oimg.astype(bool))

    def __iter__(self):
        return (img for idx, img in self.iter_networks())

    def __len__(self):
        return len(self.network_names)

