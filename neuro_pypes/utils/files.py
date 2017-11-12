# -*- coding: utf-8 -*-
"""
Helper functions to manage external files.
"""
import re
from   os          import path as op
from   glob        import glob
from   functools   import wraps


def get_vox_dims(volume):
    import nibabel as nb
    if isinstance(volume, list):
        volume = volume[0]
    nii = nb.load(volume)
    hdr = nii.header
    voxdims = hdr.get_zooms()
    return [float(voxdims[0]), float(voxdims[1]), float(voxdims[2])]


def get_data_dims(volume):
    import nibabel as nb
    if isinstance(volume, list):
        volume = volume[0]
    nii = nb.load(volume)
    hdr = nii.header
    datadims = hdr.get_data_shape()
    return [int(datadims[0]), int(datadims[1]), int(datadims[2])]


def get_affine(volume):
    import nibabel as nb
    nii = nb.load(volume)
    return nii.affine


def niftiimg_out(f):
    """ Picks a function whose first argument is an `img` or a sequence of imgs, processes its
    data and returns a numpy array. This decorator wraps this numpy array
    into a nibabel.Nifti1Image."""
    import nibabel as nib
    import nilearn.image as niimg

    @wraps(f)
    def wrapped(*args, **kwargs):
        r = f(*args, **kwargs)

        img = niimg.load_img(args[0])
        return nib.Nifti1Image(r, affine=img.get_affine(), header=img.header)

    return wrapped


def get_extension(filepath, check_if_exists=False, allowed_exts=None):
    """Return the extension of filepath.

    Parameters
    ----------
    filepath: string
        File name or path

    check_if_exists: bool

    allowed_exts: dict
        Dictionary of strings, where the key if the last part of a complex ('.' separated) extension
        and the value is the previous part.
        For example: for the '.nii.gz' extension I would have a dict as {'.gz': '.nii'}
        Default: {'.gz': '.nii'}

    Returns
    -------
    ext: str
        The extension of the file name or path
    """
    if allowed_exts is None:
        allowed_exts = {'.gz': '.nii'}

    try:
        rest, ext = op.splitext(filepath)
        if ext in allowed_exts:
            alloweds = allowed_exts[ext]
            _, ext2 = op.splitext(rest)
            if ext2 in alloweds:
                ext = ext2 + ext
    except:
        raise
    else:
        return ext


def add_extension_if_needed(filepath, ext):
    """Add the extension ext to fpath if it doesn't have it.

    Parameters
    ----------
    filepath: str
        File name or path

    ext: str
        File extension

    Returns
    -------
    filepath: str
        File name or path with extension added, if needed.
    """
    if not filepath.endswith(ext):
        filepath += ext

    return filepath


def remove_ext(filepath):
    """Removes the extension of the file.

    Parameters
    ----------
    filepath: str
        File path or name

    Returns
    -------
    filepath: str
        File path or name without extension
    """
    return filepath[:filepath.rindex(get_extension(filepath))]


def extension_duplicates(regexp_pair_list):
    """ Return a new list with the pair items of `regexp_pair_list`
    have all the '.nii' replaced by '.nii.gz'.
    This is useful for the Datasink regexp_substitutions when
    you don't know/care what extension the output image will have,
    then you put all of them like:
    "(r"/rc1[\w]+_corrected\.nii$",  "/coreg_gm.nii")"

    and then call this function to add the duplicates with modified
    extension of these same pairs.

    Parameters
    ----------
    regexp_pair_list: list of 2-tuple of str

    Returns
    -------
    mod_regexp_pair_list: list of 2-tuple of str
    """
    replace_ext = lambda x: x.replace('.nii', '.nii.gz')
    dups = [(replace_ext(pair[0]), replace_ext(pair[1]))
            for pair in regexp_pair_list
            if '.nii$' in pair[0]]
    return dups


def rename(in_files, suffix=None):
    """Rename all the files in `in_files` adding the `suffix` keeping its
    extension and basedir."""
    import os.path as path
    from nipype.utils.filemanip import (filename_to_list, split_filename,
                                        list_to_filename)
    out_files = []
    for idx, filename in enumerate(filename_to_list(in_files)):
        base, name, ext = split_filename(filename)
        if suffix is None:
            new_name = name + ('_%03d' % idx) + ext
        else:
            new_name = name + suffix + ext

        out_files.append(path.join(base, new_name))

    return list_to_filename(out_files)


def find_files_in(dirpath, file_pattern, pat_type='fnmatch'):
    """ Find files in `dirpath` without recursivity.

    Parameters
    ----------
    dirpath: str
        Folder where to search for file names.

    file_pattern: str
        File pattern to be matched.

    pat_type: str
        The type of pattern in `file_pattern`.
        Choices: 'fnmatch', 're.search', 're.match'.

    Returns
    -------
    files: List[str]
        List of paths to the files that match file_pattern.
        `dirpath` is included in each of its items.
    """
    if pat_type == 'fnmatch':
        files = glob(op.join(dirpath, file_pattern))
    elif pat_type == 're.search':
        regex = re.compile(file_pattern)
        files = [f for f in glob(op.join(dirpath, '*')) if regex.search(f)]
    elif pat_type == 're.match':
        regex = re.compile(file_pattern)
        files = [f for f in glob(op.join(dirpath, '*')) if regex.match(f)]
    else:
        raise ValueError("Expected one of ('fnmatch', 're.search' or 're.match') for "
                         "`pat_type` parameter, got {}.".format(pat_type))

    return files


def fetch_one_file(dirpath, file_pattern, file_extension=None, extra_prefix=None,
                   extra_suffix=None, pat_type='fnmatch'):
    """ Return the unique file path in dirpath that matches fnmatch file_pattern.
    Add the extra_prefix to try the search again, if the file_pattern finds more than one file.

    Parameters
    ----------
    dirpath:

    file_pattern:
        File pattern to be matched.

    file_extension:
        Extension of the file.

    extra_prefix:

    extra_suffix:

    pat_type: str
        The type of pattern in `file_pattern`.
        Choices: 'fnmatch', 're.search', 're.match'.

    Returns
    -------
    file_path: str
        The path to the unique file that matches the conditions inside `dirpath`.
        `dirpath` is included.

    Raises
    ------
    IOError
        If `dirpath` doesn't exist or the list of returned files is empty.

    ValueError
        IF the choice for `pat_type` is not valid.

    """
    if file_extension is not None:
        file_fnmatch = file_pattern + file_extension
    else:
        file_fnmatch = file_pattern

    files = find_files_in(dirpath, file_pattern, pat_type=pat_type)

    if not files:
        raise IOError('Expected at least one file that matched the '
                      'pattern {} in {}.'.format(file_fnmatch, dirpath))

    if len(files) > 1:
        if extra_prefix is None:
            extra_prefix = ''

        if extra_suffix is None:
            extra_suffix = ''

        if not extra_prefix and not extra_suffix:
            raise IOError('Found more than one file that matched the '
                          'pattern {} in {}: {}'.format(file_fnmatch, dirpath, files))

        else:
            # TODO: be careful, this might only work with fnmatch
            return fetch_one_file(dirpath, extra_prefix + file_pattern + extra_suffix,
                                  file_extension=file_extension,
                                  pat_type=pat_type)

    return files[0]


def save_object(obj, filename):
    """ Save `obj` in `filename` as a pickle."""
    import pickle

    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
