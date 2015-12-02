"""
Private helper functions
"""
import os.path as op


def flatten_list(list_of_lists):
    """ Will convert a list of lists in a list with the items inside each sub-list.

    Parameters
    ----------
    list_of_lists: list[list[object]]

    Returns
    -------
    list
    """
    if not list_of_lists:
        return []
    if isinstance(list_of_lists[0], list):
        return [l.pop() for l in list_of_lists]
    return list_of_lists


def format_pair_list(pair_list, **kwargs):
    """ Given a list of 2-tuples of str, calls format with `kwargs` for each
    item in the 2-tuples and return the formatted list of 2-tuples.

    Parameters
    ----------
    pair_list: list of 2-tuples of str

    kwargs: keyword arguments
        Arguments for the format function of each string in the 2-tuples.

    Returns
    -------
    formatted_pair_list: list of 2-tuples of str
    """
    return [(s1.format(**kwargs), s2.format(**kwargs)) for s1, s2 in pair_list]


def _check_list(str_or_list):
    """ If `str_or_list` is a list will return it as it is. If it is a str, will return a
    list with the str inside.

    Parameters
    ----------
    str_or_list: str or list

    Returns
    -------
    list
    """
    if str_or_list is None:
        return None

    if isinstance(str_or_list, list):
        return str_or_list

    if isinstance(str_or_list, str):
        return [str_or_list]

    raise ValueError('Expected a `str` or a `list`, ' \
                     'got {}.'.format(type(str_or_list)))


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
