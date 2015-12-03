"""
Private helper functions
"""


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
    return [(s[0].format(**kwargs), s[1].format(**kwargs)) for s in pair_list]


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


