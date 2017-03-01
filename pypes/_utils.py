# -*- coding: utf-8 -*-
"""
Private helper functions
"""


def grep(lines, substr):
    """ Return a list of strings from `lines` that have
    `substr` as a substring.
    """
    return [l for l in lines if substr in l]


def check_equal(lst):
    """ Return True if all items in `lst` are equal, False otherwise.
    Note that check_equal([1, True]) is True.
    """
    return not lst or lst.count(lst[0]) == len(lst)


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
        lst = []
        for l in list_of_lists:
            lst.extend(l)
        return lst
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


def get_values_map_keys(records, keyidx=0):
    """ Given a dict of str->2-tuples, e.g.:
            {'anat': [('modality', 'anat'), ('image_file', 'anat_hc.nii.gz')],
             'pet':  [('modality', 'pet'),  ('image_file', 'pet_fdg.nii.gz')],

        or

       Given a list of list of 2-tuples of str, e.g.:
            [[('modality', 'anat'), ('image_file', 'anat_hc.nii.gz')],
              ('modality', 'pet'),  ('image_file', 'pet_fdg.nii.gz')],


    Will return the unique values of each record value, in this case:
    {'modality', 'image_file'}.

    Parameters
    ----------
    values_maps_dict: Dict[str->2-tuple]

    Returns
    -------
    keys: set[str]
    """
    if not records or records is None:
        return []

    if isinstance(records, dict):
        itemset = records.values()
    elif isinstance(records, list):
        itemset = records
    else:
        raise NotImplementedError('Expected a `dict` or a `list of list` as `records, '
                                  'got {}.'.format(type(records)))

    crumb_args = set()
    for items in itemset:
        crumb_args = crumb_args.union(set([t[keyidx] for t in items]))

    return crumb_args