"""
Function utilities that use pandas.
"""

import pandas as pd


def add_table_headers(values, axis0_labels, axis1_labels, **kwargs):
    """ Return a pandas.DataFrame with the content of `values`.
    This DataFrame will have each row labeled by `axis0_labels`, and
    each column by `axis1_labels`.

    Parameters
    ----------
    values: np.ndarray
        MxN matrix
        M is the length of axis0_labels.
        N is the lenght of axis1_labels.

    axis0_labels: array or any type
        This will be used as index in the DataFrame construction.

    axis1_labels: array or any type
        This will be used as column in the DataFrame construction.

    kwargs: keyword arguments
        Additional columns to be added to the resulting DataFrame.
        The argument names are the column name and the values must be sequences of length M.

    Returns
    -------
    labeled_measure: pandas.DataFrame

    Notes
    -----
    My previous solution was more complex and I had this function
    prepared before realizing this was ridiculously simple.
    """
    df = pd.DataFrame(values, index=axis0_labels, columns=axis1_labels)

    for k, v in kwargs.items():
        if len(v) != len(axis0_labels):
            raise AttributeError('The value for argument {} should have length {} but has '
                                 'length {}.'.format(k, len(axis0_labels), len(v)))
        df[k] = v

    return df


def write_tabbed_excel(filepath, dict_df, **kwargs):
    """ Write each df in `dict_df` in a tab of the Excel spreadsheet in `filepath`.

    Parameters
    ----------
    filepath: str

    dict_df: Dict[pandas.DataFrame]

    kwargs:
        columns: list of str
            List of the columns of `df` to be written in the sheet.

        other arguments to be passed to the df.to_excel function.
    """
    from pandas import ExcelWriter

    with ExcelWriter(filepath) as writer:
        for key, df in dict_df.items():
            # write df to the tab
            df.to_excel(writer, sheet_name=key, **kwargs)