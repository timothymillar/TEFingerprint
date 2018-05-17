#! /usr/bin/env python3

import numpy as np


def append(*args):
    """
    Append multiple dtypes into a single dtype.

    :param args: numpy dtypes
    :type args: :class:`numpy.dtype`

    :return: a single numpy dtype
    :rtype: :class:`numpy.dtype`
    """
    descr = []
    for dtype in args:
        descr += dtype.descr
    return np.dtype(descr)


def remove_field(dtype, field):
    """
    Return a copy of an array without the specified field.

    :param dtype: a numpy dtype
    :type dtype: :class:`numpy.dtype`
    :param field: the name of a field in the array
    :type field: str

    :return: a numpy dtype
    :rtype: :class:`numpy.dtype`
    """
    return np.dtype([i for i in dtype.descr if i[0] != field])


def extract_field(dtype, field):
    """
    Returns the dtype for a single field from a complex numpy dtype.

    :param dtype: a numpy dtype
    :type dtype: :class:`numpy.dtype`
    :param field: the name of a field in the dtype
    :type field: str

    :return: a numpy dtype
    :type dtype: :class:`numpy.dtype`
    """
    return np.dtype([i for i in dtype.descr if i[0] == field])


def flatten_fields(dtype, prefix='', sep='_'):
    """
    Flatten a nested numpy dtype by appending nested field names.

    :param dtype: a numpy dtype
    :type dtype: :class:`numpy.dtype`
    :param prefix: optional prefix
    :type prefix: str
    :param sep: separator between appended field names
    :type sep: str

    :return: a flattened numpy dtype
    :rtype: :class:`numpy.dtype`
    """
    def flatten_descr(item, prefix=''):
        if isinstance(item, tuple):
            if isinstance(item[1], list):
                prefix += str(item[0]) + sep
                for item in flatten_descr(item[1], prefix=prefix):
                    yield item
            else:
                yield (prefix + item[0], item[1])
        elif isinstance(item, list):
            for element in item:
                for item in flatten_descr(element, prefix=prefix):
                    yield item
        else:
            pass

    return np.dtype(list(flatten_descr(dtype.descr, prefix=prefix)))


def flatten_field_names(dtype, prefix='', sep='_'):
    """
    Flatten a nested numpy dtype by appending nested field names and
    return field names only .

    :param dtype: a numpy dtype
    :type dtype: :class:`numpy.dtype`
    :param prefix: optional prefix
    :type prefix: str
    :param sep: separator between appended field names
    :type sep: str

    :return: field names of the flattened dtype
    :rtype: generator[str]
    """
    for name in flatten_fields(dtype, prefix=prefix, sep=sep).names:
        yield name


if __name__ == "__main__":
    pass
