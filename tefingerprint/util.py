#! /usr/bin/env python

import numpy as np


def append_dtypes(*args):
    """
    Append multiple dtypes into a single dtype.

    :param args: numpy dtypes
    :type args: list[:class:`numpy.dtype`]

    :return: a single numpy dtype
    :rtype: :class:`numpy.dtype`
    """
    descr = []
    for dtype in args:
        descr += dtype.descr
    return np.dtype(descr)


def remove_array_field(array, field):
    """
    Return a copy of an array without the specified field.

    :param array: a numpy array
    :type array: :class:`numpy.array`
    :param field: the name of a filed in the array
    :type field: str

    :return: a numpy array
    :rtype: :class:`numpy.array`
    """
    dtype = remove_dtype_field(array.dtype, field)
    new = np.empty(len(array), dtype=dtype)
    for field in dtype.names:
        new[field] = array[field]
    return new


def remove_dtype_field(dtype, field):
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


def extract_dtype_field(dtype, field):
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


def flatten_tuple(item):
    """
    Flatten a nested numpy element.

    :param item: a numpy element
    :type item: np.void[any] | any

    :return: A generator
    :rtype: generator[any]
    """
    if isinstance(item, tuple):
        for element in item:
            for item in flatten_tuple(element):
                yield item
    else:
        yield item


def flatten_numpy_element(item):
    """
    Flatten a nested numpy element.

    :param item: a numpy element
    :type item: np.void[any] | any

    :return: A generator
    :rtype: generator[any]
    """
    if isinstance(item, np.void) or isinstance(item, tuple):
        for element in item:
            for item in flatten_numpy_element(element):
                yield item
    else:
        yield item


def flatten_dtype(dtype, prefix='', sep='_'):
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


def flatten_dtype_fields(dtype, prefix='', sep='_'):
    """
    Flatten a nested numpy dtype by appending nested field names and return field names only .

    :param dtype: a numpy dtype
    :type dtype: :class:`numpy.dtype`
    :param prefix: optional prefix
    :type prefix: str
    :param sep: separator between appended field names
    :type sep: str

    :return: field names of the flattened dtype
    :rtype: generator[str]
    """
    for name in flatten_dtype(dtype, prefix=prefix, sep=sep).names:
        yield name


def bind_arrays(x, y):
    """
    Bind the fields of two structered numpy arrays of the same shape into a single array.

    :param x: a numpy array
    :type x: :class:`numpy.array`
    :param y: a numpy array
    :type y: :class:`numpy.array`

    :return: a numpy array with the fields and values of both x and y
    :rtype: :class:`numpy.array`
    """
    for field in x.dtype.names:
        assert field not in y.dtype.names
    assert len(x) == len(y)
    new = np.empty(len(x), dtype=np.dtype(x.dtype.descr + y.dtype.descr))
    for field in x.dtype.names:
        new[field] = x[field]
    for field in y.dtype.names:
        new[field] = y[field]
    return new


def dict_of_numpy(element):
    """
    Convert a numpy element into a dictionary.

    :param element: an element from a numpy array with named fields
    :type element: :class:`numpy.void`

    :return: a dictionary
    :rtype: dict
    """
    return {field: element[field] for field in element.dtype.names}


def quote_str(value):
    """
    Returns a string wrapped in quotes or any other value as a string

    :param value: any value
    :type value: str | any

    :return: a string that may be quoted
    :rtype: str
    """
    if isinstance(value, str):
        return '"{0}"'.format(value)
    else:
        return str(value)
