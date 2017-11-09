#! /usr/bin/env python

import numpy as np


def append_dtypes(*args):
    descr = []
    for dtype in args:
        descr += dtype.descr
    return np.dtype(descr)


def remove_array_field(array, field):
    dtype = remove_dtype_field(array.dtype, field)
    new = np.empty(len(array), dtype=dtype)
    for field in dtype.names:
        new[field] = array[field]
    return new


def remove_dtype_field(dtype, field):
    return np.dtype([i for i in dtype.descr if i[0] != field])


def extract_dtype_field(dtype, field):
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


def flatten_dtype(dtype, prefix=''):
    def flatten_descr(item, prefix=''):
        if isinstance(item, tuple):
            if isinstance(item[1], list):
                prefix += str(item[0]) + '_'
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


def flatten_dtype_fields(dtype, prefix=''):
    for name in flatten_dtype(dtype, prefix=prefix).names:
        yield name


def bind_arrays(x, y):
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
    return {field: element[field] for field in element.dtype.names}


def quote_str(value):
    if isinstance(value, str):
        return '"{0}"'.format(value)
    else:
        return str(value)