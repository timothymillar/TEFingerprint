#! /usr/bin/env python3

import numpy as np
#from tefingerprint import util
from .dtype import remove_field as _dtype_remove_field


def bind(x, y):
    """
    Bind the fields of two structered numpy arrays of the same shape
    into a single array.

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


def remove_field(array, field):
    """
    Return a copy of an array without the specified field.

    :param array: a numpy array
    :type array: :class:`numpy.array`
    :param field: the name of a filed in the array
    :type field: str

    :return: a numpy array
    :rtype: :class:`numpy.array`
    """
    dtype = _dtype_remove_field(array.dtype, field)
    new = np.empty(len(array), dtype=dtype)
    for field in dtype.names:
        new[field] = array[field]
    return new


if __name__ == "__main__":
    pass
