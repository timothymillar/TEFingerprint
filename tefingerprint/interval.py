#! /usr/bin/env python3

import numpy as np


def singular_contains(x_lower,
                      x_upper,
                      y,
                      y_lower_field='start',
                      y_upper_field='stop'):
    """
    Identify intervals contained within a single interval.

    This is the inclusive contents.

    :param x_lower: inclusive lower bound of single interval
    :type x_lower: int
    :param x_upper: inclusive upper bound of single interval
    :type x_upper: int
    :param y: array of intervals
    :type y: :class:`numpy.array`
    :param y_lower_field: field name of lower bound of intervals
    :type y_lower_field: str
    :param y_upper_field: field name of upper bound of intervals
    :type y_upper_field: str

    :return: an array of boolean values
    :rtype: :class:`numpy.array`[bool]
    """
    return np.logical_and(x_lower <= y[y_lower_field],
                          y[y_upper_field] <= x_upper)


def singular_within(x_lower,
                    x_upper,
                    y,
                    y_lower_field='start',
                    y_upper_field='stop'):
    """
    Identify intervals containing a single interval.

    This is the inclusive within.

    :param x_lower: lower bound of single interval
    :type x_lower: int
    :param x_upper: upper bound of single interval
    :type x_upper: int
    :param y: array of intervals
    :type y: :class:`numpy.array`
    :param y_lower_field: field name of lower bound of intervals
    :type y_lower_field: str
    :param y_upper_field: field name of upper bound of intervals
    :type y_upper_field: str

    :return: an array of boolean values
    :rtype: :class:`numpy.array`[bool]
    """
    return np.logical_and(y[y_lower_field] <= x_lower,
                          x_upper <= y[y_upper_field])


def singular_intersects_lower(x_lower,
                              x_upper,
                              y,
                              y_lower_field='start'):
    """
    Identify intervals who's lower bound is contained within a
    single interval.

    This is an inclusive contains.

    :param x_lower: lower bound of single interval
    :type x_lower: int
    :param x_upper: upper bound of single interval
    :type x_upper: int
    :param y: array of intervals
    :type y: :class:`numpy.array`
    :param y_lower_field: field name of lower bound of intervals
    :type y_lower_field: str

    :return: an array of boolean values
    :rtype: :class:`numpy.array`[bool]
    """
    return np.logical_and(x_lower <= y[y_lower_field],
                          y[y_lower_field] <= x_upper)


def singular_intersects_upper(x_lower,
                              x_upper,
                              y,
                              y_upper_field='stop'):
    """
    Identify intervals who's upper bound is contained within a
    single interval.

    This is an inclusive contains.

    :param x_lower: lower bound of single interval
    :type x_lower: int
    :param x_upper: upper bound of single interval
    :type x_upper: int
    :param y: array of intervals
    :type y: :class:`numpy.array`
    :param y_upper_field: field name of upper bound of intervals
    :type y_upper_field: str

    :return: an array of boolean values
    :rtype: :class:`numpy.array`[bool]
    """
    return np.logical_and(x_lower <= y[y_upper_field],
                          y[y_upper_field] <= x_upper)


def singular_intersects(x_lower,
                        x_upper,
                        y,
                        y_lower_field='start',
                        y_upper_field='stop'):
    """
    Identify intervals intersecting a single interval.

    This is an inclusive intersect.

    :param x_lower: lower bound of single interval
    :type x_lower: int
    :param x_upper: upper bound of single interval
    :type x_upper: int
    :param y: array of intervals
    :type y: :class:`numpy.array`
    :param y_lower_field: field name of lower bound of intervals
    :type y_lower_field: str
    :param y_upper_field: field name of upper bound of intervals
    :type y_upper_field: str

    :return: an array of boolean values
    :rtype: :class:`numpy.array`[bool]
    """
    return np.logical_or(np.logical_and(y[y_lower_field] <= x_upper,
                                        x_upper <= y[y_upper_field]),
                         np.logical_and(y[y_lower_field] <= x_lower,
                                        x_lower <= y[y_upper_field]))


def lengths(x, lower_field='start', upper_field='stop'):
    """
    Calculate lengths of intervals in an array.

    :param x: array of intervals
    :type x: :class:`numpy.array`
    :param lower_field: field name of lower bound of intervals
    :type lower_field: str
    :param upper_field: field name of upper bound of intervals
    :type upper_field: str

    :return: an array of interval lengths
    :rtype: :class:`numpy.array`[int]
    """
    return 1 + x[upper_field] - x[lower_field]


def length_of_contains(x, y, lower_field='start', upper_field='stop'):
    """
    Calculate sum total length of intervals from array y contained
    in each interval of array x.

    Arrays x and y must have the same field names for lower and
    upper bounds.

    :param x: array of intervals
    :param y: array of intervals
    :param lower_field: field name of lower bound of intervals
    :type lower_field: str
    :param upper_field: field name of upper bound of intervals
    :type upper_field: str

    :return: an array of summed interval lengths, equal in length
        to array x
    :rtype: :class:`numpy.array`[int]
    """
    singles = zip(x[lower_field], x[upper_field])
    masks = (singular_contains(l,
                               u,
                               y,
                               y_lower_field=lower_field,
                               y_upper_field=upper_field)
             for l, u in singles)
    return np.fromiter((np.sum(lengths(y[mask],
                                       lower_field=lower_field,
                                       upper_field=upper_field))
                        for mask in masks),
                       dtype=np.int64,
                       count=len(x))


if __name__ == "__main__":
    pass
