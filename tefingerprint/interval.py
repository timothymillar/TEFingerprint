#! /usr/bin/env python

import numpy as np


def singular_contains(x_lower, x_upper, y, y_lower_field='start', y_upper_field='stop'):
    return np.logical_and(x_lower <= y[y_lower_field], y[y_upper_field] <= x_upper)


def singular_within(x_lower, x_upper, y, y_lower_field='start', y_upper_field='stop'):
    return np.logical_and(y[y_lower_field] <= x_lower, x_upper <= y[y_upper_field])


def singular_intersects(x_lower, x_upper, y, y_lower_field='start', y_upper_field='stop'):
    return np.logical_or(np.logical_and(y[y_lower_field] <= x_upper, x_upper <= y[y_upper_field]),
                         np.logical_and(y[y_lower_field] <= x_lower, x_lower <= y[y_upper_field]))


def lengths(x, lower='start', upper='stop'):
    return 1 + x[upper] - x[lower]


def length_of_contains(x, y, lower='start', upper='stop'):
    singles = zip(x[lower], x[upper])
    masks = (singular_contains(l, u, y, y_lower_field=lower, y_upper_field=upper) for l, u in singles)
    return np.fromiter((np.sum(lengths(y[mask], lower=lower, upper=upper)) for mask in masks), dtype=np.int64, count=len(x))

