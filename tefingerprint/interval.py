#! /usr/bin/env python

import numpy as np


def singular_contains(x_lower, x_upper, y, y_lower_field='start', y_upper_field='stop'):
    return np.logical_and(x_lower <= y[y_lower_field], y[y_upper_field] <= x_upper)


def singular_within(x_lower, x_upper, y, y_lower_field='start', y_upper_field='stop'):
    return np.logical_and(y[y_lower_field] <= x_lower, x_upper <= y[y_upper_field])


def singular_intersects_lower_bound(x_lower, x_upper, y, y_lower_field='start'):
    return np.logical_and(x_lower <= y[y_lower_field], y[y_lower_field] <= x_upper)


def singular_intersects_upper_bound(x_lower, x_upper, y, y_upper_field='stop'):
    return np.logical_and(x_lower <= y[y_upper_field], y[y_upper_field] <= x_upper)


def singular_intersects(x_lower, x_upper, y, y_lower_field='start', y_upper_field='stop'):
    return np.logical_or(np.logical_and(y[y_lower_field] <= x_upper, x_upper <= y[y_upper_field]),
                         np.logical_and(y[y_lower_field] <= x_lower, x_lower <= y[y_upper_field]))


def lengths(x, lower='start', upper='stop'):
    return 1 + x[upper] - x[lower]

