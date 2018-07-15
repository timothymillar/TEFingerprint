#! /usr/bin/env python3

import numpy as np


def as_dict(element):
    """
    Convert a numpy element into a dictionary.

    :param element: an element from a numpy array with named fields
    :type element: :class:`numpy.void`

    :return: a dictionary
    :rtype: dict
    """
    return {field: element[field] for field in element.dtype.names}


def flatten(item):
    """
    Flatten a nested numpy element.

    :param item: a numpy element
    :type item: np.void[any] | any

    :return: A generator
    :rtype: generator[any]
    """
    if isinstance(item, np.void) or isinstance(item, tuple):
        for element in item:
            for item in flatten(element):
                yield item
    else:
        yield item


if __name__ == "__main__":
    pass