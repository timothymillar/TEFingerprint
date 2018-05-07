#! /usr/bin/env python3


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


if __name__ == "__main__":
    pass
