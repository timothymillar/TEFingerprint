#! /usr/bin/env python


_COLUMN_ENCODING = [('%', '%25'),
                    ('\t', '%09'),
                    ('\r', '%0D')]

_ATTRIBUTE_ENCODING = _COLUMN_ENCODING + [(';', '%3B'),
                                          ('=', '%3D'),
                                          ('&', '%26'),
                                          (',', '%2C')]


def encode_column(string):
    """
    Escape special symbols for first 8 gff3 columns

    :param string: value of column
    :type string: str

    :return: escaped string
    :rtype str:
    """
    for symbol, escape in _COLUMN_ENCODING:
        string = string.replace(symbol, escape)
    return string


def decode_column(string):
    """
    Un-escape special symbols for first 8 gff3 columns

    :param string: value of column
    :type string: str

    :return: un-escaped string
    :rtype str:
    """
    for symbol, escape in _COLUMN_ENCODING:
        string = string.replace(escape, symbol)
    return string


def encode_attribute(string):
    """
    Escape special symbols for gff3 attribute key or value

    :param string: key or value of attribute
    :type string: str

    :return: escaped string
    :rtype str:
    """
    for symbol, escape in _ATTRIBUTE_ENCODING:
        string = string.replace(symbol, escape)
    return string


def decode_attribute(string):
    """
    Un-escape special symbols for gff3 attribute key or value

    :param string: key or value of attribute
    :type string: str

    :return: un-escaped string
    :rtype str:
    """
    for symbol, escape in _ATTRIBUTE_ENCODING:
        string = string.replace(escape, symbol)
    return string


if __name__ == '__main__':
    pass
