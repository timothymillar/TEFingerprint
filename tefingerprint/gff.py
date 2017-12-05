#! /usr/bin/env python


_COLUMN_ENCODING = [('%', '%25'),
                    ('\t', '%09'),
                    ('\r', '%0D')]

_ATTRIBUTE_ENCODING = _COLUMN_ENCODING + [(';', '%3B'),
                                          ('=', '%3D'),
                                          ('&', '%26'),
                                          (',', '%2C')]


def encode_column(string):
    for symbol, escape in _COLUMN_ENCODING:
        string = string.replace(symbol, escape)
    return string


def decode_column(string):
    for symbol, escape in _COLUMN_ENCODING:
        string = string.replace(escape, symbol)
    return string


def encode_attribute(string):
    for symbol, escape in _ATTRIBUTE_ENCODING:
        string = string.replace(symbol, escape)
    return string


def decode_attribute(string):
    for symbol, escape in _ATTRIBUTE_ENCODING:
        string = string.replace(escape, symbol)
    return string


if __name__ == '__main__':
    pass
