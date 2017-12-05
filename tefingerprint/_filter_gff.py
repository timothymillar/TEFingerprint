#! /usr/bin/env python


import re
import fnmatch
from tefingerprint.gff import *


COLUMN_NAMES = ["seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes"]


def parse_feature(string):
    columns = dict(zip(COLUMN_NAMES,
                       string.strip('\n').split()))
    attributes = columns["attributes"].split(';')
    attributes = {k: v for k, v in (attribute.split('=')
                                    for attribute in attributes)}
    del columns["attributes"]

    columns = {k: decode_column(v)
               for k, v in columns.items()}
    attributes = {decode_attribute(k): decode_attribute(v)
                  for k, v in attributes.items()}

    dictionary = attributes
    dictionary.update(columns)
    return dictionary


def format_feature(dictionary):
    column_values = (encode_column(str(dictionary[field]))
                     for field in COLUMN_NAMES[0:8])
    columns = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t'
    columns = columns.format(*column_values)

    attribute_names = [field for field in dictionary.keys()
                       if field not in COLUMN_NAMES]
    attribute_names.sort()
    attributes_pairs = ((encode_attribute(str(field)),
                         encode_attribute(str(dictionary[field])))
                        for field in attribute_names)

    attributes = ':'.join('{}={}'.format(field, value)
                          for field, value in attributes_pairs)

    return columns + attributes


def _eq(x, y):
    """
    Tests a pair of strings for equivalence by attempting to convert
    them to floats and falling back to string comparison.

    :param x: string
    :type x: str
    :param y: string
    :type y: str

    :return: boolean
    :rtype: bool
    """
    try:
        return float(x) == float(y)
    except ValueError:
        return x == y


def _neq(x, y):
    """
    Tests a pair of strings for non-equivalence by attempting to convert
    them to floats and falling back to string comparison.

    :param x: string
    :type x: str
    :param y: string
    :type y: str

    :return: boolean
    :rtype: bool
    """
    try:
        return float(x) != float(y)
    except ValueError:
        return x != y


# Function dispatch based on operator passed as string
_OPERATOR_DISPATCH = {'==': _eq,
                      '=': _eq,
                      '!=': _neq,
                      '>=': lambda x, y: float(x) >= float(y),
                      '>': lambda x, y: float(x) > float(y),
                      '<=': lambda x, y: float(x) <= float(y),
                      '<': lambda x, y: float(x) < float(y)}


def parse_filter_string(string):
    operator = re.findall(">=|<=|==|!=|>|<|=", string)
    if len(operator) > 1:
        message = 'There is more than one operator in filter: "{0}"'
        raise ValueError(message.format(string))
    elif len(operator) < 1:
        message = 'Could not find a valid operator in filter: "{0}"'
        raise ValueError(message.format(string))
    else:
        operator = operator[0]
        field, value = string.split(operator)
        field = field.strip()
        value = value.strip()
        return {'field': field, 'operator': operator, 'value': value}


def apply_filter(feature, filter_, combinator):
    fields = fnmatch.filter(list(feature.keys()), filter_['field'])
    results = [_OPERATOR_DISPATCH[filter_['operator']](feature[field],
                                                       filter_['value'])
               for field in fields]
    if results == []:
        return False  # if there are no fields to match
    elif combinator == 'ANY':
        return any(results)
    elif combinator == 'ALL':
        return all(results)
    else:
        assert False


def apply_filters(feature, filters, combinator):
    results = [apply_filter(feature, filt, combinator) for filt in filters]
    if combinator == 'ANY':
        return any(results)
    elif combinator == 'ALL':
        return all(results)
    else:
        assert False
