#! /usr/bin/env python

import re
import argparse
import numpy as np


class FilterGffProgram(object):
    """Main class for the filter_gff program"""
    def __init__(self, arguments):
        """
        Init method for :class:`FilterGffProgram`.

        :param arguments: A list of commandline arguments to be parsed for the filter_gff program
        """
        self.args = self.parse_args(arguments)
        feature_strings = read_gff(self.args.input_gff[0])
        result = filter_features(feature_strings,
                                 column_filters=self.args.column_filters,
                                 attribute_filters=self.args.attribute_filters)
        print(result)

    def parse_args(self, args):
        """
        Defines an argument parser to handle commandline inputs for the filter_gff program.

        :param args: A list of commandline arguments for the filter_gff program

        :return: A dictionary like object of arguments and values for the filter_gff program
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('input_gff',
                            nargs=1,
                            help='A single gff file to be filtered')
        parser.add_argument('-c', '--column_filters',
                            nargs='*',
                            help=("List of filters to apply to standard GFF3 columns.\n"
                                  "A valid filter takes the form '<attribute><operator><value>'"
                                  "where <attribute> is the name of a GFF attribute, "
                                  "<operator> is one of '=', '==', '!=', '>=', '<=', '>' or '<' "
                                  "and the value of the GFF attribute is compared to <value> using the operator\n"
                                  "The list of filters is applied additively (i.e. a feature must meet all filters) "
                                  "and, if a feature is selected, all of it's ancestors and descendants "
                                  "will also be included in the output.\n"
                                  "Operators '=', '==' and '!=' will attempt to compare values as floating point "
                                  "numbers if possible and otherwise compare values as strings. "
                                  "Operators '>=', '<=', '>' and '<' will coerce values "
                                  "to floating point numbers before comparison."))
        parser.add_argument('-a', '--attribute_filters',
                            nargs='*',
                            help=("List of filters to apply to attributes.\n"
                                  "A valid filter takes the form '<attribute><operator><value>'"
                                  "where <attribute> is the name of a GFF attribute, "
                                  "<operator> is one of '=', '==', '!=', '>=', '<=', '>' or '<' "
                                  "and the value of the GFF attribute is compared to <value> using the operator\n"
                                  "The list of filters is applied additively (i.e. a feature must meet all filters) "
                                  "and, if a feature is selected, all of it's ancestors and descendants "
                                  "will also be included in the output.\n"
                                  "Operators '=', '==' and '!=' will attempt to compare values as floating point "
                                  "numbers if possible and otherwise compare values as strings. "
                                  "Operators '>=', '<=', '>' and '<' will coerce values "
                                  "to floating point numbers before comparison."))
        return parser.parse_args(args)


def _parse_gff_columns(string):
    return dict(zip(['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'], string.split('\t')[:-1]))


def _parse_gff_attributes(string):
    return {k: v for k, v in [item.split('=') for item in string.split('\t')[-1].split(';')]}


def _equivalence_filter(x, y):
    try:
        return float(x) == float(y)
    except ValueError:
        return x == y


def _non_equivalence_filter(x, y):
    try:
        return float(x) != float(y)
    except ValueError:
        return x != y


_FILTER_DISPATCH = {'==': _equivalence_filter,
                    '=': _equivalence_filter,
                    '!=': _non_equivalence_filter,
                    '>=': lambda x, y: float(x) >= float(y),
                    '>': lambda x, y: float(x) > float(y),
                    '<=': lambda x, y: float(x) <= float(y),
                    '<': lambda x, y: float(x) < float(y)}


def _parse_filter_string(string):
    """
    Parse a filter string to identify the attribute, operator and value

    :param string: A valid filter string in the form '<attribute><operator><value>'
    :type string: str

    :return: A dictionary with the keys 'attribute', 'operator' and 'value'
    :rtype: dict[str, str]
    """
    operator = re.findall(">=|<=|==|!=|>|<|=", string)
    if len(operator) > 1:
        raise ValueError('There is more than one operator in filter: "{0}"'.format(string))
    elif len(operator) < 1:
        raise ValueError('Could not find a valid operator in filter: "{0}"'.format(string))
    else:
        filt = {}
        operator = operator[0]
        filt['operator'] = operator
        filt['attribute'], filt['value'] = string.split(filt['operator'])
        return filt


def _matches_filter(feature, filt):
    return _FILTER_DISPATCH[filt['operator']](feature[filt['attribute']], filt['value'])


def read_gff(file):
    with open(file, 'r') as f:
        feature_strings = f.read().splitlines()
    return feature_strings


def filter_features(feature_strings, column_filters=None, attribute_filters=None):
    if column_filters is None:
        column_filters = []
    if attribute_filters is None:
        attribute_filters = []
    keep = np.empty(len(feature_strings), dtype=bool)
    keep.fill(True)
    if len(column_filters) > 0:
        columns = [_parse_gff_columns(string) for string in feature_strings]
        column_filters = [_parse_filter_string(string) for string in column_filters]
        for filt in column_filters:
            match = np.array([_matches_filter(col, filt) for col in columns], dtype=bool)
            keep = np.logical_and(keep, match)
    if len(attribute_filters) > 0:
        attributes = [_parse_gff_attributes(string) for string in feature_strings]
        attribute_filters = [_parse_filter_string(string) for string in attribute_filters]
        for filt in attribute_filters:
            match = np.array([_matches_filter(attr, filt) for attr in attributes], dtype=bool)
            keep = np.logical_and(keep, match)
    return '\n'.join(np.array(feature_strings)[keep])


