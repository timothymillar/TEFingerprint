#! /usr/bin/env python

import re
import sys
import argparse
import itertools


class FilterGffProgram(object):
    """Main class for the filter gff program"""
    def __init__(self, arguments):
        """
        Init method for :class:`FilterGffProgram`.

        :param arguments: A list of commandline arguments to be parsed for the filter_gff program
        """
        self.args = self.parse_args(arguments)
        with open(self.args.gff[0], 'r') as f:
            feature_strings = f.read().splitlines()
        result = filter_features(feature_strings,
                                 column_filters=self.args.column_filters,
                                 attribute_filters=self.args.attribute_filters)
        result = '\n'.join(result)
        print(result)

    @staticmethod
    def parse_args(args):
        """
        Defines an argument parser to handle commandline inputs for the filter_gff program.

        :param args: A list of commandline arguments for the filter_gff program

        :return: A dictionary like object of arguments and values for the filter_gff program
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('gff',
                            nargs=1,
                            help='A single gff file to be filtered')
        parser.add_argument('-c', '--column-filters',
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
        parser.add_argument('-a', '--attribute-filters',
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
    """
    Parses a string containing a single gff3 feature and returns a dictionary of the first eight columns.

    :param string: a one line string containing a single gff3 feature
    :type string: str

    :return: a dictionary containing data from the first eight gff columns
    :rtype: dict[str, str]
    """
    return dict(zip(['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'], string.split('\t')[:-1]))


def _parse_gff_attributes(string):
    """
    Parses a string containing a single gff3 feature and returns a dictionary of attributes from the attributes column.

    :param string: a one line string containing a single gff3 feature
    :type string: str

    :return: a dictionary containing attribute data from the attributes column
    :rtype: dict[str, str]
    """
    return {k: v for k, v in [item.split('=') for item in string.split('\t')[-1].split(';')]}


def _equivalence_filter(x, y):
    """
    Tests a pair of strings for equivalence by attempting to convert them to floats and
    falling back to string comparison.

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


def _non_equivalence_filter(x, y):
    """
    Tests a pair of strings for non-equivalence by attempting to convert them to floats and
    falling back to string comparison.

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


# Function dispatch based on operator passed by user
_FILTER_DISPATCH = {'==': _equivalence_filter,
                    '=': _equivalence_filter,
                    '!=': _non_equivalence_filter,
                    '>=': lambda x, y: float(x) >= float(y),
                    '>': lambda x, y: float(x) > float(y),
                    '<=': lambda x, y: float(x) <= float(y),
                    '<': lambda x, y: float(x) < float(y)}


def _parse_filter_string(string):
    """
    Parse a filter string to identify the attribute, operator and value.

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
    """
    Check if feature meats requirement of filter.

    :param feature: A dictionary of attribute-value pairs
    :type feature: dict[str, str]
    :param filt: A dictionary with keys 'attribute', 'operator' and 'value'
    :param filt: dict[str, str]

    :return:
    """
    try:
        feature[filt['attribute']]
    except KeyError:
        # attribute is not present so feature does not match filter
        return False
    else:
        # split attribute values in case it is a comma separated list
        values = feature[filt['attribute']].split(',')
        # check if any values meet filter requirement
        return any([_FILTER_DISPATCH[filt['operator']](value, filt['value']) for value in values])


def filter_features(feature_strings, column_filters=None, attribute_filters=None):
    """
    Filter a list of gff3 formatted feature strings. Filters may be applied to feature attributes or
    standardised gff3 columns (excluding the attributes column).

    A filter is a dictionary with the keys 'attribute', 'operator' and 'value' where attribute is the
    target 'attribute' to filter on, 'value' is the value to compare the attribute against and
    'operator' determines the comparison to be performed.

    If multiple filters are specified, a feature must meet the requirement of each of them to be selected.

    In the case that an attribute contains a comma-separated list of values, the feature is selected if
    any of the values meet the requirement of the filter.

    :param feature_strings: a list of gff3 formatted feature strings
    :type feature_strings: list[str]
    :param column_filters: filters to apply to first 8 gff columns
    :type column_filters: list[dict[str, str]]
    :param attribute_filters: filters to apply to attributes
    :type attribute_filters: list[dict[str, str]]

    :return: a subset of gff3 formatted feature strings
    :rtype: generator[str]
    """
    if column_filters is None:
        column_filters = []
    if attribute_filters is None:
        attribute_filters = []
    keep = (True for _ in feature_strings)
    if len(column_filters) > 0:
        columns = [_parse_gff_columns(string) for string in feature_strings]
        column_filters = [_parse_filter_string(string) for string in column_filters]
        for filt in column_filters:
            match = (_matches_filter(col, filt) for col in columns)
            keep = [k and m for k, m in zip(keep, match)]
    if len(attribute_filters) > 0:
        attributes = [_parse_gff_attributes(string) for string in feature_strings]
        attribute_filters = [_parse_filter_string(string) for string in attribute_filters]
        for filt in attribute_filters:
            match = (_matches_filter(attr, filt) for attr in attributes)
            keep = [k and m for k, m in zip(keep, match)]
    return itertools.compress(feature_strings, keep)


if __name__ == '__main__':
    FilterGffProgram(sys.argv)
