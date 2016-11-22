#! /usr/bin/env python

import sys
import argparse
import gffutils


class FilterGffProgram(object):
    """Main class for the filter_gff program"""
    def __init__(self, arguments):
        """
        Init method for :class:`FilterGffProgram`.

        :param arguments: A list of commandline arguments to be parsed for the filter_gff program
        """
        self.args = self.parse_args(arguments)

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
        parser.add_argument('-f','--filters',
                            nargs='+',
                            help=("List of filters to apply.\n"
                                  "A valid filter takes the form '<attribute><operator><value>'"
                                  "where <attribute> is the name of a GFF attribute, "
                                  "<operator> is one of '=', '!=', '==', '>=', '<=', '>' or '<' "
                                  "and the value of the GFF attribute is compared to <value> using the operator\n"
                                  "The list of filters is applied additively (i.e. a feature must meet all filters) "
                                  "and, if a feature is selected, all of it's ancestors and descendants "
                                  "will also be included in the output.\n"
                                  "Operators '=' and '!=' will compare values as strings where "
                                  "Operators '==', '>=', '<=', '>' and '<' will try to coerce values "
                                  "to floating point numbers before comparison."))
        try:
            arguments = parser.parse_args(args)
        except:
            parser.print_help()
            sys.exit(0)
        else:
            return arguments

    def _parse_filter(self, string):
        """
        Parse a filter string to identify the attribute, operator and value

        :param string: A valid filter string in the form '<attribute><operator><value>'
        :type string: str

        :return: A dictionary with the keys 'attribute', 'operator' and 'value'
        :rtype: dict[str, str]
        """
        filt = {}
        if '>=' in string:
            filt['operator'] = '>='
        elif '<=' in string:
            filt['operator'] = '<='
        elif '==' in string:
            filt['operator'] = '=='
        elif '!=' in string:
            filt['operator'] = '!='
        elif '>' in string:
            filt['operator'] = '>'
        elif '<' in string:
            filt['operator'] = '<'
        elif '=' in string:
            filt['operator'] = '='
        filt['attribute'], filt['value'] = string.split(filt['operator'])
        return filt

    def run(self):
        """
        Run the filter_gff program with parameters specified in an instance of :class:`FilterGffProgram`.
        Imports the target gff file, subsets it by specified filters, and prints subset to stdout.
        """
        db = GffFilterDB(gffutils.create_db(self.args.input_gff[0],
                                            dbfn=':memory:',
                                            keep_order=True,
                                            merge_strategy='merge',
                                            sort_attribute_values=True))
        filters = [self._parse_filter(string) for string in self.args.filters]
        db.filter_by_attributes(filters)
        print(db)


class GffFilterDB(object):
    """Subset a gff file using a list of filters."""
    def __init__(self, db):
        """
        Init method for :class:`GffFilterDB`.

        :param input_gff: GFF file to be read into database for sub-setting.
        :type input_gff: str
        """
        self.db = db

    def __str__(self):
        return '\n'.join([str(feature) for feature in self.db.all_features()])

    def descendants(self, feature):
        """
        Recursively find all descendants of a feature.

        :param feature: Feature to find descendants of
        :type feature: :class:`gffutils.feature.Feature`

        :return: A generator of descendant features
        :rtype: generator[:class:`gffutils.feature.Feature`]
        """
        for child in self.db.children(feature):
            yield child
            self.descendants(child)

    def ancestors(self, feature):
        """
        Recursively find all ancestors of a feature.

        :param feature: Feature to find ancestors of
        :type feature: :class:`gffutils.feature.Feature`

        :return: A generator of ancestors features
        :rtype: generator[:class:`gffutils.feature.Feature`]
        """
        for child in self.db.parents(feature):
            yield child
            self.ancestors(child)

    def matches_filter(self, feature, filt):
        """
        Ascertains whether a feature meets the requirements of a single filter.
        A filter is a dictionary with the keys 'attribute', 'operator' and 'value' and values are strings.
        The 'value' of 'attribute' will be compared against the value of the gff features attribute of the same name.
        The 'operator' determines the manor of comparison.
        Operators '=' and '!=' will compare values as strings.
        Operators '==', '>=', '<=', '>' and '<' will coerce values to floats before comparison.

        :param feature: the feature to be tested
        :type feature: :class:`gffutils.feature.Feature`
        :param filt: A dictionary with the keys 'attribute', 'operator' and 'value'
        :type filt: dict[str, str]

        :return: Boolean value indicating whether the feature meets the filter criteria
        :rtype: bool
        """
        if filt['operator'] == '=':
            return feature.attributes[filt['attribute']][0] == filt['value']
        elif filt['operator'] == '!=':
            return feature.attributes[filt['attribute']][0] != filt['value']
        elif filt['operator'] == '==':
            return float(feature.attributes[filt['attribute']][0]) == float(filt['value'])
        elif filt['operator'] == '>=':
            return float(feature.attributes[filt['attribute']][0]) >= float(filt['value'])
        elif filt['operator'] == '<=':
            return float(feature.attributes[filt['attribute']][0]) <= float(filt['value'])
        elif filt['operator'] == '>':
            return float(feature.attributes[filt['attribute']][0]) > float(filt['value'])
        elif filt['operator'] == '<':
            return float(feature.attributes[filt['attribute']][0]) < float(filt['value'])

    def matches_filters(self, feature, filters):
        """
        Ascertains whether a feature meets the requirements of each of a list of filters.

        :param feature: the feature to be tested
        :type feature: :class:`gffutils.feature.Feature`
        :param filters: A list of dictionaries with the keys 'attribute', 'operator' and 'value'
        :type filters: list[dict[str, str]]

        :return: Boolean value indicating whether the feature meets all of the filter criteria
        :rtype: bool
        """
        if {filt['attribute'] for filt in filters}.issubset(feature.attributes.keys()):
            return all([self.matches_filter(feature, filt) for filt in filters])
        else:
            return False

    def relative_matches_filters(self, feature, filters):
        """
        Ascertains whether a feature or any of its relatives meets the requirements of each of a list of filters.

        :param feature: the feature to have itself and all ancestors and descendants tested
        :type feature: :class:`gffutils.feature.Feature`
        :param filters: A list of dictionaries with the keys 'attribute', 'operator' and 'value'
        :type filters: list[dict[str, str]]

        :return: Boolean value indicating whether the feature or any relative meets all of the filter criteria
        :rtype: bool
        """
        return any([self.matches_filters(feature, filters),
                    any([self.matches_filters(f, filters) for f in self.descendants(feature)]),
                    any([self.matches_filters(f, filters) for f in self.ancestors(feature)])])

    def filter_by_attributes(self, filters):
        """
        For every feature in the database of :class:`GffFilterDB`, checks if the feature or any of its relatives meets
        the requirements of each of a list of filters. If not, the feature is dropped from the database.

        :param filters: A list of dictionaries with the keys 'attribute', 'operator' and 'value'
        :type filters: list[dict[str, str]]
        """
        for feature in self.db.all_features():
            if self.relative_matches_filters(feature, filters):
                pass
            else:
                self.db.delete(feature)

if __name__ == '__main__':
    pass
