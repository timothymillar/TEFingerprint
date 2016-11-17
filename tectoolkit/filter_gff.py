#! /usr/bin/env python

import sys
import argparse
import gffutils


class FilterGffProgram(object):
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)
        self.run()

    def parse_args(self, args):
        """

        :param args:
        :return:
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('input_gff',
                            nargs=1,
                            help='A single gff file to be filtered')
        parser.add_argument('-f','--filters',
                            nargs='+',
                            help='Filters to apply')
        try:
            arguments = parser.parse_args(args)
        except:
            parser.print_help()
            sys.exit(0)
        else:
            return arguments

    def _parse_filter(self, string):
        filt = {}
        if '>=' in string:
            filt['operator'] = '>='
        elif '<=' in string:
            filt['operator'] = '<='
        elif '==' in string:
            filt['operator'] = '=='
        elif '>' in string:
            filt['operator'] = '>'
        elif '<' in string:
            filt['operator'] = '<'
        elif '=' in string:
            filt['operator'] = '='
        filt['attribute'], filt['value'] = string.split(filt['operator'])
        return filt

    def run(self):
        gff_db = GffFilterDB(self.args.input_gff)
        filters = [self._parse_filter(string) for string in self.args.filters]
        gff_db.filter_by_attributes(filters)
        for feature in gff_db.db.all_features():
            print(feature)


class GffFilterDB(object):
    """"""
    def __init__(self, input_gff):
        self.db = gffutils.create_db(input_gff,
                                     dbfn=':memory:',
                                     keep_order=True,
                                     merge_strategy='merge',
                                     sort_attribute_values=True)

    def descendants(self, feature):
        for child in self.db.children(feature):
            yield child
            self.descendants(child)

    def ancestors(self, feature):
        for child in self.db.parents(feature):
            yield child
            self.ancestors(child)

    def _apply_filter(self, feature, filt):
        if filt['operator'] == '=':
            return feature.attributes[filt['attribute']][0] == filt['value']
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

    def selected(self, feature, filters):
        if {filt['attribute'] for filt in filters}.issubset(feature.attributes.keys()):
            return all([self._apply_filter(feature, filt) for filt in filters])
        else:
            return False

    def relative_selected(self, feature, filters):
        return any([self.selected(feature, filters),
                    any([self.selected(f, filters) for f in self.descendants(feature)]),
                    any([self.selected(f, filters) for f in self.ancestors(feature)])])

    def filter_by_attributes(self, filters):
        for feature in self.db.all_features():
            if self.relative_selected(feature, filters):
                pass
            else:
                self.db.delete(feature)
