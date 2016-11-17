#! /usr/bin/env python

import gffutils


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

    def selected(self, feature, **kwargs):
        if set(kwargs.keys()).issubset(feature.attributes.keys()):
            return all([feature.attributes[k][0] == kwargs[k] for k in kwargs.keys()])
        else:
            return False

    def relative_selected(self, feature, **kwargs):
        return any([self.selected(feature, **kwargs),
                    any([self.selected(f, **kwargs) for f in self.descendants(feature)]),
                    any([self.selected(f, **kwargs) for f in self.ancestors(feature)])])

    def filter_by_attributes(self, **kwargs):
        for feature in self.db.all_features():
            if self.relative_selected(feature, **kwargs):
                pass
            else:
                self.db.delete(feature)
