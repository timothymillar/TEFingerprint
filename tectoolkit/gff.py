#! /usr/bin/env python


class GffFeature(object):
    """"""
    def __init__(self,
                 seqid,
                 source='.',
                 ftype='.',
                 start='.',
                 end='.',
                 score='.',
                 strand='.',
                 phase='.',
                 **kwargs):
        self.seqid = seqid
        self.source = source
        self.type = ftype
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = kwargs
        self.children = []

    def attribute_names(self):
        return set(self.attributes.keys())

    def _parse_attributes(self, attributes):
        if attributes is not None:
            return ';'.join(tuple('{0}={1}'.format(key, value) for key, value in self.attributes.items()))
        else:
            return'.'

    def add_children(self, *args):
        for child in args:
            assert isinstance(child, GffFeature)
            assert "ID" in self.attributes
            child.attributes["Parent"] = self.attributes["ID"]
            self.children.append(child)

    def __str__(self):
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(self.seqid,
                               self.source,
                               self.type,
                               self.start,
                               self.end,
                               self.score,
                               self.strand,
                               self.phase,
                               self._parse_attributes(self.attributes))

    def _str_nested(self):
        return [self.__str__()] + [child._str_nested() for child in self.children]

    def _flatten(self, item):
        """

        :param lst:
        :return:
        """
        if isinstance(item, list):
            for element in item:
                for item in self._flatten(element):
                    yield item
        else:
            yield item

    def __format__(self, code):
        assert code in {'nested', 'single'}
        if code == "nested":
            return '\n'.join(self._flatten(self._str_nested()))
        else:
            return self.__str__()


if __name__ == '__main__':
    pass