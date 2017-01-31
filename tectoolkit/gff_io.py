#! /usr/bin/env python


class GffFeature(object):
    """"""
    def __init__(self,
                 seqid='.',
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
        self.tags = kwargs
        self.children = []

    def attribute_names(self):
        return set(self.tags.keys())

    def _parse_attributes(self, attributes):
        if attributes is not None:
            return ';'.join(tuple('{0}={1}'.format(key, self.tags[key]) for key in sorted(self.tags.keys())))
        else:
            return'.'

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
                               self._parse_attributes(self.tags))

if __name__ == '__main__':
    pass