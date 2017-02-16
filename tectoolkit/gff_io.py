#! /usr/bin/env python


class GffFeature(object):
    """
    Simple class to output GFF3 formatted features.
    All fields required by GFF3 will default to '.'.
    Addition key word arguments are treated as named attributes.
    Attributes are automatically sorted alphabetically when converted to a string.

    :param seqid:
    :type seqid: str
    :param source:
    :type source: str
    :param ftype:
    :type ftype: str
    :param start: 1-based inclusive lower extent of the features location
    :type start: int | str
    :param end: 1-based inclusive upper extent of the features location
    :type end: int | str
    :param score:
    :type score: str
    :param strand: String indicating strandedness of feature ('+', '-' or '.')
    :type strand: str
    :param phase:
    :type phase: str
    """
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
        """
        Returns the set of attribute names used in the feature.

        :return: Set of attribute names used in the feature
        :rtype: set[str]
        """
        return set(self.tags.keys())

    def _parse_attributes(self):
        """
        Sorts and parses attributes into GFF3 format.

        :return: a string of attribute-value pairs
        :rtype: str
        """
        if self.tags is not None:
            return ';'.join(tuple('{0}={1}'.format(key, self.tags[key]) for key in sorted(self.tags.keys())))
        else:
            return'.'

    def __str__(self):
        """
        String method for :class:`GffFeature`.
        Returns a GFF3 formatted string representation of the feature.

        :return: GFF3 formatted string
        :rtype: str
        """
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(self.seqid,
                               self.source,
                               self.type,
                               self.start,
                               self.end,
                               self.score,
                               self.strand,
                               self.phase,
                               self._parse_attributes())

if __name__ == '__main__':
    pass