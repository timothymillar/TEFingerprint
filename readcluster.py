import numpy as np


class ReadCluster(object):
    """"""
    def __init__(self, reference, name, strand, tips):
        self.reference = reference
        self.name = name  # need to escape regex
        self.strand = strand
        self.start = np.min(tips)
        self.end = np.max(tips)
        self.count = len(tips)
        self.mean = np.mean(tips)
        self.median = np.meadian(tips)
        self.mode = np.bincount(tips).argmax()

    def __id(self):
        return '_'.join((self.reference,
                         self.name,
                         self.strand,
                         self.start))

    def __attributes(self):
        template = 'ID={0};Name={1};count={2};mean={3};median={4};mode={5}'
        return template.format(self.__id(),
                               self.name,
                               self.count,
                               self.mean,
                               self.median,
                               self.mode)

    def __str__(self):
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(self.reference,
                               '.',
                               '.',
                               self.start,
                               self.end,
                               '.',
                               self.strand,
                               '.',
                               self.attributes)
