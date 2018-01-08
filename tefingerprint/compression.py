#! /usr/bin/env python

import bz2
import gzip
import os.path


def zopen(fname, *args, **kwargs):
    if os.path.isfile(fname):
        f = open(fname, *args, **kwargs)
        token = f.read(3)
        f.seek(0)
        if token == b'\x1f\x8b\x08':
            return gzip.GzipFile(fileobj=f)
        elif token == b'\x42\x5a\x68':
            return bz2.BZ2File(f)
        else:
            return f
    else:
        if fname.endswith('gz'):
            return gzip.open(fname, *args, **kwargs)
        elif fname.endswith('bz'):
            return bz2.open(fname, *args, **kwargs)
        else:
            return open(fname, *args, **kwargs)

if __name__ == "__main__":
    pass
