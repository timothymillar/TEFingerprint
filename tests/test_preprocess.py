#! /usr/bin/env python
import os
import pytest
import subprocess
from tefingerprint import preprocess


@pytest.mark.parametrize('query,answer',
                         [('AGATC', 'GATCT'),
                          ('aGtcn', 'NGACT')])
def test_reverse_complement(query, answer):
    assert preprocess.reverse_complement(query) == answer

