"""
TransNet
========

    TransNet is a Python package for working with transcriptional regulatory
    networks.
"""

from __future__ import absolute_import

import sys

import chipseq
import genome

if sys.version_info[:2] < (2,6):
    m = "Python version 2.6 or greater is required for TransNet (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

__author__ = "Matthew Peterson"

