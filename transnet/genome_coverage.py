"""Classes describing coverage along the genome"""
__author__ = "Matthew Peterson"


class GenomeCoverage(object):
    """The coverage along a genome"""

    def __init__(self, infile):
        self.chromosomes = {}


    def _read_bu_wig(self, sequence_name):
        """Reads a 'wig' file .  Note that this is not the same as the wiggle
        standard.  This file is a tab-delimited list of

        position\ttotal coverage\treverse coverage\tforward coverage

        :param bu_wig_handle: handle to be read from
        :type bu_wig_handle: file
        :param sequence_name: chromosome name
        :type sequence_name: string
        """
        pass

    def _read_swig(self, swigfile_handle, on_disk = False):
        """Reads in a SWIG file.  This file consists of delimited lines
        consisting of

        sequence\tposition\treverse coverage\tforward coverage
        """
        for line in swigfile_handle():

class SequenceCoverage(object):
        """The coverage along a sequence/chromosome"""
        pass

class DiskSequenceCoverage(SequenceCoverage):
    """Disk-based sequence coverage, for larger chromosomes"""
    pass

