#!/usr/bin/env python
"""Classes describing coverage along the genome"""
__author__ = "Matthew Peterson"

class SequenceCoverage(object):
    """The coverage along a genome"""

    def __init__(self, infile):
        self._coverage = {} # Would this work better as a dict of dicts?

    def _read_bu_wig(self, handle, sequence_name):
        """Reads a 'wig' file .  Note that this is not the same as the wiggle
        standard.  This file is a tab-delimited list of

        position\ttotal coverage\treverse coverage\tforward coverage

        :param bu_wig_handle: handle to be read from
        :type bu_wig_handle: file
        :param sequence_name: chromosome name
        :type sequence_name: string
        """
        self.chromosomes[sequence_name] = () 
        for line in handle:
            position, total, reverse, forward = line.rstrip("\r\n").split()
            self.coverage[(sequence, position)] = (int(reverse), int(forward))

    def _read_swig(self, swigfile_handle, on_disk = False):
        """Reads in a SWIG file.  This file consists of delimited lines
        consisting of

        sequence\tposition\treverse coverage\tforward coverage
        """
        for line in swigfile_handle():
            sequence, position, reverse, forward = line.rstrip("\r\n").split()
            self._coverage[(sequence, position)] = (int(reverse), int(forward))

    def get_coverage(self, sequence, position):
        """
        Get the coverage at a given position

        :param sequence: the sequence
        :type sequence: string
        :param position: position on the sequence
        :type position: int
        """
        return self._coverage[(sequence, position)]

    def get_coverage_as_array(self, length):
        """
        Returns the coverage for a given chromosome as a numpy.array.
        """
        coverage = np.array(zeros(length))
        for i in range(0,length):
            coverage[i] = self.coverage(sequence, i)

class DiskBasedGenomeCoverage(object):
    """
    Disk-based Genome Coverage, using pytables.  Used to allow for the loading
    of larger genomes without worrying about memory limitations.
    """
    # TODO: Look into MemMaps and PyTables.
    pass

