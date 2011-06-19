#!/usr/bin/env python

from transnet.chipseq.chip_peak import ChipPeak

def parse(peaks_handle, chromosome = "Genome"):
    """
    Reads in a set of peaks
    """
    for line in peaks_handle:
        tokens = line.rstrip("\r\n").split("\t")
        peak = LogNormalPeak(chromosome, int(tokens[0]), int(tokens[1]),
                             int(tokens[2]), int(tokens[3]))
        yield peak

def read(peaks_handle, chromosome = "Genome"):
    """
    Returns a list of peaks.
    """
    peaks = parse(peaks_handle, chromosome)
    return list(peaks)

class LogNormalPeak(ChipPeak):
    """
    A peak called from a log-normal distribution. (Describe procedure here)
    """
    def __init__(self, chromosome, start, stop, height, shift):
        self.chromosome = chromosome
        self.chrom_start = start
        self.chrom_end = stop
        self.height = height
        self.shift = shift

    def score(self):
        return self.height
