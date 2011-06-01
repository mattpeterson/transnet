"""
Parser for the SICER ChIP-seq program.  This is designed to parse
the summary file generated as output
"""
from transnet.interval import Interval

__author__ = 'petersmw'

def parse(handle):
    """
    Returns an iterator of SicerPeaks

    Parameters:
    - `handle`: A file handle to the SICER output file to be parsed
    """
    handle = iter(handle)

    # Skip blank lines
    for line in handle:
        tokens = line.rstrip("\r\n").split()
        peak = SicerPeak(tokens[0], int(tokens[1]), int(tokens[2]),
                         int(tokens[3]), float(tokens[4]), float(tokens[5]),
                         float(tokens[6]), float(tokens[7]))
        yield peak

def read(handle):
    """
    Returns a list, if we need to keep this around

    Parameters:
    - `handle`: A file handle to the SICER output file to be read
    """
    return list(parse(handle))

class SicerPeak(Interval):
    """
    A region of binding identified by SICER
    """
    def __init__(self, chromosome, start, stop, island_read_count,
                 control_read_count, p_value, fold_change, fdr):
        """
        Constructor
        """
        super(SicerPeak, self).__init__(chromosome, start, stop)
        self.island_read_count = island_read_count
        self.control_read_count = control_read_count
        self.p_value = p_value
        self.fold_change = fold_change
        self.fdr = fdr

