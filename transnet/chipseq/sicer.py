"""
Parser for the SICER ChIP-seq program.  This is designed to parse
the summary file generated as output
"""
from transnet.chipseq.chip_peak import ChipPeak

__author__ = 'Matthew Peterson'

def parse(handle, method = "sicer"):
    """
    Returns an iterator of SicerPeaks

    Parameters:
    - `handle`: A file handle to the SICER output file to be parsed
    """

    peak_class_map = {"sicer" : SicerPeak, "sicer_rb" : SicerRBPeak}

    if method not in peak_class_map:
        raise KeyError("Unsupported method. Select one of 'sicer' or "
                       "sicer_rb")

    handle = iter(handle)

    # Skip blank lines
    for line in handle:
        peak = peak_class_map[method](line)
        yield peak

def read(handle):
    """
    Returns a list, if we need to keep this around

    Parameters:
    - `handle`: A file handle to the SICER output file to be read
    """
    return list(parse(handle))

class SicerPeak(ChipPeak):
    """
    A region of binding identified by SICER
    """
    def __init__(self, chromosome, start, stop, island_read_count,
                 control_read_count, p_value, fold_change, fdr):
        """
        Creates a new SicerPeak object

        Parameters:
        - `chromosome`: The chromosome/sequence the peak is found on
        - `start`: The start position (zero-based) on the chromosome
        - `stop`: The stop position (zero-based) on the chromosome
        - `island_read_count`: Number of reads in the peak/island
        - `control_read_count`: Number of reads on control
        - `p_value`: Probability that region is enriched
        - `fold_change`: Fold enrichment over control
        - `fdr`: False discovery rate
        """
        super(SicerPeak, self).__init__(chromosome, start, stop)
        self.island_read_count = island_read_count
        self.control_read_count = control_read_count
        self.p_value = p_value
        self.fold_change = fold_change
        self.fdr = fdr

class SicerRBPeak(ChipPeak):
    """
    A region of binding identified by SICER (no background)
    """
    def __init__(self, input_line):
        """
        Creates a new SicerRBPeak

        Reads a line from the output file from SICER_rb.sh - This is of the
        format

        chromosome start stop score

        :param input_line: line from output file
        :type input_line: string
        """
        tokens = input_line.rstrip("\r\n").split()
        super(SicerRBPeak, self).__init__(tokens[0], int(tokens[1]),
                                        int(tokens[2]))
        self.score = float(tokens[3]) 

