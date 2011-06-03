"""
Classes and methods for scoring a ChIP-seq experiment using a Poisson
background model.
"""
from transnet.interval import Interval

__author__ = "Matthew Peterson"

def parse(handle, chromosome = "Genome"):
    """
    Parse a set of peaks.  Returns an iterator

    Parameters:
    - `handle`: The handle to be read
    - `chromosome`: The chromosome
    """
    for line in handle:
        tokens = line.rstrip("\r\n").split("\t")
        peak = PoissonPeak(chromosome, int(tokens[1]), int(tokens[2]),
                           int(tokens[3]), float(tokens[4]), int(tokens[5]))
        yield peak

def read(handle, chromosome = "Genome"):
    """
    Read in a set of peaks from an output file.

    Parameters:
    - `handle`: The handle of the file to be read
    - `chromosome`: The
    """
    return list(parse(handle, chromosome))

class PoissonPeak(Interval):
    """
    A peak generated by the MATLAB or Python scripts for scoring against
    a Poisson background
    """
    def __init__(self, chromosome, start, stop, height, mean_p, shift):
        """
        Create a new Poisson-identified Peak
        
        Parameters:
        - `chromosome`: Chromosome/sequence peak is on
        - `start`: Position (1-based) on the chromosome of start of peak
        - `stop`: Position (1-based) on the chromosome of peak end
        - `mean_p`: Mean p-value of the region
        - `shift`: Max of cross_correlation(forward, reverse)
        """
        super(PoissonPeak, self).__init__(chromosome, start - 1, stop - 1)
        self.height = height
        self.mean_pval = mean_p
        self.shift = shift
