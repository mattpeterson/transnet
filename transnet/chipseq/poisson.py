"""
Classes and methods for scoring a ChIP-seq experiment using a Poisson
background model.
"""
from transnet.chipseq.chip_peak import ChipPeak

__author__ = "Matthew Peterson"

def parse(handle, chromosome = "Genome"):
    """
    Parse a set of peaks.  Returns an iterator

    :param handle: The handle to be read
    :type handle: file
    :param chromosome: the chromosome the peaks are found on
    :type chromosome: string
    """
    for line in handle:
        tokens = line.rstrip("\r\n").split("\t")
        peak = PoissonPeak(chromosome, int(tokens[0]), int(tokens[1]),
                           int(tokens[2]), float(tokens[3]), int(tokens[4]))
        yield peak

def read(handle, chromosome = "Genome"):
    """
    Read in a set of peaks from an output file.

    Parameters:
    - `handle`: The handle of the file to be read
    - `chromosome`: The
    """
    return list(parse(handle, chromosome))

def score(ip_coverage, bg_coverage=None):
    """Scores an experiment against a background lane.  Will scale the
    background coverage to match that of the IP lane

    :param ip_coverage: coverage in IP lane
    :type ip_coverage: GenomeCoverage
    """
    raise NotImplementedError("Scoring has yet to be implemented.")

def _score_no_bg(ip_coverage, num_std, min_length):
    raise NotImplementedError("Scoring has yet to be implemented.")


def _score_vs_bg(ip_coverage, bg, min_length):
    pass

class PoissonPeak(ChipPeak):
    """
    A peak generated by the MATLAB or Python scripts for scoring against
    a Poisson background
    """
    def __init__(self, chromosome, start, stop, height, mean_p, shift):
        """
        Create a new Poisson-identified Peak

        :param chromosome: chromosome
        :type chromosome: string
        :param start: start position on chromosome
        :type start: int
        :param stop: stop position on chromosome
        :type stop: int
        :param mean_p: mean p-value in region
        :type mean_p: float
        :param shift: shift between forward-reverse peaks
        :type shift: int
        """
        super(PoissonPeak, self).__init__(chromosome, start - 1, stop - 1)
        self.height = height
        self.mean_pval = mean_p
        self.shift = shift

    def score(self):
        return self.height


