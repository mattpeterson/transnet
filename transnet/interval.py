"""
Class descrbing an interval on the genome.
"""
__author__ = 'Matthew Peterson'

class Interval(object):
    """
    An interval along the genome.  Superclass for any feature mapped to
    the Genome (ChIP peaks, genes, intergenic regions, etc.)
    """
    def __init__(self, chromosome = "", start=1, stop = 1):
        """
        Create a new Interval.
        """
        if start < 0:
            raise ValueError("Start index cannot be less than zero.")

        self.chromosome = chromosome
        self.chrom_start = start
        self.chrom_end = stop

    def overlaps(self, other):
        """
        Tests whether or not this interval overlaps another.
        """
        if self.chromosome != other.chromosome:
            return False

        if self.chrom_start <= other.chrom_end and self.chrom_end >= other.chrom_start:
            return True

        return False

    def get_genome_features(self, genome, genic=True, intergenic=True):
        """
        Get overlapping features and classifications
        """
        hits = []

        for f in genome.features(True):
            if self.overlaps(f):
                hits.append(f)

        return hits

    def _filter_hits(self):
        hits = set()
        return hits

    def __str__(self):
        """
        String representation of interval
        """
        return "%s:%d-%d" % (self.chromosome, self.chrom_start, self.chrom_end)

