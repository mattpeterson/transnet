"""
Class descrbing an interval on the genome.
"""
__author__ = 'petersmw'

class Interval(object):

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
        Test whether this interval overlaps another on the genome.
        """
        if self.chromosome != other.chromosome:
            return False

        if self.chrom_start <= other.chrom_end and self.chrom_end >= other.chrom_start:
            return True

        return False

    def get_genome_features(self, genome):
        """
        Get overlapping features and classifications
        """
        hits = []
        
#        for g in genome.features(True):
#            if self.overlaps(g):
#                hits.append(g)

        for i in genome.intergenic_regions:
            print i.chromosome, i.chrom_start, i.chrom_end
            if self.overlaps(i):
                print i.chromosome, i.chrom_start, i.chrom_end
                hits.append(i)

        if hits:
            print(",".join([h.locus for h in hits]))

    def _filter_hits(self):
        hits = set()
        return hits

    def __str__(self):
        """
        String representation of interval
        """
        return "%s:%d-%d" % (self.chromosome, self.chrom_start, self.chrom_end)
    

  