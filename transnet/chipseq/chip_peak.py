"""Base class describing a ChIP-seq Peak"""
from transnet.interval import Interval
from transnet.genome import IntergenicRegion

class ChipPeak(Interval):
    """An enriched region, or "peak" identified by ChIP-seq analysis."""
    def __init__(self, chromosome, start, stop):
        """
        Create a new ChIP-seq peak

        :param chromosome: the chromosome the peak is on
        :type chromosome: string
        :param start: the start position on the chromosome
        :type start: int
        :param stop: the stop position on the chromosome
        :type stop: int
        """
        super(ChipPeak, self).__init__(chromosome, start, stop)

    def score():
        """Returns the 'score' for a peak.  Not implemented in this abstract
        class, there is no concept of what the 'score' is"""
        raise NotImplementedError("score() is not implemented for the base" +
                                  "Class")

    def get_regulated_genes(self, annotation, filter_hits = True):
        """Gets a set of regulated genes, as well as their classification.
        Perhaps 'regulated' is not the correct word - 'implicated' may be
        better.

        Classifications include:
        * US - Upstream
        * DS - Downstreams
        * G - Genic

        :param annotation: the annotation used to call genes
        :type annotation: Genome
        :return a list of tuples of the form (Gene, classification)
        :rtype tuple
        """
        hits = set()
        upstream = set()
        downstream = set()
        genic = set()

        for f in annotation.features(True):
            if self.overlaps(f):
                hits.add(f)

        if filter_hits:
            self._filter_hits(hits)

        for h in hits:
            if h.__class__.__name__ == "Gene":
                genic.add(h)
            elif h.__class__.__name__ == "IntergenicRegion":
                # Now check to see which genes are up and downstream
                if h.left_gene.strand == "-":
                    upstream.add(h.left_gene)
                else:
                    downstream.add(h.left_gene)

                if h.right_gene.strand == "+":
                    upstream.add(h.right_gene)
                else:
                    downstream.add(h.right_gene)
        
        # Add upstream and downstream genes for genic peaks
        genes = list(annotation.features(False))
        for g in genic:
            gene_idx = genes.index(g)
            
            previous_idx = gene_idx - 1
            next_idx = gene_idx + 1
            
            # Bad news for linear chromosomes. Fix this.
            if next_idx == len(gene_idx):
                next_idx = 0
            if previous_idx == -1:
                previous_idx = len(gene_idx) - 1
            
            if genes[previous_idx].strand == "-":
                upstream.add(g)
            else:
                downstream.add(g)
            
            if genes[next_idx].strand == "+":
                upstream.add(g)
            else:
                downstream.add(g)

        return upstream, downstream, genic

    def _filter_hits(self, hits):
        """Remove genic regions also implicated by intergenic overlaps.

        :param hits: set of hits returned by get_regulated_regions
        :type hits: set
        """
        to_remove = set()

        for h in hits:
            if h.__class__ == IntergenicRegion:
                if h.left_gene in hits:
                    to_remove.add(h.left_gene)
                if h.right_gene in hits:
                    to_remove.add(h.right_gene)

        for h in to_remove:
            hits.remove(h)
    
    def to_bed(self, name="peak"):
        """
        Writes the interval as a BED field
        """
        return "%s\t%d\t%d\t%s\t%f" % (self.chromosome, self.chrom_start, 
                                       self.chrom_end, self.name, self.score())
