"""
Classes describing Genomes and genes
"""
from operator import attrgetter
from transnet.interval import Interval
from collections import defaultdict

__author__ = 'Matthew Peterson'

def read(handle, format, **kwargs):
    """
    Reads a genome from a BED file
    """
    if format == "broad":
        if "mapping" in kwargs:
            mapping = kwargs["mapping"]
        else:
            mapping = None
        (genes, intergenic_regions) = _read_broad_summary(handle, mapping)
    elif format == "bed":
        (genes, intergenic_regions) = _read_bed(handle)
    else:
        raise ValueError("Invalid file type.")

    genome = Genome()
    genome.intergenic_regions = intergenic_regions

    for g in genes:
        genome.gene_dict[g.locus] = g

    return genome

def _read_bed(handle):
    """
    Reads a genome annotation from a BED-formatted file
    """
    genes = []
    
    handle = iter(handle)
    for line in handle:
        tokens = line.rstrip("\r\n").split()
        chromosome = tokens[0]
        start = int(tokens[1])
        stop = int(tokens[2])
        locus = tokens[3]
        genes.append(Gene(chromosome, start, stop, locus))

    genes = _sort_genes(genes)
    intergenic_regions = _create_intergenic_regions(genes)
    

    return genes, intergenic_regions

def _read_broad_summary(handle, mapping=None):
    """
    Reads a genome from the Broad Institute's genome summaries
    """
    genes = []

    handle = iter(handle)
    # Skip first line
    next(handle)
    for line in handle:
        tokens = line.rstrip("\r\n").split("\t")
        locus = tokens[0]
        start = int(tokens[4])
        stop = int(tokens[5])
        strand = tokens[6]
        name = tokens[7]
        chromosome = tokens[8]
        if mapping:
            if chromosome not in mapping:
                raise ValueError("Chromosome '%s' not found in mapping." %
                                 chromosome)
            
            chromosome = mapping[chromosome]

        g = Gene(chromosome, start, stop, locus, strand)
        g.name = name

        genes.append(g)

    intergenic_regions = _create_intergenic_regions(genes)
    
    return genes, intergenic_regions

def _sort_genes(genes):
    return sorted(genes, key=attrgetter('chromosome', 'chrom_start'))

def _create_intergenic_regions(genes):
    _intergenic_regions = []
    for i in range(1, len(genes)):
        if genes[i].chromosome != genes[i-1].chromosome:
            continue
        _intergenic_regions.append(IntergenicRegion(genes[i-1], genes[i]))
    return _intergenic_regions

class Gene(Interval):
    """
    A gene.
    """
    annotations = defaultdict(list)

    def __init__(self, chromosome, start, stop, locus, strand = "+"):
        super(Gene, self).__init__(chromosome, start, stop)
        self.locus = locus
        self.strand = strand

    def __str__(self):
        return self.locus


class IntergenicRegion(Interval):
    """
    An intergenic region.
    """
    def __init__(self, first_gene, second_gene):
        """

        """
        self.left_gene = first_gene
        self.right_gene = second_gene
        super(IntergenicRegion, self).__init__(property(self.get_chrom),
                                               property(self.get_start),
                                               property(self.get_end))
        if first_gene.chromosome != second_gene.chromosome:
            raise ValueError("Genes must be on same chromosome.")

    def get_chrom(self):
        """
        Method for chromosome property
        """
        return self.left_gene.chromosome

    def get_start(self):
        """
        Method for chrom_start property
        """
        return self.left_gene.chrom_end

    def get_end(self):
        """
        Method for chrom_end property
        """
        return self.right_gene.chrom_end

    def get_locus(self):
        """
        Method for locus property
        """
        return self.left_gene.locus + "-" + self.right_gene.locus

    locus = property(get_locus)

    def __str__(self):
        return self.locus

class Genome(object):
    """
    An entire genome

    TODO: want this backed by a dictionary for random access
    """
    def __init__(self):
        self.gene_dict = {}
        self.intergenic_regions = []

    def add_annotation(self, key, mapping_dict):
        """
        Adds an annotation (For example, GO, PFAM, etc) to a genome.
        """

        # TODO: implement check to see if key already exists.
        for gene, annotation in mapping_dict:
            if gene not in genes:
                raise ValueError("Gene %s not in Genome" % gene)

            self.gene_dict[gene].annotation[key].append(annotation)

    def features(self, intergenic=False):
        for g in self.gene_dict:
            yield self.gene_dict[g]

        if intergenic:
            for i in self.intergenic_regions:
                yield i
