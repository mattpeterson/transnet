#!/usr/bin/env python

import sys

def write_summary_html(peaks, out_dir):
    """
    Creates a summary HTML file for an experiment.

    :param peaks: list of peaks
    :type peaks: list
    """
    pass
    #for peak in peaks:

def write_text_summary(peaks, annotation, out_handle=sys.stdout):
    """
    Writes a text summary of the ChIP-seq experiment.  Produces a tab-delimited
    file with the following columns:

    Chromosome - the chromosome where the peak was found
    Start - Start position of the peak on the chromosome
    Stop - Stop position of the peak on the chromosome
    Height- the "score" of the peak
    Upstream - A list of any genes this peak was found to be upstream of
    Downstream - A list of any genes this peak was found to be downstream of
    Genic - A list of genes that this peak was found inside

    :param peaks: a list of peaks
    :param annotation: A Genome object, used to call genes
    """
    out_handle.write("Chromosome\tStart\tStop\tHeight\tUpstream\t" +
                     "Downstream\tGenic\n")

    for peak in peaks:
        (upstream, downstream, genic) = peak.get_regulated_genes(annotation)

        genic_string = ",".join([g.locus for g in genic])
        us_string = ",".join([g.locus for g in upstream])
        ds_string = ",".join([g.locus for g in downstream])

        out_handle.write("%s\t%d\t%d\t%f\t%s\t%s\t%s\n" % (peak.chromosome,
                                                           peak.chrom_start,
                                                           peak.chrom_end,
                                                           peak.score(),
                                                           us_string,
                                                           ds_string,
                                                           genic_string))

