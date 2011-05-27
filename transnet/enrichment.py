# Tests for enrichment
__author__ = "Matthew Peterson"

def read_mapping(handle):
    """
    Reads a mapping from genes to
    """
    handle = iter(handle)
    for line in handle:
        tokens = line.rstrip("\r\n").split("\t")
        

class EnrichmentTest(object):
    def _read_mapping(handle):
        pass
