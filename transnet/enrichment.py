# Tests for enrichment
__author__ = "Matthew Peterson"

def read_mapping(handle, key="function"):
    """
    Reads a mapping from genes to 
    """
    features = {}
    
    handle = iter(handle)
    for line in handle:
        tokens = line.rstrip("\r\n").split("\t")
        features[tokens[0]] = tokens[1]
    
    return features

class EnrichmentTest(object):
    def _read_mapping(handle):
        pass
