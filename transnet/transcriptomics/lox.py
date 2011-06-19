#!/usr/bin/env python
"""
Provides classes for working with the output from LOX

Zhang, Z., F. Lopez-Giraldez and J.P. Townsend, 2010. 
LOX: inferring Level Of eXpression from diverse methods of census sequencing. 
Bioinformatics 26: 1918-1919.
"""

def read(lox_handle, read_pvals = False):
    pass

class LOXMeasurement(object):
    
    _scale_value = None
    
    def __init__(self, level, lower_confidence, upper_confidence):
        self.expression_level = level
        self.lower_confidence = lower_confidence
        self.upper_confidence = upper_confidence
    
    def scale(self, multiplier):
        """
        Scales an experiment by an absolute expression level.  Can be useful
        for presenting
        """
        self.expression_level = self.expression_level * multiplier
        self.lower_confidence = self.lower_confidence * multiplier
        self.upper_confidence = self.upper_confidence * multiplier

class LOXExperiment(object):
    """
    Describes a LOX experiment
    """
    def __init__(self, lox_handle):
        self.measurements = []
        self.indices = []

    def _parse_first_line(self, handle):
        """
        Gets all of the information from the header line.  Used to read
        the data, and initializes the measurements and indices 
        """
        tokens = handle.readline().rstrip("\r\n").split("\t")
        num_experiments = len(tokens[2:]) / 3
        
        for i in range(2, 2 + num_experiments):
            self.measurements.append({})
            self.indices.append(tokens[i])
    
    def __getitem__(self, key):
        key_idx = 
    
    def experiments(self):
        return self.measurements.keys()
