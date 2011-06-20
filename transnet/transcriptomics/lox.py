#!/usr/bin/env python
"""
Provides classes for working with the output from LOX

Zhang, Z., F. Lopez-Giraldez and J.P. Townsend, 2010. 
LOX: inferring Level Of eXpression from diverse methods of census sequencing. 
Bioinformatics 26: 1918-1919.
"""

def read(lox_handle, read_pvals = False):
    """Reads in a LOX output.
    
    :param lox_handle: A handle to the LOX output
    :type lox_handle: file
    :param read_pvals: Read the p-values for differential expression from the
                       directory the LOX output is in.  If True, will look for
                       the files in the directory, and read them into the
                       pvals attribute
    """
    experiment = LOXExperiment(lox_handle)
    
    # TODO: insert pval logic here
    
    return experiment

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
    Describes a LOX experiment.  That is, the results from a run of LOX.  This
    is dependent on the input set.
    """
    def __init__(self, lox_handle):
        self.measurements = []
        self.experiments = []
        self._parse_first_line(lox_handle)
        self._parse_lox_output(lox_handle)
    
    def _parse_lox_output(self, handle):
        """
        Reads the rest of the values from the LOX file and loads in the
        measurement data
        
        param: handle: The file handle to read from
        type: handle: file 
        """
        for line in handle:
            tokens = line.rstrip("\r\n").split("\t")
            # Create a new measurement for each 
            num_experiments = len(self.experiments)
            for i in range(0, num_experiments):
                value = float(tokens[2 + i])
                lower_bound = float(tokens[2 + i + num_experiments])
                upper_bound = float(tokens[2 + i + 2*num_experiments])
                measurement = LOXMeasurement(value, lower_bound, upper_bound)
                self.measurements[i][tokens[0]] = measurement

    def _parse_first_line(self, handle):
        """
        Gets all of the information from the header line.  Used to read
        the data, and initializes the measurements and indices 
        """
        tokens = handle.readline().rstrip("\r\n").split("\t")
        num_experiments = len(tokens[2:]) / 3
        
        for i in range(2, 2 + num_experiments):
            self.measurements.append({})
            self.experiments.append(tokens[i])
    
    def __getitem__(self, key):
        """
        Gets measurements for a given experiment.
        
        :param key: The experiment to get
        :type key: string
        """
        key_idx = self.experiments.index(key)
        return self.measurements[key_idx]
