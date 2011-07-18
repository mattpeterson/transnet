#!/usr/bin/env python
"""
Provides classes for working with the output from LOX

Zhang, Z., F. Lopez-Giraldez and J.P. Townsend, 2010. 
LOX: inferring Level Of eXpression from diverse methods of census sequencing. 
Bioinformatics 26: 1918-1919.
"""

from os import path
import math
import re

def read(lox_file, read_pvals = False):
    """Reads in a LOX output.
    
    :param lox_handle: A handle to the LOX output
    :type lox_handle: file
    :param read_pvals: Read the p-values for differential expression from the
                       directory the LOX output is in.  If True, will look for
                       the files in the directory, and read them into the
                       pvals attribute
    """
    with open(lox_file) as lox_handle:
        experiment = LOXExperiment(lox_handle)
    
    results_dir = path.dirname(lox_file)
    if results_dir == "":
        results_dir = "."
    
    filename = path.basename(lox_file)
    (filename, extension) = path.splitext(filename)
    
    if read_pvals:
        # Iterate through each experiment, and load in the p-values
        for e in experiment.experiments:
            experiment.add_pvals(e, results_dir + "/" + filename + "." + e + 
                                 ".pvalue")

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
        for presenting values in terms of approximate absolute expression,
        rather than using relative expression
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
        self.pvalues = {}
        self._parse_first_line(lox_handle)
        self._read_lox_output(lox_handle)
    
    def _read_lox_output(self, handle):
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
    
    def scale(self, rpkm_dict):
        """ 
        Scales the set of experiments in order to get some measure of
        absolute expression.
        """
        pass
    
    def add_pvals(self, experiment, filename):
        """
        Adds a set of pvalues to the experiment.  Reads in the corresponding
        .pvalues file from LOX, and adds them to the 
        """
        with open(filename) as filehandle:
            pairs = self._parse_pval_header(filehandle.readline())
            for line in filehandle:
                tokens = line.rstrip("\r\n").split("\t")
                locus = tokens[0]
                for i, t in enumerate(tokens[2:]):
                    value = float(t)
                    (from_exp, to_exp) = pairs[i]
                    self.pvalues[(from_exp, to_exp, locus)] = value
    
    def get_diff_expressed(self, experiment_1, experiment_2, 
                           min_change, threshold):
        """
        Finds differentially expressed genes between two experiments.
        """
        diff_expressed = set()
        
        # Get the indices for experiments first, so we don't have to do the 
        # key search each time
        idx_1 = self.experiments.index(experiment_1)
        idx_2 = self.experiments.index(experiment_2)
        
        genes = self.measurements[0].keys()
        if experiment_1 not in self.experiments or \
           experiment_2 not in self.experiments:
            raise KeyError("Experiment not in dataset.")
        
        for g in genes:
            
            # TODO: Let's put this in a get_ratio function
            val_exp_1 = self.measurements[idx_1][g].expression_level
            if val_exp_1 == 0:
                val_exp_1 = 0.01
            
            val_exp_2 = self.measurements[idx_2][g].expression_level
            if val_exp_2 == 0:
                val_exp_2 = 0.01
            
            ratio = math.fabs(math.log(val_exp_1 / val_exp_2, 2))
            
            pval = self.pvalues[(experiment_1, experiment_2, g)]
            if pval <= threshold or (1.0 - pval) <= threshold:
                if ratio > min_change:
                    diff_expressed.add((g, pval, ratio))
        
        return diff_expressed
    
    def get_pval(experiment_1, experiment_2, gene):
        pass
    
    def _parse_pval_header(self, line):
        """
        Reads in the header from the LOX p-value file and makes sure the names
        match what is expected in the file.
        """
        pairs = []
        
        # RegEx pattern for matching each header.
        pattern = re.compile('^P\((.+)>(.+)\)$')
        
        tokens = line.rstrip("\r\n").split("\t")
        
        for t in tokens[2:]:
            pairs.append(pattern.match(t).groups())
        
        return pairs
    
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
