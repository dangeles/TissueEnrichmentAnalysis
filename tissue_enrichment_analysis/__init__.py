# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:48:18 2016

Tissue Enrichment Analysis (TEA).

classes: none yet (v 0.11)

functions:
enrichment_analysis
plot_enrichment_results

@author: davidangeles
"""

__author__ = "David Angeles-Albores"
__copyright__ = "Copyright 2016, WormBase"
__credits__ = ["David Angeles-Albores", "Raymond Y. Lee", "Juancarlos Chan"]
__license__ = "MIT"
__version__ = "0.13.13"
__maintainer__ = "David Angeles-Albores"
__email__ = "dangeles@caltech.edu"
__status__ = "Production"

from .hypergeometricTests import enrichment_analysis
from .hypergeometricTests import plot_enrichment_results
from .hypergeometricTests import fetch_dictionary
