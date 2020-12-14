#!/usr/bin/env python3
"""
Parses the DMRG output file and saves spectral weights to files. 
"""
from DMRG_parser_functions import save_variables_simple

saved_sets, saved = save_variables_simple("get_spectral_weights", "get_spectral_weights")


if saved!=0:
	print("Saved {} energies in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 