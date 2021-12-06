#!/usr/bin/env python3

"""
Parses the DMRG output file and saves entropies to files. 
"""

from DMRG_parser_functions import save_variables_simple

saved_sets, saved = save_variables_simple("get_entropies", "get_entropies")

if saved!=0:
	print("Saved {} entropies in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 