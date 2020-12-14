#!/usr/bin/env python3
"""
Parses the DMRG output file and saves ns to files. This is useful in the cases where the range of computed n-sectors changes during the iteration. 
"""
from DMRG_parser_functions import save_variables_simple

saved_sets, saved = save_variables_simple("get_exciteds", "get_exciteds")

if saved!=0:
	print("Saved {} exciteds in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 