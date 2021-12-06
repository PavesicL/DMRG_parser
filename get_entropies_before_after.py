#!/usr/bin/env python3

"""
Parses the DMRG output file and saves entropies to files. 
"""

from DMRG_parser_functions import save_variables_simple

for which in ["before", "after"]:
	print(f"saving {which}")

	saved_sets, saved = save_variables_simple("get_entropies_before_after", "get_entropies_before_after", which)

	if saved!=0:
		print("Saved {} {} entropies in {} sets.".format(saved, which, saved_sets))
		
	else:
		print("No files found in given range!") 