#!/usr/bin/env python3
"""
Parses the DMRG output file and saves the doublet overlaps. 
"""

from DMRG_parser_functions import save_variables_simple

for which in ["BCSL", "BCSR", "OS"]:
	print("saving " + which)
	saved, saved_sets = save_variables_simple("get_ZBA_overlaps_" + which, "get_ZBA_overlaps_" + which)

	if saved!=0:
		print("Saved {} overlaps in {} sets.".format(saved, saved_sets))
		
	else:
		print("No files found in given range!") 