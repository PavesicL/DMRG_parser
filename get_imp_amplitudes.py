#!/usr/bin/env python3
"""
Parses the DMRG output file and saves impurity occupations of the nGS sector to files. 
"""

from DMRG_parser_functions import save_variables_simple

for which in ["zero", "up", "down", "two"]:
	print("saving " + which)
	saved, saved_sets = save_variables_simple("get_imp_amplitudes_" + which, "get_imp_amplitudes_" + which)

	if saved!=0:
		print("Saved {} occupancies in {} sets.".format(saved, saved_sets))
		
	else:
		print("No files found in given range!") 