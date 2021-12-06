#!/usr/bin/env python3
"""
Parses the DMRG output file and saves energies to files. 
"""

from DMRG_parser_functions import save_variables_simple

for which in ["u", "v"]:
	print("saving " + which)
	saved_sets, saved = save_variables_simple("get_sc_amplitudes", "get_sc_amplitudes", which)

	if saved!=0:
		print("Saved {} amplitudes in {} sets.".format(saved, saved_sets))
		
	else:
		print("No files found in given range!") 