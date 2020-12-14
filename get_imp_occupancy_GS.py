#!/usr/bin/env python3
"""
Parses the DMRG output file and saves impurity occupations of the nGS sector to files. 
"""

from DMRG_parser_functions import save_variables_simple

saved, saved_sets = save_variables_simple("get_imp_occupancy_GS", "get_imp_occ_nSz", [0, 0])

if saved!=0:
	print("Saved {} occupancies in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 