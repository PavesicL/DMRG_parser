#!/usr/bin/env python3
"""
Parses the DMRG output file and saves impurity occupations of the nGS sector to files. 
"""

from DMRG_parser_functions import save_variables_simple

import sys

if len(sys.argv) != 3:
	print("Set which matrix element you want!")
	print("Usage: " + str(sys.argv[0]) + " i, j")
	exit()

i = int(sys.argv[1])
j = int(sys.argv[2])


saved, saved_sets = save_variables_simple(f"get_charge_susceptibility_{i}_{j}", "get_charge_susceptibility", [i, j])

if saved!=0:
	print("Saved {} occupancies in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 