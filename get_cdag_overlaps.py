#!/usr/bin/env python3
"""
Parses the DMRG output file and saves impurity occupations of the nGS sector to files. 
"""

from DMRG_parser_functions import save_variables_simple

import sys

if len(sys.argv) != 7:
	print("Give the quantum numbers of the states!")
	print("Usage: " + str(sys.argv[0]) + " n1 Sz1 i1 n2 Sz2 i2")
	exit()

n1 = int(sys.argv[1])
Sz1 = float(sys.argv[2])
i1 = int(sys.argv[3])
n2 = int(sys.argv[4])
Sz2 = float(sys.argv[5])
i2 = int(sys.argv[6])

saved_sets, saved = save_variables_simple(f"get_cdag_overlaps", "get_cdag_overlaps", [n1, Sz1, i1, n2, Sz2, i2])

if saved!=0:
	print("Saved {} cdag overlaps in {} sets.".format(saved, saved_sets))	
else:
	print("No files found in given range!") 