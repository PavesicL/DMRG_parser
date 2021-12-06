#!/usr/bin/env python3


from DMRG_parser_functions import save_variables_simple

saved_sets, saved = save_variables_simple("get_right_spin_spin_sum", "get_right_spin_spin_sum")

if saved!=0:
	print("Saved {} occupancies in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 