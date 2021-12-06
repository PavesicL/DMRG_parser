#!/usr/bin/env python3


from DMRG_parser_functions import save_variables_simple

for which in ["1", "2"]:
	print(f"Saving {which}")

	saved_sets, saved = save_variables_simple("get_channels_energy", "get_channels_energy", which)

if saved!=0:
	print("Saved {} energies in {} sets.".format(saved, saved_sets))
	
else:
	print("No files found in given range!") 