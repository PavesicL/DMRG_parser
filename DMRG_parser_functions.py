#!/usr/bin/env python3

"""
Functions used within the DMRG parser scripts. 
"""

import os
import sys
import re
import numpy as np

#####################################################

def readNameFileParsing(file, observable, noParam=False):
	"""
	Reads the nameFile and returns the path to which to save the files and the parameter over which to save the sweeps.
	If noParam, there is no parameter to sweep over, so no param is defined (for example get_all_occupancies, saves occs
	for one file only). 
	"""

	obsCheck = False
	
	with open(file, "r") as f:
		
		if noParam:
			for line in f:

				a = re.search(observable+"{", line)
				b = re.search("}"+observable, line)

				c = re.search("path\s*=\s*(.*)", line)
		
				if line[0] == "#":	#this line is a comment
					pass
				
				if a:
					obsCheck=True
				if b:
					obsCheck=False
				
				if c and obsCheck:
					path = c.group(1)
		
			return path			
		
		else:
			for line in f:

				a = re.search(observable+"{", line)
				b = re.search("}"+observable, line)

				c = re.search("path\s*=\s*(.*)", line)
				d = re.search("sweep\s*(.*)", line)

				if line[0] == "#":	#this line is a comment
					pass
				
				if a:
					obsCheck=True
				if b:
					obsCheck=False
				
				if c and obsCheck:
					path = c.group(1)
				if d and obsCheck:
					param = d.group(1)

		return param, path			


def get_occs(n, result_file):
	"""
	Gets level occupancies, impurity occupancy and energy for a sector with n particles.
	"""

	with open(result_file,  "r") as resF:
		sector=False
		for line in resF:
			a = re.search("impindex=(\d+)", line)
			b = re.search("RESULTS FOR THE SECTOR WITH {0} PARTICLES:".format(n), line)
			c = re.search("site occupancies", line)
			d = re.search("Ground state energy = (\-*\d+\.*\d*)", line)

			if a:
				impindex = int(a.group(1))

			if b:
				sector=True

			if sector and d:
				E = float(d.group(1))	

			if sector and c:

				occupancies = re.findall("\d*\.\d+|\d+", line)

				occupancies = [float(i) for i in occupancies]
				occ = occupancies[:impindex-1] + occupancies[impindex:]
				impocc = occupancies[impindex-1]

				sector=False
				break

	return occ, impocc, E