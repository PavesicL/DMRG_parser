#!/usr/bin/env python3

"""
Functions used within the DMRG parser scripts. 
"""

import os
import sys
import re
import numpy as np

#####################################################

def readNameFileParsing(file, observable):
	"""
	Reads the nameFile and returns the path to which to save the files and the parameter over which to save the sweeps. 
	"""

	obsCheck = False
	
	with open(file, "r") as f:

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