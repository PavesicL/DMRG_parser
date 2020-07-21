#!/usr/bin/env python3
"""
Parses the DMRG output file and saves the dimensions of matrices during sweeps, saves data for each n-sector to a separate file. 
"""

import os
import sys
import re
import numpy as np
import itertools

sys.path.insert(1, "/home/pavesic/git_repos/arc_manipulator/")
from arc_functions import readNameFile
from DMRG_parser_functions import readNameFileParsing, get_occs
##################################

#get params, job name, path to save and the sweep parameter
name, paramList = readNameFile("nameFile")
regexName, _ = readNameFile("nameFile", regex=True)
savepath = readNameFileParsing("nameFile", "get_all_dimensions", noParam=True)

outputName = "output"	#subject to change - maybe allow user input here

#get the input values and save them to paramVals
if len(sys.argv) != 1 + len(paramList):
	printstring = ""
	for i in paramList:
		printstring += i[0] + " "
	print("usage: " + sys.argv[0] + " " + printstring)
	exit()

else:
	paramVals = sys.argv[1:]


result_file = os.getcwd() + "/results/" + name.format(*paramVals) + "/" + outputName

if not os.path.isfile(result_file):
	result_file += ".txt"
if not os.path.isfile(result_file):
	print("The output file is not output or output.txt; or does not exist!")


#OPEN THE OUTPUT FILE AND FIND THE GROUND STATE SECTOR
with open(result_file,  "r") as resF:
	
	dimlist=[]
	ns=[]

	inSector=False	
	sweep, HS = 0, 0
	for line in resF:
		
		if inSector: 
			oldsweep=sweep

		a = re.search("Sweeping in the sector with (\d+) particles.", line)	#to determine the sector
		b = re.search("Sweep=([0-9]+), HS=([1-2]), Bond=([0-9]+)/[0-9]+", line)										
		c = re.search("States kept: dim=([0-9]+)", line)

		if a:
			n = int(a.group(1))
			ns.append(n)		#save ns to this list
			inSector=True
			dimlist.append([])	#append a new empty list, save the dimensions to it

		if b and inSector:
			sweep = int(b.group(1))
			HS = int(b.group(2))
			bond = int(b.group(3))
			if sweep!=oldsweep:
				dimlist[-1].append([])

		#if the line contains info about dimension, append it to the correct sub-list
		if c and inSector:
			dimlist[-1][sweep-1].append(c.group(1))


for i in range(len(dimlist)):

	with open(savepath.format(*paramVals, ns[i]), "w+") as ff:
		for j in range(len(dimlist[i][0])):
			dimstring = "{0}	".format(j)
			for k in range(len(dimlist[i])):
				dimstring += str(dimlist[i][k][j])
				dimstring += "	"
			dimstring+="\n"	

			ff.writelines(dimstring)

	print("Saved to {}".format(savepath.format(*paramVals, ns[i])))
