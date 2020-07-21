#!/usr/bin/env python3
"""
Parses the DMRG output file and saves all energy level occupations to files. 
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
savepath = readNameFileParsing("nameFile", "get_one_all_occupancies", noParam=True)

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
	E_n = []
	for line in resF:
		b = re.search("n = ([0-9]+).*E.*= (-*[0-9]+.[0-9]*)", line) #this is specifically for energy in the standard DMRG output
		if b:
			n = int(b.group(1))
			E = float(b.group(2))
			
			E_n.append([E, n])		
	
	E_n = sorted(E_n)	#sort by E
	nGS = E_n[0][1]	

print("ENERGIES")
print(E_n)

#OPEN THE OUTPUT FILE AND FIND THE impIndex AND SAVE THE OCCUPANCIES FOR THE GROUND STATE AND THE TWO SPECTROSCOPICALLY AVAILABLE STATES
print("nGS is ", nGS)
occ, impocc, Egs = get_occs(nGS, result_file)
occp, impoccp, Ep = get_occs(nGS+1, result_file)
occm, impoccm, Em = get_occs(nGS-1, result_file)

#CREATE THE DATA TO SAVE
toSave = np.array([occ, occp, occm])
toSave = np.transpose(toSave)

#SAVE THE IMPPURITY OCCUPATION AND THE TOTAL OCCUPATION TO THE HEADER
headSAMPLE = "#n = {0}, impocc = {1}, nSC = {2}, E={3}\n"
headGS = headSAMPLE.format(nGS, impocc, sum(occ), Egs)
headp = headSAMPLE.format(nGS+1, impoccp, sum(occp), Ep)
headm = headSAMPLE.format(nGS-1, impoccm, sum(occm), Em)

head = headGS + headp + headm

#SAVE THE OCCUPANCIES TO A TEXT FILE
print("The ground state has {0} particles.".format(nGS))
print(head)
np.savetxt(fname=savepath.format(*paramVals), X=toSave, delimiter="	", header=head)
print("Saved to {}".format(savepath.format(*paramVals)))