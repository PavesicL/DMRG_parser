#!/usr/bin/env python3
"""
Parses the DMRG output file and saves all energy level occupations to files. 
THIS SCRIPT HAS A BIT OF A DIFFERENT APPROACH TO GETTING ALL FILES; WE'LL SEE WHAT WORKS BETTER IN THE LONG RUN.
"""

import os
import sys
import re
import numpy as np
import itertools
import subprocess

sys.path.insert(1, "/home/pavesic/git_repos/arc_manipulator/")
from arc_functions import readNameFile
from DMRG_parser_functions import readNameFileParsing, get_occs, get_nupndn
##################################

##################################
#get params, job name, path to save and the sweep parameter
name, paramList = readNameFile("nameFile")
regexName, _ = readNameFile("nameFile", regex=True)
savepath = readNameFileParsing("nameFile", "lanczos_get_energies", noParam=True)

outputName = "output.txt"	#subject to change - maybe allow user input here
result_dir = os.getcwd() + "/results"

#THIS SAVES THE LIST OF ALL FOLDERS IN /results TO endOfPipe
# define the ls command
ls = subprocess.Popen(["ls", "-p", result_dir],stdout=subprocess.PIPE,)
# define the grep command
grep = subprocess.Popen(["grep", "/$"],stdin=ls.stdout,stdout=subprocess.PIPE,)
# read from the end of the pipe (stdout)
endOfPipe = grep.stdout
allFolders = [i.decode("utf8") for i in endOfPipe]

#get all unique values of each parameter and save them to uniqueParamList.
uniqueParamList = [[] for i in range(len(paramList))]		
for i in range(len(uniqueParamList)):
	param = paramList[i][0]	

	for line in allFolders:
		
		folder = os.path.join(result_dir, line)
		#print(folder)
		#print(param)
		a = re.search(param+"(-*[0-9]+\.*[0-9]*)", folder)
		if a:
			paramval = float(a.group(1))
			uniqueParamList[i].append(paramval)

allParamCombinations = [[ uniqueParamList[k][i] for k in range(len(uniqueParamList))] for i in range(len(uniqueParamList[0]))]


saved=0
for ii in range(len(allParamCombinations)):

	result_file = os.getcwd() + "/results/" + name.format(*allParamCombinations[ii]) + "/" + outputName

	Es, ns = [], []
	nextLine, sector = False, False
	with open(result_file,  "r") as f:
		inSector=False
		for line in f:

			if sector and nextLine:
				#Es.append([float(i) for i in line.split()])
				Es.append(line)
				nextLine, sector = False, False

			a = re.search("Results in the sector with (\d+) particles:", line)
			b = re.search("ENERGIES:", line)
			if a:
				n = int(a.group(1))
				ns.append(n)
				sector=True
			if b:
				nextLine=True	
		
	#SAVE THE IMPURITY OCCUPATION AND THE TOTAL OCCUPATION TO THE HEADER
	head = "lines from sectors with n: "
	for n in ns:
		head += "{} ".format(n)

	with open(savepath.format(*allParamCombinations[ii]), "w") as out:
		out.writelines(head+"\n")
		for E in Es:
			out.writelines(E)

	#np.savetxt(fname=savepath.format(*allParamCombinations[ii]), X=Es, delimiter="	", header=head)
	#print("Saved to {}".format(savepath.format(*allParamCombinations[ii])))
	saved+=1

#SAVE THE OCCUPANCIES TO A TEXT FILE
print("Saved {0} sets of energies.".format(saved))