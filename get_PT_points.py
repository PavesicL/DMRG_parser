#!/usr/bin/env python3
"""
Parses the output file from findPT and save the Gamma_PT to file. 
"""

import os
import sys
import re
import numpy as np
import itertools

sys.path.insert(1, "/home/pavesic/git_repos/arc_manipulator/")
from arc_functions import readNameFile
from DMRG_parser_functions import readNameFileParsing
##################################

#get params, job name, path to save and the sweep parameter
name, paramList = readNameFile("nameFile")
regexName, _ = readNameFile("nameFile", regex=True)
sweepParam, savepath = readNameFileParsing("nameFile", "get_PT")

outputName = "resultFile.txt"	#subject to change - maybe allow user input here
outputName = "output.txt"

#determine which param in paramList is the sweep one -- paramList[whichParam][0] is the parameter for which to save the energy sweeps
for i in range(len(paramList)):
	param = paramList[i][0]
	if param == sweepParam:
		whichParam = i 
		break

result_dir = os.getcwd() + "/results"

#get all unique values of each parameter and save them to uniqueParamList. Range is len(paramList)-1 as one parameter sweep is always saved in the files
uniqueParamList = [[] for i in range(len(paramList)-1)]
for i in range(len(uniqueParamList)):

	param = paramList[i][0]

	if param is not sweepParam:

		for subdir, dirs, files in os.walk(result_dir):
			for direc in dirs:	#iterate over all folders
				
				folder = os.path.join(subdir, direc)

				a = re.search(param+"([0-9]+\.*[0-9]*)", folder)
				if a:
					paramval = float(a.group(1))
					uniqueParamList[i].append(paramval)

	uniqueParamList[i] = np.unique(uniqueParamList[i])	#save only unique values			




#get GammaPT from all folders
saved=0
EList=[]
for subdir, dirs, files in os.walk(result_dir):
	for direc in dirs:	#iterate over all folders
		gammafound = False	#reset the parameter for each folder
		folder = os.path.join(subdir, direc)

		aa = re.search(regexName, folder)	#for each folder, check if it is of the job

		if aa:
			#THIS IS A DIRTY TRICK, AS regaxName ALSO MATCHES gamma(0.0/gridlog), SO ALL SUBFOLDERS ARE MATCHED ALSO
			try:
				paramVals = [float(aa.group(i+1)) for i in range(len(paramList))] 	#values of the parameters are saved in this list	
			except ValueError:
				continue

			

			tuki = os.popen("tail -n 3 {0}/output.txt".format(folder)).read().splitlines() #read the last 3 lines of the file
			for line in tuki:
				a = re.search("Phase transition at Gamma = (\d+\.*\d*)", line)
				if a:
					GammaPT = float(a.group(1))	
					gammafound = True 	

			if gammafound:		
				EList.append(paramVals+[GammaPT])	#this is now a list of values of all parameters (in the same order as in paramList), with impOcc at the end
				saved+=1			
			
			#GammaPT not found	
			else:
				print("GammaPT not found!", folder)		

			"""
			#read the output file
			with open(folder+"/"+outputName, "r") as resF:
				#save the occupancy of the ground state sector

				for line in resF:

					a = re.search("Phase transition at Gamma = (\d\.+\d+)", line)
					if a:
						GammaPT = float(a.group(1))	
						break
				try:		
					EList.append(paramVals+[GammaPT])	#this is now a list of values of all parameters (in the same order as in paramList), with impOcc at the end
					saved+=1			
				
				#GammaPT not found	
				except NameError:
					print(folder)

				"""

#split EList into lists for each unique value of the parameters given in uniqueParamList
allParamCombinations = list(itertools.product(*uniqueParamList))
allParamCombinations = [list(i) for i in allParamCombinations]	#transform all elements to lists, not tuples
#splitElist is just a reshape of EList, with splitElist[i] being the energy sweep over the sweep parameter and with other parameters given as allParamCombinations[i]
splitElist = [[] for i in range(len(allParamCombinations))]	
for i in range(len(EList)):
	for j in range(len(allParamCombinations)):

		paramVals = EList[i][:-1]	#all values in EList except the ns and Es
		paramVals = paramVals[:whichParam]+paramVals[whichParam+1:]	#take out the sweep parameter

		if paramVals == allParamCombinations[j]:
			splitElist[j].append([EList[i][whichParam], EList[i][-1]])	#append the value of the sweep parameter and impOcc

#sort the lists by the value of the sweep param
for i in range(len(splitElist)):
	splitElist[i] = sorted(splitElist[i])		


#save every splitElist[i] 
for i in range(len(splitElist)):
	if len(splitElist[i])>0:
		print("Saved to {0}".format(savepath.format(*allParamCombinations[i])))
		np.savetxt(fname=savepath.format(*allParamCombinations[i]), X=splitElist[i], delimiter="	")

if saved!=0:
	print("Saved {0} sets of occupancies.".format(saved))
	
else:
	print("No files found in given range!")