#!/usr/bin/env python3
"""
Parses the DMRG output file and saves energies to files. 
"""

import os
import sys
import re
import numpy as np
import itertools

#import the arc_functions module
sys.path.insert(1, "/home/pavesic/git_repos/arc_manipulator/")
from arc_functions import readNameFile
from DMRG_parser_functions import readNameFileParsing

##################################

#get params, job name, path to save and the sweep parameter
name, paramList = readNameFile("nameFile")
regexName, _ = readNameFile("nameFile", regex=True)
sweepParam, savepath = readNameFileParsing("nameFile", "get_energies")

outputName = "output.txt"	#subject to change - maybe allow user input here

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




#get energies from all folders
saved=0
EList=[]
for subdir, dirs, files in os.walk(result_dir):
	for direc in dirs:	#iterate over all folders
		folder = os.path.join(subdir, direc)

		aa = re.search(regexName, folder)	#for each folder, check if it is of the job

		if aa:
			#THIS IS A DIRTY TRICK, AS regaxName ALSO MATCHES gamma(0.0/gridlog) FOR SOME REASON, SO ALL SUBFOLDERS ARE MATCHED ALSO
			try:
				paramVals = [float(aa.group(i+1)) for i in range(len(paramList))] 	#values of the parameters are saved in this list	
			except ValueError:
				continue


			#tuki = os.popen("tail -n 20 {0}/output.txt".format(folder)).read().splitlines() #read the last 20 lines of the file

			#open the output file
			with open(folder+"/"+outputName, "r") as resF:

				n_E = []
				for line in resF:
				
					a = re.search("n = ([0-9]+).*E.*= (-*[0-9]+.[0-9]*)", line) #this is specifically for energy in the standard DMRG output
					if a:
						n = int(a.group(1))
						E = float(a.group(2))
						
						n_E.append([n, E])		
				
				n_E = sorted(n_E)	#sort by ns
				Es = np.transpose(n_E)[1]
				ns = np.transpose(n_E)[0]
			
				EList.append(paramVals+[ns, Es])	#this is now a list of values of all parameters (in the same order as in paramList), a list of ns and a list of Es
				saved+=1			

#split EList into lists for each unique value of the parameters given in uniqueParamList
allParamCombinations = list(itertools.product(*uniqueParamList))
allParamCombinations = [list(i) for i in allParamCombinations]	#transform all elements to lists, not tuples
#splitElist is just a reshape of EList, with splitElist[i] being the energy sweep over the sweep parameter and with other parameters given as allParamCombinations[i]
splitElist = [[] for i in range(len(allParamCombinations))]	
for i in range(len(EList)):
	for j in range(len(allParamCombinations)):

		paramVals = EList[i][:-2]	#all values in EList except the ns and Es
		paramVals = paramVals[:whichParam]+paramVals[whichParam+1:]	#take out the sweep parameter

		if paramVals == allParamCombinations[j]:
			splitElist[j].append([EList[i][whichParam], EList[i][-2], EList[i][-1]])	#append the value of the sweep parameter, ns and Es

#sort the lists by the value of the sweep param
for i in range(len(splitElist)):
	splitElist[i] = sorted(splitElist[i])		

#save the energies to a file
for i in range(len(splitElist)):

	if len(splitElist[i])>0:
		ns = splitElist[i][0][1]	#this assumes that ns do not change during the sweep - will not work if the sweep parameter is n0

		if sweepParam == "n0":
			print("ns will change during the sweep, this will not work correctly!")

		nsstring = "	".join([str(int(i)) for i in ns])
		head = "gamma	" + nsstring

		#create and save a list with elements: [sweepParam, E1, E2, ...]
		sweepEList = []

		for j in range(len(splitElist[i])):
			sweepEList.append([splitElist[i][j][0]] + splitElist[i][j][2].tolist())

		np.savetxt(fname=savepath.format(*allParamCombinations[i]), X=sweepEList, delimiter="	", header=head)

if saved!=0:
	print("Saved {0} sets of energies.".format(saved))
	
else:
	print("No files found in given range!")