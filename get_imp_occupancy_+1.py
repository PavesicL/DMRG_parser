#!/usr/bin/env python3
"""
Parses the DMRG output file and saves the impurity occupation of the ground state to files. 
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
sweepParam, savepath = readNameFileParsing("nameFile", "get_imp_occupancy_p1")

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
j=-1
for i in range(len(uniqueParamList)+1):

	param = paramList[i][0]	
	if param != sweepParam:
		j+=1

		for subdir, dirs, files in os.walk(result_dir):
			for direc in dirs:	#iterate over all folders
				
				folder = os.path.join(subdir, direc)

				a = re.search(param+"([0-9]+\.*[0-9]*)", folder)
				if a:
					paramval = float(a.group(1))
					uniqueParamList[j].append(paramval)


uniqueParamList = [np.unique(i) for i in uniqueParamList]	#save only unique values				

#get energies from all folders
saved=0
EList=[]
for subdir, dirs, files in os.walk(result_dir):
	for direc in dirs:	#iterate over all folders
		folder = os.path.join(subdir, direc)

		aa = re.search(regexName, folder)	#for each folder, check if it is of the job
		bb = re.search("gridlog", folder)	#dirty fix to not include the gridlog subfolder, generated by ARC

		if aa and not bb:
			#THIS IS A DIRTY TRICK, AS regaxName ALSO MATCHES gamma(0.0/gridlog) FOR SOME REASON, SO ALL SUBFOLDERS ARE MATCHED ALSO
			try:
				paramVals = [float(aa.group(i+1)) for i in range(len(paramList))] 	#values of the parameters are saved in this list	
			except ValueError:
				continue


			#tuki = os.popen("tail -n 20 {0}/output.txt".format(folder)).read().splitlines() #read the last 20 lines of the file

			#open the output file
			with open(folder+"/"+outputName, "r") as resF:

				#find the impurity index
				for line in resF:
					a = re.search("impindex=(\d+)", line)
					if a:
						impindex = int(a.group(1))
						break
			
			with open(folder+"/"+outputName, "r") as resF:
				#find the ground state sector
				E_n = []
				for line in resF:
					b = re.search("n = ([0-9]+).*E.*= (-*[0-9]+.[0-9]*)", line) #this is specifically for energy in the standard DMRG output
					if b:
						n = int(b.group(1))
						E = float(b.group(2))
						
						E_n.append([E, n])		
				
				E_n = sorted(E_n)	#sort by E
				nGS = E_n[0][1]	

			with open(folder+"/"+outputName, "r") as resF:
				#save the occupancy of the ground state sector
				sector=False
				for line in resF:

					a = re.search("RESULTS FOR THE SECTOR WITH {0} PARTICLES".format(nGS+1), line)
					if a:
						sector=True

					if sector:
						b = re.search("site occupancies", line)
						if b:
							occupancies = re.findall("\d*\.\d+|\d+", line)

							occupancies = [float(i) for i in occupancies]

							sector=False
							
							impOcc = occupancies[impindex-1]
							break

				EList.append(paramVals+[impOcc])	#this is now a list of values of all parameters (in the same order as in paramList), with impOcc at the end
				saved+=1			

#split EList into lists for each unique value of the parameters given in uniqueParamList
allParamCombinations = list(itertools.product(*uniqueParamList))
allParamCombinations = [list(i) for i in allParamCombinations]	#transform all elements to lists, not tuples
#splitElist is just a reshape of EList, with splitElist[i] being the energy sweep over the sweep parameter and with other parameters given as allParamCombinations[i]
print(allParamCombinations)
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