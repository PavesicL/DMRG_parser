#!/usr/bin/env python3
"""
Reads the output of the DMRG calculation and saves files of [n0, PTgamma].
"""

import sys
import os
import shutil
import re
import math

# SET THE NAME OF THE JOBS AND THE PATH TO SAVE DATA ##############################################
jobname = "U0.1_Ec([0-9]+\.*[0-9]*)_n0([0-9]+\.*[0-9]*)_PT$"
pathtosave = 
###################################################################################################

def close(num1, num2, prec=1e-6):

	if np.abs(num1 - num2)<prec:
		return True
	else:
		return False	


print("Getting all energies...")

result_dir = os.getcwd() + "/results"


N=800

#get all Ecs
EcList=[]
for subdir, dirs, files in os.walk(result_dir):
	for direc in dirs:	#iterate over all folders

		folder = os.path.join(subdir, direc)

		aa = re.search("Ec([0-9]+\.*[0-9]*)", folder)

		if aa:
			Ec = float(aa.group(1))
			EcList.append(Ec)
EcList = np.unique(EcList)

#get all n0s
n0List=[]
for subdir, dirs, files in os.walk(result_dir):
	for direc in dirs:	#iterate over all folders

		folder = os.path.join(subdir, direc)

		aa = re.search("n0([0-9]+\.*[0-9]*)", folder)

		if aa:
			n0 = float(aa.group(1))
			n0List.append(n0)
n0List = np.unique(n0List)

#get energies from all folders
saved=0
PTlist=[]
for subdir, dirs, files in os.walk(result_dir):
	for direc in dirs:	#iterate over all folders
		folder = os.path.join(subdir, direc)

		aa = re.search("U0.1_Ec([0-9]+\.*[0-9]*)_n0([0-9]+\.*[0-9]*)_PT$", folder)

		if aa:
			Ec = float(aa.group(1))
			n0 = float(aa.group(2))

			tuki = os.popen("tail -n 20 {0}/output.txt".format(folder)).read().splitlines()

			for line in tuki:

				a = re.search("Phase transition at Gamma = ([0-9]+\.*[0-9]*)", line)
				if a:
					gamma = float(a.group(1))	
	
			PTlist.append([Ec, n0, gamma])
			saved+=1			




#PTlist is a list of triplets [Ec, n0, gamma]. Split it into a list of lists of [n0, gamma] for each Ec
PTlist = sorted(PTlist)	#sort Elist by Ec
#check to which Ec the list corresponds and append [n0, gamma] there.
splitPTlist = [[] for i in range(len(EcList))]
for i in range(len(PTlist)):
	for k in range(len(EcList)):

		if close(PTlist[i][0], EcList[k]):
			splitPTlist[k].append([PTlist[i][1], PTlist[i][2]]) 

#sort each list by n0
for i in range(len(splitPTlist)):
	splitPTlist[i] = sorted(splitPTlist[i])


#save to file
count=0
for i in range(len(EcList)):

	np.savetxt(fname=pathtosave, X=splitPTlist[i], delimiter="	")
	
	np.savetxt(fname="../../data_DMRG/paper1/fig2/PTpoints_n0_U0.1_N{0}_Ec{1}.txt".format(N, EcList[i]), X=splitPTlist[i], delimiter="	")

		
	count += 1

			
if saved!=0:
	print("Saved {0} sets of energies, for {1} Ecs.".format(saved, count))
	
else:
	print("No files found in given range!")

