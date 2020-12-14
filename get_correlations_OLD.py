#!/usr/bin/env python3
"""
Parses the DMRG output file and saves charge correlations for each level to files. 
"""

import os
import sys
import re
import numpy as np
import itertools
import subprocess


sys.path.insert(1, "/home/pavesic/git_repos/arc_manipulator/")
from arc_functions import readNameFile
from DMRG_parser_functions import readNameFileParsing, get_occs
##################################
#Choose which correlations to save
#The input argument should be a combination of letters representing the correlations one wants to save. 
#eg. to save charge and pairing correlations, the argument should be cp or pc. To save all, a or csp should both work.

if len(sys.argv)<2:
	print("usage: " + sys.argv[0] + " a(all) OR c(charge) s(spin) p(pairing) h(hopping)")
	exit()
else:
	chargeCorrelationSave, spinCorrelationSave, pairCorrelationSave, hoppingCorrelationSave = False, False, False, False
	if re.search("c", sys.argv[1]) or sys.argv[1]=="a":
		chargeCorrelationSave=True
	if re.search("s", sys.argv[1]) or sys.argv[1]=="a":
		spinCorrelationSave=True
	if re.search("p", sys.argv[1]) or sys.argv[1]=="a":
		pairCorrelationSave=True
	if re.search("h", sys.argv[1]) or sys.argv[1]=="a":
		hoppingCorrelationSave=True

##################################
#get params, job name, path to save and the sweep parameter
name, paramList = readNameFile("nameFile")
regexName, _ = readNameFile("nameFile", regex=True)
savepath = readNameFileParsing("nameFile", "get_correlations", noParam=True)

outputName = "output"	#subject to change - maybe allow user input here
result_dir = os.getcwd() + "/results"

"""
#get all unique values of each parameter and save them to uniqueParamList. Range is len(paramList)-1 as one parameter sweep is always saved in the files
uniqueParamList = [[] for i in range(len(paramList))]
for i in range(len(uniqueParamList)):
	param = paramList[i][0]

	for subdir, dirs, files in os.walk(result_dir):
		for direc in dirs:	#iterate over all folders
			
			folder = os.path.join(subdir, direc)

			a = re.search(param+"([0-9]+\.*[0-9]*)", folder)
			if a:
				paramval = float(a.group(1))
				uniqueParamList[i].append(paramval)


uniqueParamList = [np.unique(i) for i in uniqueParamList]	#save only unique values	
allParamCombinations = list(itertools.product(*uniqueParamList))
allParamCombinations = [list(i) for i in allParamCombinations]	#transform all elements to lists, not tuples
"""
#THIS SAVES THE LIST OF ALL FOLDERS IN /results TO endOfPipe
# define the ls command
ls = subprocess.Popen(["ls", "-p", result_dir],stdout=subprocess.PIPE,)
# define the grep command
grep = subprocess.Popen(["grep", "/$"],stdin=ls.stdout,stdout=subprocess.PIPE,)
# read from the end of the pipe (stdout)
endOfPipe = grep.stdout
allFolders = [i.decode("utf8") for i in endOfPipe]

#get all unique values of each parameter and save them to uniqueParamList. Range is len(paramList)-1 as one parameter sweep is always saved in the files
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


#For each unique combination of parameters, save the correlations for each n
saved=0
for ii in range(len(allParamCombinations)):

	result_file = os.getcwd() + "/results/" + name.format(*allParamCombinations[ii]) + "/" + outputName

	#CHECK IF THE OUTPUT FILE EXISTS
	if not os.path.isfile(result_file):
		result_file += ".txt"
	if not os.path.isfile(result_file):
		print("The output file is not output or output.txt; or does not exist!")

	ns=[]	#a list of ns in the same order as the correlations are saved
	chargeCorr, spinCorr, pairCorr, hopCorrUp, hopCorrDn = [], [], [], [], []
	with open(result_file, "r") as f:
		inSector=False
		spinCorrs=False
		for line in f:

			a = re.search("RESULTS FOR THE SECTOR WITH (\d+) PARTICLES", line)
			b = re.search("charge correlation", line)
			bb = re.search("charge correlation tot", line)

			c = re.search("pair correlation", line)
			cc = re.search("pair correlation tot", line)
			
			d = re.search("spin correlations:", line)
			e = re.search("SzSz correlations:", line)
			f = re.search("S\+S- correlations:", line)
			g = re.search("S-S\+ correlations:", line)

			h = re.search("hopping spin up", line)
			l = re.search("hopping spin down", line)

			if a: 
				inSector=True	
				n=a.group(1)
				ns.append(n)

			#######################################################
			#CHARGE CORRELATIONS

			if chargeCorrelationSave and inSector and b and not bb:		
				chCorr = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				chargeCorr.append([float(i) for i in chCorr])

				if not pairCorrelationSave and not spinCorrelationSave:
					inSector=False 	#if this is the only corr being saved, turn off inSector here

			#######################################################
			#FOR SPIN IT IS A BIT MORE COMPLICATED
			if spinCorrelationSave and inSector and d:	
				spinCorrs=True

			if spinCorrelationSave and inSector and spinCorrs and e:
				SzSzCorr = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				SzSzCorr = np.array([float(i) for i in SzSzCorr])

			if spinCorrelationSave and inSector and spinCorrs and f:
				SpSmCorr = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				SpSmCorr = np.array([float(i) for i in SpSmCorr])

			if spinCorrelationSave and inSector and spinCorrs and g:
				SmSpCorr = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				SmSpCorr = np.array([float(i) for i in SmSpCorr])
				
				if not pairCorrelationSave:
					inSector=False	#if pair corrs are not being saved, turn off inSector here

				spinCorr.append(SzSzCorr + 0.5*(SpSmCorr + SmSpCorr))			

			#######################################################
			#PAIR CORRELATIONS

			if pairCorrelationSave and inSector and c and not cc:		
				pCorr = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				pairCorr.append([float(i) for i in pCorr])
				
				if not hoppingCorrelationSave:
					inSector=False	#if hopping corrs are not being saved, turn off inSector here

			#######################################################
			#HOPPING CORRELATIONS

			if hoppingCorrelationSave and inSector and h:
				hCorrUp = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				hopCorrUp.append([float(i) for i in hCorrUp])	

			if hoppingCorrelationSave and inSector and l:
				hCorrDn = re.findall("-*\d*\.\d+e[-\+]\d+|-*\d*\.\d+|-*\d+", line)
				hopCorrDn.append([float(i) for i in hCorrDn])	
				inSector=False	#this is the last correlation, turn off inSector

	head=""
	for i in range(len(ns)):
		head+=str(ns[i]) + " "				
	head=head[:-1]	

	#SAVE THE CORRELATIONS TO FILES
	if chargeCorrelationSave:
		chargeCorr = np.array(chargeCorr).transpose()
		np.savetxt(X=chargeCorr, fname=savepath.format("charge", *allParamCombinations[ii]), delimiter="	", header=head)
		
	if spinCorrelationSave:
		spinCorr = np.array(spinCorr).transpose()
		
		#print("AA")
		#print(len(spinCorr))
		#print(spinCorr)
		np.savetxt(X=spinCorr, fname=savepath.format("spin", *allParamCombinations[ii]), delimiter="	", header=head)

	if pairCorrelationSave:
		pairCorr = np.array(pairCorr).transpose()
		np.savetxt(X=pairCorr, fname=savepath.format("pair", *allParamCombinations[ii]), delimiter="	", header=head)		
	
	if hoppingCorrelationSave:
		hopCorrUp = np.array(hopCorrUp).transpose()
		hopCorrDn = np.array(hopCorrDn).transpose()

		np.savetxt(X=hopCorrUp, fname=savepath.format("hopping_up", *allParamCombinations[ii]), delimiter="	", header=head)		
		np.savetxt(X=hopCorrDn, fname=savepath.format("hopping_dn", *allParamCombinations[ii]), delimiter="	", header=head)		

	saved+=1

print("Saved {0} sets of correlations.".format(saved))