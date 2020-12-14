#!/usr/bin/env python3

"""
Functions used within the DMRG parser scripts. 
"""

import os
import sys
import re
import math
import numpy as np
import h5py
import argparse

###################################################################################################

def readNameFileParsing(file, observable, noParam=False):
	"""
	Reads the nameFile and returns the path to which to save the files and the parameter over which to save the sweeps.
	If noParam, there is no parameter to sweep over, so no param is defined (for example get_all_occupancies, saves occs
	for one file only). 
	"""

	obsCheck = False
	param = None

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
				if line[0]=="#":
					continue

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
		try: 		
			return param, path			
		except UnboundLocalError:
			return 0, 0

def readNameFile(file, regex=False):
	"""
	Reads the nameFile and returns a jobname and a list of parameters and their types.
	If regex = true, instead of a python formatted string, return the name as a regex expression string.  
	INPUT:
	file - relative path to the file with the information about jobname and parameters.
	regex - wheter to return the name as a regex expression string
	OUTPUT:
	name/regexname - a jobname, a python formatted string or a regex type string
	paramsList - a list of all parameters and their types, paramsList[i] = [param, paramtype]
	"""

	paramsCheck = False
	paramsList = []

	with open(file, "r") as f:

		for line in f:

			a = re.search("name\s*=\s*(.*)", line)	
		
			b = re.search("params\s*{", line) 
			c = re.search("}\s*endparams", line) 
			
			d = re.search("\s*(\S*)\s*(\S*)\s*(\S*)", line)
			e = re.search("\s*(\S*)\s*(\S*)", line)
			
			if line[0] == "#":	#this line is a comment
				pass
			if a:	#this line has the name, save it
				name = a.group(1)
			if b:	#we are in the part of the file with the params info
				paramsCheck=True
				continue
			if c:	#we are past the part of the file with the params info
				paramsCheck=False

			if paramsCheck and d:	#parse parameters
				param, paramtype, takeOut = d.group(1), d.group(2), d.group(3)
				paramsList.append([param, paramtype, takeOut])

			elif paramsCheck and e:	#parse parameters
				param, paramtype = e.group(1), e.group(2)
				paramsList.append([param, paramtype])


	if regex:
		regexname = re.sub("{[0-9]+}", "(-*[0-9]+\.*[0-9]*)", name)	#replace all instances of {number} in the name with ([0-9]+.*[0-9]*), which matches floats and ints
		return regexname, paramsList
	
	else:
		return name, paramsList

###################################################################################################

def save_variables_simple(which_var, which_func, params=[]):
	"""
	This is a general function. Given the name of the variable to save, it handles the reading of
	the parameters, iteration through the results and reshaping and saving to the files. 
	Parameters are given in the nameFile file.
	The actual parsing is done in the function called by which_var; defined in this file too (eg. get_energies()).
		
	quantity is an optional argument which gives the name of the quantity, and HAS TO MATCH ITS NAME IN THE HDF5 FILE. It
			is passed to the hdf5 parsing function and is used to determine the path to the number in the hdf5 file.

	Currently (OCT 2020), the first function called is the hdf5 one, if it fails, the text file parsing one is called.
	"""

	#get params, job name, path to save and the sweep parameter
	name, paramList = readNameFile("nameFile")
	regexName, _ = readNameFile("nameFile", regex=True)
	sweepParam, savepath = readNameFileParsing("nameFile", which_var)

	outputName = "output"	#subject to change - maybe allow user input here

	if name == 0 and paramList == 0:
		if not re.match("print.*", which_var):
			print("Name and paramList were not found. Check if nameFile has correct data!")
			exit()

	#determine which param in paramList is the sweep one -- paramList[whichParam][0] is the parameter for which to save the sweeps
	for i in range(len(paramList)):
		param = paramList[i][0]
		if param == sweepParam:
			whichParam = i 
			break


	#take out the parameters that are labeled as take out
	tmpl=[]
	for pp in paramList:
		if len(pp)>=3:
			if pp[2].lower()=="true" or pp[2].lower()=="takeout":
				pass
			else:
				tmpl.append(pp[:2])
		else:
			tmpl.append(pp[:2])
	paramList = tmpl		

	#if whichParam is not defined (because sweepParam == None), set it to -1 
	try:
		whichParam==1
	except:
		whichParam=-1


	result_dir = os.getcwd() + "/results"

	#LOAD finishedJobs.txt AND CHECK ONLY THOSE FOLDERS
	#finishedJobs=[]
	#with open("finishedJobs.txt", "r") as fff:
	#	for line in fff:
	#		finishedJobs.append(line[:-1])	#TAKE AWAY THE \n SYMBOL WITH [:-1]

	#get results from all folders and save them to Elist
	EList=[]
	saved=0
	noth5, h5 = 0, 0
	for subdir, dirs, files in os.walk(result_dir):
		for direc in dirs:	#iterate over all folders
			folder = os.path.join(subdir, direc)
			
			aa = re.match(regexName, direc)	#for each folder, check if it is of the job


			if aa:# and "DONE" in os.listdir(folder):	#ALSO CHECK WHETER THE FOLDER BELONGS TO A FINISHED JOB (there is a file called DONE in its folder), TO AVOID PARSING UNFINISHED OUTPUT FILES. 
													#THIS ONLY WORKS FOR SLURM, AS IN ARC ALL DOWNLOADED JOBS ARE ALREADY FINISHED. 
													#SO THIS THING NEEDS A SWITCH IF THE JOBS ARE FROM SLURM.

				#finishedJobs.remove(direc)	

				paramVals = [float(aa.group(i+1)) for i in range(len(paramList))] 	#values of the parameters are saved in this list	

				result_file_h5 = folder+"/solution.h5"

				result_file = folder+"/"+outputName
				#CHECK IF THE OUTPUT FILE EXISTS
				if not os.path.isfile(result_file):
					result_file += ".txt"
					#result_file = result_file[:-4]
				if not os.path.isfile(result_file):
					print("The output file is not output or output.txt; or does not exist! Job: "+direc)
					continue

				#TRY CALLING THE HDF5 PARSING FUNCTION. IF IT FAILS (useful if there is no .h5 file, typically in old jobs), FALLBACK TO TEXT PARSING. 
				try:
					states = get_all_states(result_file_h5)					
					this_module = sys.modules[__name__]
					Es = getattr(this_module, which_func)(result_file_h5, states)
					h5+=1
						
				except:
					if which_var != "print_times":
						#print times is not saved in the .h5 files, so the fall back will happen every time
						print("Falling back to the text parsing function.")	
						print("file: "+result_file)
						print()

					#This magic calls the function with the name which_func, with the argument result_file and *params
					#This is here in order to be able to pass a name of the function as an argument, and call it here (eg. "get_energies" is passed as a param, so in this line the
					#function get_energies() is called). 
					this_module = sys.modules[__name__]
					Es = getattr(this_module, which_func+"_text")(result_file, *params)

					noth5+=1 	#just to check how many files were saved not in the h5 way
				
				if len(Es)!=0:
					EList.append(paramVals+[Es])	#this is now a list of values of all parameters (in the same order as in paramList), and corresponding Es
					saved+=1

	#THIS SORTS EList BY ALL NON-SWEEP PARAMS FIRST, AND THEN BY THE SWEEP PARAM LAST. 
	#THE RESULT IS THE LIST OF Es, THAT CAN BE SLICED WHERE THE SET OF PARAMETERS (EXCEPT THE SWEEP ONE) CHANGES 
	EList = sorted(EList, key = lambda x : ([x[i] for i in range(whichParam)] + [x[i] for i in range(whichParam+1,len(paramList))] + [x[whichParam]]) )

	if which_var == "print_times":
		print("getting walltimes...")
		onlyTimes = [i[-1][0] for i in EList]
		np.savetxt(fname="times", X=onlyTimes)
		return None

	else:
		previousParamVals = np.array([EList[0][j] for j in range(whichParam)] + [EList[0][j] for j in range(whichParam+1,len(paramList))])
		tempList=[]
		saved_sets=0
		for i in range(len(EList)):
			currentParamVals = np.array([EList[i][j] for j in range(whichParam)] + [EList[i][j] for j in range(whichParam+1,len(paramList))])
			
			#when the parameters change, save what you accumulated in previous steps and continue with new params
			compareList = [i for i, j in zip(previousParamVals, currentParamVals) if i==j]
			
			if len(compareList)!=len(currentParamVals):
			
				#np.savetxt(fname=savepath.format(*previousParamVals), X=tempList, delimiter="	")
				with open(savepath.format(*previousParamVals), "w") as txt_file:
				    for line in tempList:
				    	txt_file.write("	".join([str(i) for i in line]) + "\n") # works with any number of elements in a line
				
				saved_sets+=1
				
				#reset the list 
				tempList = []
				
			if whichParam==-1:
				tempList.append(EList[i][-1])
			#accumulate values of the sweep parameter and of observables in a temp list
			else:
				tempList.append([EList[i][whichParam]]+list(EList[i][-1]))

			previousParamVals = currentParamVals	#overwrite the values of parameters	with the new ones

		#now save the last set:
		#print(tempList, savepath.format(*previousParamVals))
				
		#print(tempList)
		with open(savepath.format(*previousParamVals), "w") as txt_file:
			for line in tempList:
				txt_file.write("	".join([str(i) for i in line]) + "\n") # works with any number of elements in a line
					
		#np.savetxt(fname=savepath.format(*previousParamVals), X=tempList, delimiter="	")	
		saved_sets+=1
		print("{} were saved from the hdf5 file, while {} were NOT.".format(h5, noth5))

		return saved_sets, saved

###################################################################################################
#UTILITY FUNCTIONS

def read_impindex(file):
	with open(file, 'r') as f:
		for line in f:
			a = re.search("impindex=(\d+)", line)
			if a:
				impindex = int(a.group(1))

				return impindex

def get_impindex(file):
	"""
	Gets the impurity index from the text output file. This is 1-based, eg. impindex=1 is the first site!
	"""

	#Assuming the .h5 file is solution.h5, this transforms file (the path to solution.h5) to the path to file output or output.txt
	a = re.match("(.*\/).*", file)
	impFile = a.group(1) + "output"

	try:
		impind = read_impindex(impFile)
	except:
		impind = read_impindex(impFile + ".txt")

	return impind	

def get_all_quantities_h5(file, quantity):
	"""
	Opens the DMRG output and recovers the quantity, saved in the regexPath in the HDF5 file.
	Works if the quantity is a single number!
	Used for energy, 
	"""



	#READ ALL PATH IN THE FILE USING h5dump
	os.system("h5dump -n "+file+" > h5content.temp")

	hf = h5py.File(file, "r")

	#ITERATE OVER ALL CONTENTS OF THE h5 FILE, FINDING ONLY THE RELEVANT PATHS
	ns, Szs, exciteds, Es = [], [], [], []
	with open("h5content.temp", "r") as contentFile:
		for line in contentFile:
			a = re.search("dataset\s+(/(\d+)/(-?\d+.?\d*)/(\d)/"+quantity+")", line)
			if a:
				fullPath = a.group(1)
				n = int(a.group(2))	
				Sz = float(a.group(3))
				excited = float(a.group(4))

				E = np.array(hf.get(fullPath))[()]	#get the number, transform it to np.array, which is 0-d. The get the only element of the array uiong [()]
				
				ns.append(n)
				Szs.append(n)
				exciteds.append(excited)
				Es.append(E)	

	sortedEs = [E for _, _, _, E in sorted(zip(ns, Szs, exciteds, Es), key = lambda x : (x[0], x[1], x[2]))]

	os.system("rm h5content.temp")

	return sortedEs

def get_quantity_h5(file, state, quantity):
	"""
	Given the .h5 file, the state quantum numbers and the name of the quantity, returns it.
	"""
	hf = h5py.File(file, "r")

	n, Sz, excited, path = state

	#path = "/{}/{}/{}/{}/".format(n, Sz, excited, quantity)

	E = np.array(hf.get(path+"/"+quantity))

	#if E is a zero dimensional quantity, take it out of the list
	if E.shape == ():
		E = E[()]

	return E		

def get_all_states(file):
	"""
	Returns all states that are saved in the h5 file given.
	"""
	os.system("h5dump -n "+file+" > h5content.temp")

	states=[]
	with open("h5content.temp", "r") as contentFile:
		for line in contentFile:
			a = re.search("(/(\d+)/(-?\d+.?\d*)/(\d)/)", line)
			if a:
				fullPath = a.group(1)
				n = int(a.group(2))	
				Sz = float(a.group(3))
				excited = float(a.group(4))

				states.append([n, Sz, excited, fullPath])

	states = np.unique(states, axis=0) #axis=0 is important here! Without it the array is flattened.
	states = np.array(sorted(states, key = lambda x : (x[0], x[1], x[2]))) #sort this to get sorted results later
		
	os.system("rm h5content.temp")
	return states

def put_impindex_to_front(a, impindex):
	"""
	Given the list a, take the element at impindex and put it as the first element of the list.
	"""
	a = a.tolist()

	impval = a.pop(impindex-1)	#pop = delete + return
	a = [impval] + a

	a = np.array(a)

	return a

###################################################################################################
#hdf5 parsing functions
def print_times(file, states):
	
	raise Exception("Times are not saved in the h5 files!")	
	return None

def get_energies(file, states):
	"""
	Get energies from the .h5 output file.
	"""
	Es = []
	for s in states:
		Es.append(get_quantity_h5(file, s, "E"))

	return Es	

def get_ns(file, states):
	"""
	Get sorted ns from the .h5 file.
	"""

	return [int(i) for i in states[:,0]]

def get_Szs(file, states):
	"""
	Get sorted ns from the .h5 file.
	"""
	return [float(i) for i in states[:,1]]

def get_exciteds(file, states):
	"""
	Get sorted ns from the .h5 file.
	"""
	return [float(i) for i in states[:,2]]

def get_imp_occupancies(file, states):
	"""
	Get impurity occupancies from all states.
	"""

	impindex = get_impindex(file)
		
	occs=[]
	for s in states:
		occupancies = get_quantity_h5(file, s, "site_occupancies")
		imp_occupancy = occupancies[impindex-1]

		occs.append(imp_occupancy)

	return occs	

def get_all_occupancies(file, states):
	"""
	Get all occupancies from all states.
	"""
	impindex = get_impindex(file)
		
	occs=[]
	for state in states:		
		occs.append(get_quantity_h5(file, state, "site_occupancies"))

	return occs	

def get_spectral_weights(file, states):
	"""
	Gets all spectral weights. States is passed just because it is passed everywhere else. 
	"""

	path = "/weights/{0}/{1}/{2}/"

	w0pup = get_single_spectral_weight(file, 0, 1, "up")
	w0pdn = get_single_spectral_weight(file, 0, 1, "dn")
	w0mup = get_single_spectral_weight(file, 0, -1, "up")
	w0mdn = get_single_spectral_weight(file, 0, -1, "dn")

	w1pup = get_single_spectral_weight(file, 1, 1, "up")
	w1pdn = get_single_spectral_weight(file, 1, 1, "dn")
	w1mup = get_single_spectral_weight(file, 1, -1, "up")
	w1mdn = get_single_spectral_weight(file, 1, -1, "dn")

	return [w0up, w0dn, w0mup, w0mdn, w1pup, w1pdn, w1mup, w1mdn]

def get_single_spectral_weight(file, excited, pm, updn):
	"""
	Reads a single spectral weight. If it is not found, returns zero, to avoid problems with NaNs etc.
	"""

	path = "/weights/{0}/{1}/{2}/".format(excited, pm, updn)

	try:
		hf = h5py.File(file, "r")
		w = np.array(hf.get(path))[()]
		
		if math.isnan(w):
			return 0
		else:
			return w

	except:
		return 0

def get_correlation(file, state, which_corr):
	"""
	General function to recover correlations. The impurity-impurity correlation is saved as the first element no matter what.
	"""
	impindex = get_impindex(file)

	if which_corr == "spin_correlation":
		mp = get_quantity_h5(file, state, which_corr+"/mp")
		mpImp = get_quantity_h5(file, state, which_corr+"_imp/mp")

		pm = get_quantity_h5(file, state, which_corr+"/pm")
		pmImp = get_quantity_h5(file, state, which_corr+"_imp/pm")

		zz = get_quantity_h5(file, state, which_corr+"/zz")
		zzImp = get_quantity_h5(file, state, which_corr+"_imp/zz")

		#mp = np.array([[mpImp]+mp])
		#pm = np.array([pmImp]+[pm])
		#zz = np.array([zzImp]+[zz])
		mp = np.insert(mp, 0, mpImp)
		pm = np.insert(pm, 0, pmImp)
		zz = np.insert(zz, 0, zzImp)

		corrs = zz + 0.5*(mp + pm)

	else:	
		corrs = get_quantity_h5(file, state, which_corr)
		corrs = put_impindex_to_front(corrs, impindex)
	
	return corrs

def get_charge_correlations(file, states):
	"""
	Gets charge correlations. 
	"""
	corrs = get_correlation(file, states, "charge_correlation")
	return corrs

def get_pair_correlations(file, states):
	"""
	Gets pair correlations. 
	"""
	corrs = get_correlation(file, states, "pair_correlation")
	return corrs

def get_spin_correlations(file, states):
	"""
	Gets spin correlations. 
	"""
	corrs=[]
	for state in states:
		corrs.append(get_correlation(file, state, "spin_correlation"))
	
	return corrs

def get_hopping_correlations(file, states):
	"""
	Gets hooping correlations. 
	"""
	corrs = get_correlation(file, states, "hopping_correlation")
	return corrs

###################################################################################################
###################################################################################################
##           THESE ARE OLD FUNCTIONS, USING THE TEXT OUTPUT TO PARSE THE QUANTITIES              ##	
###################################################################################################
###################################################################################################
#FUNCTIONS THAT GET THE PATH TO output AS THE INPUT AND RETURN A LIST OF OBSERVABLES

def print_times_text(file):
	with open(file, "r") as resF:
		for line in resF:
			lastline = line		
	
	a = re.search("Wall time: (\d+)", lastline)
	if a:
		time = int(a.group(1))
	return [time]

def get_energies_text(file, returnStates=False):
	"""
	Opens the DMRG output and recovers the energies.
	"""

	ns, iis, Szs, Es = [], [], [], []	
	with open(file, "r") as resF:
		for line in resF:
			a = re.search("n = ([0-9]+).*Sz = (-?\d*.?\d*)\s+i = (\d+)\s+E\s=\s(-?\d*.?\d*).*DeltaE.*", line) #this is specifically for energy in the standard DMRG output
			b = re.search("n = ([0-9]+).*Sz = (-?\d*.?\d*).*?E\s=\s(-?\d*.?\d*)", line) #this is specifically for energy in the standard DMRG output
			
			if a:
				n = int(a.group(1))
				Sz = float(a.group(2))
				i = int(a.group(3))
				E = float(a.group(4))

				ns.append(n)
				iis.append(i)
				Szs.append(Sz)
				Es.append(E)

			elif b:
				n = int(b.group(1))
				Sz = float(b.group(2))
				E = float(b.group(3))
				
				ns.append(n)
				Szs.append(Sz)
				Es.append(E)		
		
			#THIS IS IF THE OUTPUT FILE IS OLD AND DOES NOT HAVE Szs
			else:
				c = re.search("n = ([0-9]+).*E.*= (-?\d*.?\d*)", line) #this is specifically for energy in the standard DMRG output
				if b:
					n = int(c.group(1))
					E = float(c.group(2))
					
					ns.append(n)
					Es.append(E)					

	if len(Es) == 0:
		print("NO ENERGIES FOUND IN " + str(file))
		return []

	else:
		if len(iis)==0:
			if len(Szs)==0:
				sortedEs = [E for _, E in sorted(zip(ns, Es), key = lambda x : x[0] )]
				allstates =  sorted(zip(ns, Es), key = lambda x : x[0] )
			else:
				sortedEs = [E for _, _, E in sorted(zip(ns, Szs, Es), key = lambda x : (x[0], x[1]) )]
				allstates =  sorted(zip(ns, Szs, Es), key = lambda x : (x[0], x[1]) )
		else:	
			sortedEs = [E for _, _, _, E in sorted(zip(ns, Szs, iis, Es), key = lambda x : (x[0], x[1], x[2]) )]
			allstates =  [[x, y, z, w] for x, y, z, w in sorted(zip(ns, Szs, iis, Es), key = lambda x : (x[0], x[1], x[2]) ) ]

	if returnStates:
		return allstates
	else:
		return sortedEs

def get_ns_text(file):
	"""
	Opens the DMRG output and recovers the energies.
	"""

	ns=[]	
	with open(file, "r") as resF:
		for line in resF:			
			a = re.search("n = ([0-9]+).*E = .*", line) 
			if a:
				n = int(a.group(1))
				ns.append(n)			
				
	if len(ns) == 0:
		print("NO ns FOUND IN " + str(file))
		return []

	else:
		ns = sorted(ns)	#sort by ns
					
	return ns

def get_Szs_text(file):
	"""
	Opens the DMRG output and recovers the energies.
	"""

	ns, Szs = [], []	
	with open(file, "r") as resF:
		for line in resF:			
			a = re.search("n = ([0-9]+).*Sz = (-?\d*.?\d*).*E.*= (-?\d*.?\d*)", line) #this is specifically for energy in the standard DMRG output
			if a:
				n = int(a.group(1))
				Sz = float(a.group(2))

				ns.append(n)			
				Szs.append(Sz)	

	if len(Szs) == 0:
		print("NO Szs FOUND IN " + str(file))
		return []

	else:
		sortedSzs = [Sz for _, Sz in sorted(zip(ns, Szs), key = lambda x : (x[0], x[1]) )]
					
	return sortedSzs

def get_exciteds_text(file):

	ns, Szs, iis = [], [], []	
	with open(file, "r") as resF:
		for line in resF:
			a = re.search("n = ([0-9]+).*Sz = (-?\d*.?\d*).*i = (\d+).*E.*= (-?\d*.?\d*)", line) #this is specifically for energy in the standard DMRG output
			if a:
				n = int(a.group(1))
				Sz = float(a.group(2))
				i = int(a.group(3))

				ns.append(n)			
				Szs.append(Sz)	
				iis.append(i)
	
	if len(iis) == 0:
		print("NO is FOUND IN " + str(file))
		return []

	else:
		sortediis = [i for _, _, i in sorted(zip(ns, Szs, iis), key = lambda x : (x[0], x[1], x[2]) )]
					
	return sortediis

def get_spectral_weights_text(file):

	with open(file, "r") as resF:
		#save the occupancy of the ground state sector
		wpup, wpdn, wmup, wmdn = 0, 0, 0, 0
		for line in resF:
			b = re.search("weight w\+ up: (\-*\d+\.*\d*)", line)
			c = re.search("weight w\+ dn: (\-*\d+\.*\d*)", line)
			d = re.search("weight w\- up: (\-*\d+\.*\d*)", line)
			e = re.search("weight w\- dn: (\-*\d+\.*\d*)", line)

			if b:
				wpup = float(b.group(1))
			if c:
				wpdn = float(c.group(1))
			if d:
				wmup = float(d.group(1))
			if e:
				wmdn = float(e.group(1))

	return [wpup, wpdn, wmup, wmdn]			

def get_imp_occ_nSz_text(file, nchange, Szchange):
	"""
	Opens the DMRG output and recovers the impurity occupation of the sector with n = nGS +  nchange and Sz = SzGS + Szchange.
	"""

	#find the impurity index
	with open(file, "r") as resF:
		for line in resF:
			a = re.search("impindex=(\d+)", line)
			if a:
				impindex = int(a.group(1))
				break

	with open(file, "r") as resF:
		
		#find the ground state sector
		ns, Szs, Es = [], [], []
		for line in resF:
			a = re.search("n = ([0-9]+).*Sz = (-?\d*.?\d*)\s+i = (\d+)\s+E\s=\s(-?\d*.?\d*).*DeltaE.*", line) #this is specifically for energy in the standard DMRG output
			b = re.search("n = ([0-9]+).*Sz = (-?\d*.?\d*).*?E\s=\s(-?\d*.?\d*)", line) #this is specifically for energy in the standard DMRG output

			if a:
				n = int(a.group(1))
				Sz = float(a.group(2))
				i = int(a.group(3))
				E = float(a.group(4))

				ns.append(n)
				iis.append(i)
				Szs.append(Szs)
				Es.append(E)

			elif b:
				n = int(b.group(1))
				Sz = float(b.group(2))
				E = float(b.group(3))
				
				ns.append(n)
				Szs.append(Szs)
				Es.append(E)	

	if len(Es) == 0:
		print("NO ENERGIES FOUND IN " + str(file))
		return []
	else:
		if len(Szs)==0:
			sortedEsNs = [[E, n] for E, n in sorted(zip(Es, ns), key = lambda x : x[0] )]
			nGS = sortedEsNs[0][1]
			SzGS = None

		else:
			sortedEsNsSzs = [[E, Sz, n] for E, Sz, n in sorted(zip(Es, Szs, ns), key = lambda x : (x[0], x[1]) )]
			nGS = sortedEsNsSzs[0][2]
			SzGS = sortedEsNsSzs[0][1]

	saved=False

	#save the impurity occupancy for GS
	with open(file, "r") as resF:
		#save the impurity occupancy of the GS sector
		sector=False
		for line in resF:
			#if SzGS:
			#	a = re.search("RESULTS FOR THE SECTOR WITH {0} PARTICLES, Sz {1}".format(int(nGS+nchange), SzGS+Szchange), line)
			#else:
			#	a = re.search("RESULTS FOR THE SECTOR WITH {0} PARTICLES".format(int(nGS+nchange)), line)
			
			a = re.search("RESULTS FOR THE SECTOR WITH {0} PARTICLES".format(int(nGS+nchange)), line)
			if a:
				sector=True

			if sector:
				b = re.search("site occupancies", line)
				if b:
					occupancies = re.findall("\d*\.\d+|\d+", line)

					occupancies = [float(i) for i in occupancies]

					impOcc = occupancies[impindex-1]
					sector=False
					saved=True
					break

	if not saved:
		print("THE SECTOR WITH nGS+{}={} PARTICLES WAS NOT FOUND! File: {}".format(nchange, nGS+nchange, file))				
		return []
	else:		
		return [impOcc]				

def get_entropies_text(file):
	"""
	Opens the DMRG output and recovers the entropies.
	"""

	entropies, ns = [], []
	with open(file, "r") as resF:
		for line in resF:
			
			a = re.search("RESULTS FOR THE SECTOR WITH (\d+) PARTICLES:", line)
			b = re.search("Entanglement entropy across impurity bond b=\d, SvN = (-?\d+\.*\d*)", line)

			if a: 
				n=int(a.group(1))
				ns.append(n)

			#READ OUT THE ENTROPY
			if b:
				SvN = float(b.group(1))
				entropies.append(SvN)

	if len(entropies)==0:
		return []		

	else:	
		sortedEntropies = [S for _, S in sorted(zip(ns, entropies))]			
		return sortedEntropies

def get_imp_occupancies_text(file):
	"""
	Gets impurity occupation from the text output file.
	"""
	impindex = get_impindex(file)

	states = get_energies_text(file, returnStates=True)

	impoccs = []
	with open(file, "r") as resF:
		for line in resF:
			a = re.search("RESULTS FOR THE SECTOR WITH (\d+) PARTICLES, Sz (\d+.?\d?), state (\d+):", line)
			b = re.search("site occupancies = ", line)
			
			if a:
				n = int(a.group(1))
				Sz = float(a.group(2))
				ii = int(a.group(3))

			if b:
				occupancies = re.findall("\d*\.\d+|\d+", line)

				occupancies = [float(i) for i in occupancies]
				occ = occupancies[:impindex-1] + occupancies[impindex:]
				impocc = occupancies[impindex-1]
			
				impoccs.append([n, Sz, ii, impocc])


	sortedimpoccs = [n for _, _, _, n in sorted(impoccs, key = lambda x : (x[0], x[1], x[2]))]			
		
	return sortedimpoccs
	
###################################################################################################
#ADDITIONAL AUXILIARY FUNCTIONS

def get_occs_text(n, result_file, verbose=True):
	"""
	Gets level occupancies, impurity occupancy and energy for a sector with n particles.
	"""
	if verbose:
		print("n is ", n)
	
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

def get_nupndn_text(n, result_file):
	"""
	Reads the line with info about nup and ndn occupancies.
	"""	
	with open(result_file,  "r") as resF:
		sector=False
		for line in resF:
			a = re.search("RESULTS FOR THE SECTOR WITH {0} PARTICLES:".format(n), line)
			b = re.search("impurity nup ndn =", line)
			

			if a:
				sector=True

			if sector and b:
				nupndn = re.findall("\d*\.\d+|\d+", line)
				nup, ndn = nupndn[0], nupndn[1]

				sector=False
				break

	return nup, ndn


