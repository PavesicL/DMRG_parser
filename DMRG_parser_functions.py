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

class parameter:
 
	def __init__(self, name, sweeptype):
		self.name = name
		self.sweeptype = sweeptype

		self.value = None

###################################################################################################

def readNameFileParsing(file, observable, noParam=False):
	"""
	Reads the nameFile and returns the path to which to save the files and the parameter over which to sweep.
	If noParam, there is no parameter to sweep over, so no sweep param is defined (for example get_all_occupancies, saves occs
	for one calculation only). 
	"""

	obsCheck = False
	param = None

	with open(file, "r") as f:
		for line in f:
			line = line.strip()

			if len(line)==0:
				continue

			if line[0]=="#":	#comment
				continue

			a = re.search(observable+"{", line)
			b = re.search("}"+observable, line)

			c = re.search("path\s*=\s*(.*)", line)
			d = re.search("\s*sweep\s*(.*)", line)
			
			if a:
				obsCheck=True
			if b:
				obsCheck=False
			
			if c and obsCheck:
				path = c.group(1)
			if d and obsCheck:
				param = d.group(1)

	try: 		
		return param, path.strip()			
	
	except UnboundLocalError:
		return 0, 0


def getParamsNameFile(file):
	"""
	Reads the nameFile, and returns a generic jobname and a dictionary of parameters, 
	where key is the parameter name and value is its sweeptype (sweep, case or relation).
	"""

	paramList = []

	paramsCheck = False
	with open(file, "r") as f:
		for line in f:
			line = line.strip()	#strip the leading and trailing whitespace
			if len(line)==0:
				continue

			a = re.search("name\s*=\s*(.*)", line)
		
			b = re.search("params\s*{", line) 
			c = re.search("}\s*endparams", line) 
			
			d = re.search("(\w*)\s*(\w*)", line)
			
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
				param, sweeptype = d.group(1), d.group(2)
				paramList.append(parameter(param, sweeptype))

		return name, paramList
	
	return None, None

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
	name, paramList = getParamsNameFile("nameFile")	#each key to paramsDict is a name of the parameter, with the value being its sweeptype
	regexName = re.sub("{[0-9]+}", "([+-]?[0-9]+(?:\.?[0-9]*(?:[eE][+-]?[0-9]+)?)?)", name)	#replace all instances of {number} in the name with ([0-9]+.*[0-9]*), which matches floats and ints
	sweepParam, savepath = readNameFileParsing("nameFile", which_var)

	outputName = "output"	#subject to change - maybe allow user input here

	if name == None and paramList == None:
		print("Name and paramDict were not found. Check if nameFile has correct data!")
		exit()


	whichParam = determineWhichParam(sweepParam, [p.name for p in paramList])		

	doneJobs = findFinishedJobs(os.getcwd())	#a list of all jobs found in finishedJobs.txt and savedJobs.txt


	#get results from all folders and save them to EList
	EList=[]
	saved=0
	noth5, h5 = 0, 0

	result_dir = os.getcwd() + "/results/"

	for job in doneJobs:
		folder = result_dir + job
		aa = re.match(regexName, job)

		updateParamValues(paramList, aa)

		result_file_h5 = folder+"/solution.h5"
		result_file = folder+"/"+outputName

		#CHECK IF THE OUTPUT FILES EXIST
		if not os.path.isfile(result_file):
			result_file += ".txt"
			#result_file = result_file[:-4]
		if not os.path.isfile(result_file):
			print("The output file is not output or output.txt; or does not exist! Job: "+job)
			continue

		Es, A, B = getDataHDF5(job, result_file_h5, result_file, which_func, params)

		h5+=A
		noth5+=B

		if len(Es)!=0:
			EList.append([p.value for p in paramList]+[Es])	#this is now a list of values of all parameters (in the same order as in paramList), and corresponding Es
			saved+=1
	
	EList = sortByParamValues(EList, whichParam, paramList)

	saved_sets = sliceAndSave(EList, paramList, whichParam, savepath)
	print("{} were saved from the hdf5 file, while {} were NOT.".format(h5, noth5))

	return saved_sets, saved

###################################################################################################
#SUBFUNCTIONS, USED IN THE MAIN FUNTIONS ABOVE

def determineWhichParam(sweepParam, listOfParams):
	"""
	Determines which parameter is the sweep parameter.
	"""

	whichParam=None
	for i in range(len(listOfParams)):
		param = listOfParams[i]
		if param == sweepParam:
			whichParam = i 
			break

	return whichParam

def findFinishedJobs(folder):
	"""
	Returns a list of jobs that are in finishedJobs.txt or savedJobs.txt.
	"""
	finished = readFileToList(folder + "/finishedJobs.txt")
	saved = readFileToList(folder + "/savedJobs.txt")

	return finished + saved

def updateParamValues(paramList, regex):
	i=0
	for p in paramList:
		i+=1
		p.value = float(regex.group(i))

	return None

def getDataHDF5(job, result_file_h5, result_file, which_func, params=[]):
	"""
	Calls the desired parsing function and saves the recovered data.
	"""
	
	h5, noth5 = 0, 0

	#TRY CALLING THE HDF5 PARSING FUNCTION. IF IT FAILS (useful if there is no .h5 file, typically in old jobs), FALLBACK TO TEXT PARSING. 
	try:
		"""
		This magic calls the function with the name which_func, with the argument result_file and *params
		This is here in order to be able to pass a name of the function as an argument, and call it here (eg. "get_energies" is passed as a param, so in this line the
		function get_energies() is called). 
		"""

		states = get_all_states(result_file_h5)							
		this_module = sys.modules[__name__]

		if params==[]:
			Es = getattr(this_module, which_func)(result_file_h5, states)
		else:
			Es = getattr(this_module, which_func)(result_file_h5, states, params)

		h5+=1

	except:
		this_module = sys.modules[__name__]
		try:
			Es = getattr(this_module, which_func+"_text")(result_file)
			noth5+=1 	#to check how many files were saved not from the h5 files

		except:
			Es=[]
			print("Not found! file: "+job)
	
	return Es, h5, noth5

def sortByParamValues(EList, whichParam, paramList):
	"""
	Sort EList by all non sweep parameters first, and then by the sweep parameter.
	Allows EList to be sliced where the parameters except the sweep one change.	
	"""
	if whichParam == None:
		EList = sorted(EList, key = lambda x : ([x[i] for i in range(len(EList[0])-1)]) )
	else:
		EList = sorted(EList, key = lambda x : ([x[i] for i in range(whichParam)] + [x[i] for i in range(whichParam+1,len(paramList))] + [x[whichParam]]) )

	return EList	

def getParamsToDisregard(paramList, savepath):
	"""
	The parameters that have to be taken into account here are only the ones that are mentioned in the savepath; and the sweep one.
	Others can be discounted. For example the relation parameters; gamma2 = gamma1; sweep is along gamma1, no need to have gamma2 anywhere.
	This function gets the indeces of these parameters.
	"""
	indexList = []
	i=-1
	for param in paramList:
		i+=1
		name = param.name
		a = re.search("_"+name+"{}", savepath)
		if a == None:
			indexList.append(i)

	return indexList

def sliceAndSave(EList, paramList, whichParam, savepath):
	"""
	Reshapes the EList so that the values can be saved accoring to the sweep parameter and saves them.
	"""

	disregardIndeces = getParamsToDisregard(paramList, savepath)	#disregards parameters that are not mentioned in the savepath, and the sweep parameter.

	tempList = []
	saved_sets = 0
	for i in range(len(EList)-1):

		currentParamValues = takeOutParams(EList[i], disregardIndeces)
		nextParamValues = takeOutParams(EList[i+1], disregardIndeces)

		if whichParam==None:
			tempList.append(EList[i][-1])
		else: 
			tempList.append([EList[i][whichParam]]+list(EList[i][-1]))

		#If the next parameter values are different to current, save and reset the tempList  
		if newParamBatch(nextParamValues[:-1], currentParamValues[:-1]):# and i>0:
			saveToFile(tempList, savepath, currentParamValues)
			saved_sets +=1
			tempList=[]

	#append and save the last case
	if whichParam==None:
		tempList.append(EList[-1][-1])
	else: 	
		tempList.append([EList[-1][whichParam]]+list(EList[-1][-1]))
	
	saveToFile(tempList, savepath, nextParamValues)
	saved_sets += 1

	return saved_sets

###################################################################################################
#UTILITY FUNCTIONS

def readFileToList(file):
	ll=[]

	try:
		with open(file, "r") as f:
			for line in f:
				ll.append(line.strip())
	except FileNotFoundError:
		pass

	return ll

def valuesWithoutSweep(i, paramList, whichParam):
	if whichParam == None:
		return np.array([EList[i][j] for j in range(len(paramList))])

	else:
		return np.array([EList[i][j] for j in range(whichParam)] + [EList[i][j] for j in range(whichParam+1,len(paramList))])

def newParamBatch(oldSet, newSet):
	for i in range(len(oldSet)):
		if oldSet[i] != newSet[i]:
			return True
	return False

def takeOutParams(valueList, takeOutList):
	"""
	Takes out elements from the valueList at indexes given in takeOutList. Works also if takeOutList == [].
	"""
	for i in sorted(takeOutList, reverse=True):
		newList = np.array([valueList[j] for j in range(i)] + [valueList[j] for j in range(i+1,len(valueList))])
		valueList = newList
	
	return valueList	

def saveToFile(tempList, savepath, parameterValues):
	with open(savepath.format(*parameterValues), "w") as txt_file:
		for E in tempList:
			txt_file.write("	".join([str(i) for i in E]) + "\n") # works with any number of elements in a line
			
def fix_string_Sz(Sz):

	Sz = float(Sz)
	if Sz == 0:
		Sz = "0"
	if Sz == 0.5:
		Sz = "0.5"
	if Sz == -0.5:
		Sz = "-0.5"
	if Sz == 1:
		Sz = "1"
	if Sz == -1:
		Sz = "-1"

	return Sz		

###################################################################################################
#UTILITY PARSING FUNCTIONS

def read_p_param(param, file):
	"""
	Reads a parameter from the parameter class - these are the ones printed out at the start of output.txt.
	"""

	with open(file, "r") as f:
		for line in f:
			a = re.search("params." + param + " = (.*)", line)
			if a:
				res = a.group(1)
				return res

	print(f"Param {param} not found in file {file}!")
	return None

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

def get_total_spin(file, states):
	"""
	Get total spin from all states.
	"""
	
	ll=[]
	for s in states:
		SS = get_quantity_h5(file, s, "S2")

		ll.append(SS)

	return ll

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

	elif which_corr == "hopping_correlation":
		corrsUp = get_quantity_h5(file, state, "hopping/up")
		corrsDn = get_quantity_h5(file, state, "hopping/dn")
		
		#THERE IS NO IMP TERM HERE!
		
		#corrs = np.abs(corrsUp) + np.abs(corrsDn)
		corrs = corrsUp + corrsDn

	else:
		corrs = get_quantity_h5(file, state, which_corr)
				#NO IMPURITY TERM HERE ALSO!
	
	return corrs

def get_charge_correlations(file, states):
	"""
	Gets charge correlations. 
	"""
	corrs=[]
	for state in states:
		corrs.append(get_correlation(file, state, "charge_correlation"))
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
	Gets hopping correlations. 
	"""
	corrs=[]
	for state in states:
		corrs.append(get_correlation(file, state, "hopping_correlation"))
	
	return corrs

def get_imp_amplitudes_zero(file, states):
	"""
	Gets impurity amplitudes.
	"""
	amps = []
	for s in states:
		amp = get_quantity_h5(file, s, "imp_amplitudes/0")
		amps.append(amp)
	return amps

def get_imp_amplitudes_up(file, states):
	"""
	Gets impurity amplitudes.
	"""
	amps = []
	for s in states:
		amp = get_quantity_h5(file, s, "imp_amplitudes/up")
		amps.append(amp)
	return amps

def get_imp_amplitudes_down(file, states):
	"""
	Gets impurity amplitudes.
	"""
	amps = []
	for s in states:
		amp = get_quantity_h5(file, s, "imp_amplitudes/down")
		amps.append(amp)
	return amps

def get_imp_amplitudes_two(file, states):
	"""
	Gets impurity amplitudes.
	"""
	amps = []
	for s in states:
		amp = get_quantity_h5(file, s, "imp_amplitudes/2")
		amps.append(amp)
	return amps

def get_charge_susceptibility(file, states, which):
	"""
	Gets charge susceptibilities. 
	Reads how many excited states are computed in a given result and saves all the combinations.
	For all unique combinations of (n, Sz) gets all charge susceptibility combinations, without double counting.
	"""

	i, j = which

	nSz = np.unique([[s[0], s[1]] for s in states], axis=0)	#these are all pairs of n, Sz; taking out the excited states.

	css = []
	for n, Sz in nSz:
			
		Sz = fix_string_Sz(Sz)	

		cs = get_single_charge_susceptibility(file, n, Sz, i, j)
		css.append(cs)

	return css	

def get_single_charge_susceptibility(file, n, Sz, i, j):
	"""
	Reads a singlet charge susc number. If it is not found returns zero.
	"""

	path = f"/charge_susceptibilty/{n}/{Sz}/{i}/{j}/"

	try:
		hf = h5py.File(file, "r")
		cs = np.array(hf.get(path))[()]
		
		if math.isnan(cs):
			return 0
		else:
			return cs
	except:
		return 0

def get_which_overlaps(file, states, which):
	overlaps = []
	for s in states:
		overlaps.append(get_quantity_h5(file, s, which))
	return overlaps

def get_ZBA_overlaps_BCSL(file, states):
	"""
	This is from the Lanczos code!
	"""
	o = get_which_overlaps(file, states, "doublet_overlaps_BCSL")
	return o

def get_ZBA_overlaps_BCSR(file, states):
	"""
	This is from the Lanczos code!
	"""
	o = get_which_overlaps(file, states, "doublet_overlaps_BCSR")
	return o

def get_ZBA_overlaps_OS(file, states):
	"""
	This is from the Lanczos code!
	"""
	o = get_which_overlaps(file, states, "doublet_overlaps_OS")
	return o



###################################################################################################
###################################################################################################
##           THESE ARE OLD FUNCTIONS, USING THE TEXT OUTPUT TO PARSE THE QUANTITIES              ##	
###################################################################################################
###################################################################################################
#FUNCTIONS THAT GET THE PATH TO output AS THE INPUT AND RETURN A LIST OF OBSERVABLES

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
				if c:
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

	#states = get_energies_text(file, returnStates=True)


	impoccs = []
	with open(file, "r") as resF:
		for line in resF:
			a = re.search("RESULTS FOR THE SECTOR WITH (\d+) PARTICLES, Sz (\d+.?\d?), state (\d+):", line)
			olda = re.search("RESULTS FOR THE SECTOR WITH (\d+) PARTICLES, Sz (\d+.?\d?):", line)
			oldera = re.search("RESULTS FOR THE SECTOR WITH (\d+) PARTICLES:", line)
			b = re.search("site occupancies = ", line)
				

			if a:
				n = int(a.group(1))
				Sz = float(a.group(2))
				ii = int(a.group(3))

			if olda:
				n = int(olda.group(1))
				Sz = float(olda.group(2))
				ii = 0

			if oldera:
				n = int(oldera.group(1))
				Sz = 0
				ii = 0

			if b:
				occupancies = re.findall("\d*\.\d+|\d+", line)

				occupancies = [float(i) for i in occupancies]
				occ = occupancies[:impindex-1] + occupancies[impindex:]
				impocc = occupancies[impindex-1]

				impoccs.append([n, Sz, ii, impocc])

	sortedimpoccs = [nimp for _, _, _, nimp in sorted(impoccs, key = lambda x : (x[0], x[1], x[2]))]					

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


