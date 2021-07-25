#!/usr/bin/env python3

import os
import sys
import re
import subprocess

###################################################################################################

def get_overlap(SITES1, MPS1, SITES2, MPS2):

	#run the calcOverlap process
	overlap = subprocess.run(['/home/pavesic/git_repos/tensor/calcOverlap', SITES1, MPS1, SITES2, MPS2], capture_output=True)

	#this is the returned bitstring, decoded
	res = overlap.stdout.decode('ascii')
	
	#capture the overlap number
	overlap = re.search("overlap: (.*)", res).group(1)
	overlap = float(overlap)

	return overlap


def get_param(param, jobname):
	"""
	Get gamma from the name of the job.
	"""

	param = re.search(f"{param}(\d+.?\d*)", jobname).group(1)
	param = float(param)

	return param


###################################################################################################
#take n, Sz, ii as input

if len(sys.argv) != 6:
	print(f"usage: {sys.argv[0]} overlapwith, param, n, Sz, i")
	print("overlapwith - folder name of the job to compute all overlaps with")
	print("param - which parameter to sweep and save in the output file")
	print("n, Sz, i - which sector to take the overlaps with")
	exit()

overlapwith = sys.argv[1]
param = sys.argv[2]
n = int(sys.argv[3])
Sz = float(sys.argv[4])
i = int(sys.argv[5])

###################################################################################################
#setup the base state
MPS0 = overlapwith + f"/MPS_n{n}_S{Sz}_i{i}"
SITES0 = overlapwith + f"/SITES_n{n}_S{Sz}_i{i}"


results = sorted(os.listdir("results/"))
#run the overlaps calculation on all folders
os.chdir("results")

overlaps, params = [], []
for case in results:

	paramValue = get_param(param, case)
	
	MPSi = case + f"/MPS_n{n}_S{Sz}_i{i}"
	SITESi = case + f"/SITES_n{n}_S{Sz}_i{i}"

	overlap = get_overlap(SITES0, MPS0, SITESi, MPSi)

	print(paramValue, overlap)

	overlaps.append(overlap)
	params.append(paramValue)

#save params and overlaps to file
os.chdir("..")

with open(f"overlaps_n{n}_Sz{Sz}_i{i}", "w") as file:
	for j in range(len(params)):
		file.write(f"{params[j]}	{overlaps[j]}\n")










