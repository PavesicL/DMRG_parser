#!/usr/bin/env python3
"""
Parses the DMRG output file and prints out the times. 
"""

from DMRG_parser_functions import findFinishedJobs
import numpy as np
import os
import re

def toH(seconds):
	return round(seconds/3600, 1)

def get_walltime(file):
	with open(file, "r") as resF:
		for line in resF:
			lastline = line		
	
	a = re.search("Wall time: (\d+)", lastline)
	if a:
		time = int(a.group(1))
		return time
	else:
		print(file)
		return 0

###########################################################################		

jobs = findFinishedJobs(os.getcwd())

timeList = []
for job in jobs:
	folder = os.getcwd() + "/results/" + job

	output = folder + "/output"

	try:
		time = get_walltime(output)
	except FileNotFoundError:
		time = get_walltime(output+".txt")

	if time!=0:
		timeList.append(time)

avg, std, mi, ma = np.average(timeList), np.std(timeList), np.min(timeList), np.max(timeList)

print("Read walltime from {} jobs.".format(len(timeList)))
print("Average time: {} h, std: {} h".format(toH(avg), toH(std)))
print("min: {} h, max: {} h".format(toH(mi), toH(ma)))