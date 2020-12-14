#!/usr/bin/env python3
"""
Parses the DMRG output file and prints out the times. 
"""

from DMRG_parser_functions import save_variables_simple
import numpy as np

def toH(seconds):
	return seconds/3600

save_variables_simple("print_times", "print_times")

timeList = np.loadtxt("times")

avg, std, mi, ma = np.average(timeList), np.std(timeList), np.min(timeList), np.max(timeList)

print("Average time: {} h, std: {} h".format(toH(avg), toH(std)))
print("min: {} h, max: {} h".format(toH(mi), toH(ma)))