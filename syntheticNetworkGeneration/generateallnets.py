#***************************************************************************************************
#*
#*	The code will generate sets of 121 synthetic networks faor varying values of Alpha , P and P1.
#*	
#*	FUNCTIONS:
#*		1. callForNonWeighted : Generates networks by adding/removing edges
#*   			                and saves them in ./netsForDtDmDb/_networks/params/nonweighted/
#*		2. callForWeighted :    Generates networks by inc/dec edges weights
#*   			                and saves them in ./netsForDtDmDb/_networks/params/weighted/
#*****************************************************************************************************
import os, sys, time
import subprocess
import networkx as nx
import random
import pickle
import random
import math
from sklearn.metrics import *
from collections import defaultdict
import matplotlib.pyplot as plt
from multiprocessing import Pool
from copy import deepcopy
from generateNets2 import getSeries as getSeriesNonweighted
from generateNets3 import getSeries as getSeriesWeighted


alphalist = [0.3,0.5,0.8]
plist =[0.3,0.5,0.8]
p1list =[0.3,0.5,0.8]

def callForNonWeighted():
	Basefilepath = "./netsForDtDmDb/_networks/baseNetworks/"
	filepath = "./netsForDtDmDb/_networks/params/nonweighted/"
	count =0

	for alpha in alphalist:
		file1 = "network_" + str(alpha) + "_"
		directory1 = filepath+ "aplha" + str(alpha) 
		if not os.path.exists(directory1):
			os.makedirs(directory1)
		for p in plist:
			file2 = file1 + str(p) + "_0.05_"
			directory2 = directory1 + "/P" + str(p)
			if not os.path.exists(directory2):
				os.makedirs(directory2)
			for p1 in p1list:
				file3 = file2 + str(p1) + "_0.0_10"
				directory3 = directory2 + "/P1_" + str(p1)
				if not os.path.exists(directory3):
					os.makedirs(directory3)
				basefile = Basefilepath + file3
				newfilepath  = directory3 + "/" 			
				getSeriesNonweighted(basefile , newfilepath + file3 +"_")
				count+=1
				print("{3}. Done for alpha: {0} P: {1} P1: {2}".format(alpha,p,p1,count))


def callForWeighted():
	Basefilepath = "./netsForDtDmDb/_networks/baseNetworks/"
	filepath = "./netsForDtDmDb/_networks/params/weighted/"
	count =0

	for alpha in alphalist:
		file1 = "network_" + str(alpha) + "_"
		directory1 = filepath+ "aplha" + str(alpha) 
		if not os.path.exists(directory1):
			os.makedirs(directory1)
		for p in plist:
			file2 = file1 + str(p) + "_0.05_"
			directory2 = directory1 + "/P" + str(p)
			if not os.path.exists(directory2):
				os.makedirs(directory2)
			for p1 in p1list:
				file3 = file2 + str(p1) + "_0.0_10"
				directory3 = directory2 + "/P1_" + str(p1)
				if not os.path.exists(directory3):
					os.makedirs(directory3)
				basefile = Basefilepath + file3
				newfilepath  = directory3 + "/" 			
				getSeriesWeighted(basefile , newfilepath + file3 +"_")
				count+=1
				print("{3}. Done for alpha: {0} P: {1} P1: {2}".format(alpha,p,p1,count))


callForNonWeighted()
callForWeighted()
