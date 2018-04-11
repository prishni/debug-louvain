import networkx as nx
import random
import pickle
import math
from sklearn.metrics import *
from mixmod_wt.Status import Status
from mixmod_wt.mixmod_wt_correctingImplementation_without_half import __modularity
#from mixmod_wt.modularity_LFR_our import getModularityQ as __modularity
from collections import defaultdict
import matplotlib.pyplot as plt
from mixmod_wt.aux import _get_com_wise_nodes, printsomeinformation, read_raw_network
from multiprocessing import Pool
from union import parallelimplementation,runformanynetworks
from stabilityCode import runformanynetworks as runStability

def main():
	networkpath = "./syntheticNetworkGeneration/diggnetwork/"
	pathtosave = "./syntheticNetworkGeneration/diggnetwork/"
	weighted=0
	networklist = ['DiggNewFormat']
	modfilename = pathtosave+"testingnmi"
	detectedCommuFile = pathtosave + "digg_detected"
	tmp = [weighted,modfilename,networkpath,detectedCommuFile]
	networklist.append(tmp)
	runformanynetworks(networklist)
	#modfilename = pathtosave+"mod_digg_unionStability"
	#tmp = [weighted,modfilename,networkpath]
	#networklist.append(tmp)
	#runStability(networklist)

if __name__ == '__main__':
	main()
