import networkx as nx
import random
import pickle,os
import math
from sklearn.metrics import *
from mixmod_wt.Status import Status
from mixmod_wt.mixmod_wt_correctingImplementation_without_half import __modularity
#from mixmod_wt.modularity_LFR_our import getModularityQ as __modularity
from collections import defaultdict
import matplotlib.pyplot as plt
from mixmod_wt.aux import _get_com_wise_nodes, printsomeinformation, read_raw_network
from multiprocessing import Pool

def read(filename):
	dtcom = {}
	count=1
	cmufile =open(filename)
	for line in cmufile:
		line = set(int(each)+1 for each in line.split())
		dtcom[count] = line
		count+=1
	return dtcom

def computegtmod(filename,weighted):
	fnetwork = 0
	#with open(filename+'_ml_network.pickle') as handle:
	#	fnetwork = pickle.load(handle)
	ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = read_raw_network(filename,weighted)
	
	status = Status()
	status.layer=layer
	status.node_l=node_l
	status.node_c=node_c
	status.top=top
	status.bot=bot
	status.edge_l=edge_l
	status.edge_c=edge_c
	status.couple = couple
	status.mu = mu
	mod = __modularity(commu, status, ml_network)
	return mod,commu

def compute_nmi(gtcom,dtcom):
  num_nodes = 200
  true_labels = [None]*num_nodes
  pred_labels = [None]*num_nodes

  linenum = 1
  for commu in gtcom:
    for c in list(gtcom[commu]):
        true_labels[c-1] = linenum
    linenum+=1

  linenum =1
  for d in dtcom:
    for node in list(dtcom[d]):
    	pred_labels[node-1] = linenum
    linenum+=1


  #Normalised mutual information. Function present in sklearn.metrics. Do not forget to import.
  nmi = normalized_mutual_info_score(true_labels, pred_labels) 
  return nmi

def _calculate_dtmod(filename,dtcom,weighted):
	ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = read_raw_network(filename,weighted)
	
	status = Status()
	status.layer=layer
	status.node_l=node_l
	status.node_c=node_c
	status.top=top
	status.bot=bot
	status.edge_l=edge_l
	status.edge_c=edge_c
	status.couple = couple
	status.mu = mu
	mod = __modularity(dtcom, status, ml_network)
	return mod

def calculatemod(networkpath,cmuFolder, modfilename,weighted):
	cmulist = os.listdir(cmuFolder)
	for eachfile in cmulist:
		dtcom = read(cmuFolder+eachfile)
		netfile = eachfile.split('_merge')[0]
		print(netfile)
		requirednet = networkpath+netfile
		gtmod,gtcom = computegtmod(requirednet,weighted)
		dtmod = _calculate_dtmod(requirednet,dtcom,weighted)
		nmi = compute_nmi(gtcom,dtcom)
		modfile = open(modfilename, 'a')
		modfile.write(netfile+ ":	"+ str(gtmod)+"  "+str(dtmod)+"  "+str(nmi)+"\n")
		modfile.close()

def main():
	weighted =0
	cmuFolder = "./CompMod/temphetero/"
	networkpath = "./syntheticNetworkGeneration/netsForcomparingBaseline/_networks/"
	pathtosave = "./syntheticNetworkGeneration/netsForcomparingBaseline/results/"
	modfilename = "compmod"
	f=open(pathtosave+modfilename,'w')
	f.write("Network 			GroundTruth    Detected   NMI\n")
	f.close()
	calculatemod(networkpath, cmuFolder,pathtosave+modfilename,weighted)

if __name__ == '__main__':
	main()