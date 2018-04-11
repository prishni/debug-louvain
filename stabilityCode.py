import networkx as nx
import random,math
import pickle,os,sys
from sklearn.metrics import *
from mixmod_wt.Status import Status
from mixmod_wt.mixmod_wt_correctingImplementation_without_half import __modularity
from collections import defaultdict
from mixmod_wt.aux import _get_com_wise_nodes, printsomeinformation, read_raw_network
from multiprocessing import Pool
from pprint import pprint
from copy import deepcopy
from union import getSeries as unionGetSeries
from union import computegtmod,compute_nmi,partition_at_level,__one_level,__renumber
from union import induced_graph_multilayer,_get_commu_dict
stabilitymodif = True

def compute_nmi(gtcom,dtcom):
  num_nodes = 200
  true_labels = [None]*num_nodes
  pred_labels = [None]*num_nodes

  linenum = 1
  for commu in gtcom:
    for c in list(gtcom[commu]):
        true_labels[c-1] = linenum
    linenum+=1
  
  for d in dtcom:
    for node in list(dtcom[d]):
        pred_labels[node-1] = d[node]

  #Normalised mutual information. Function present in sklearn.metrics. Do not forget to import.
  nmi = normalized_mutual_info_score(true_labels, pred_labels) 
  return nmi

def _stability(initialpart, networkname,graph, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu,weighted):
	checkstabilitymodif = False

	current_graph = graph.copy()
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

	status.init(current_graph,initialpart)
	old_status = status.copy()
	
	status_list = list()
	level_count = 0

	mod = __modularity(_get_commu_dict(status.node2com), status, graph)
	#print "Modularity before First Iteration ",mod

	__one_level(graph, old_status, current_graph, status, status_list, level_count)
	partition = __renumber(status.node2com)
	status_list.append(partition)
	current_graph,part,status = induced_graph_multilayer(partition, current_graph,status)

	mod1 = __modularity(_get_com_wise_nodes(part), status, current_graph)

	p = _get_com_wise_nodes(partition_at_level(status_list, len(status_list)-1))
	new_mod = __modularity(p, old_status, graph)
	
	A = nx.adjacency_matrix(current_graph)	
	status.init(current_graph,initialpart)
	
	while True :
		level_count+=1
		modif=__one_level(graph, old_status, current_graph, status, status_list, level_count, 1)
		
		partition = __renumber(status.node2com)
		status_list.append(partition)

		new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), old_status, graph)  
		
		if modif==False:
			checkstabilitymodif = False		
			break
		checkstabilitymodif = True
		mod = new_mod
		current_graph, part,status = induced_graph_multilayer(partition, current_graph,status)
		status.init(current_graph, part)
		
	return status_list, mod ,checkstabilitymodif

def stability(part, networkname,maxmod,weighted):
	bestpart = deepcopy(part)
	global stabilitymodif
	itr =0
	while stabilitymodif == True and itr<2:
		graph, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = read_raw_network(networkname,weighted)
		
		initial_part = deepcopy(part)
		dendogram, mod ,stabilitymodif = _stability(part, networkname,graph, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu,weighted)
		
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
		mod_old = mod
		part = partition_at_level(dendogram, len(dendogram)-1)
		commu =_get_com_wise_nodes(part)
		mod = __modularity(commu, status, graph)
		if(mod>maxmod):
			maxmod = mod
			bestpart = deepcopy(part)
		itr+=1
		if(compute_nmi(initial_part, part)):
			stabilitymodif = False
	return maxmod,bestpart

def runformanynetworks(args):
	output = []
	networklist  = args[:-1]
	tmp = args[-1]
	weighted = args[-1][0]
	modfilename = args[-1][1]
	networkpath = args[-1][2]
	p = args[-1][3]
	for network in networklist:
		str2 = networkpath + str(network)
		dtmod, dtcom = unionGetSeries(str2,weighted)
		gtmod,gtcom = computegtmod(str2,weighted)
		output.append((gtmod, dtmod))
		
		dtcom1 = _get_com_wise_nodes(partition_at_level(dtcom, len(dtcom)-1))
		nmi = compute_nmi(gtcom,dtcom1)
		#dtcom = _get_com_wise_nodes(partition_at_level(dtcom, len(dtcom)-1))
		part = partition_at_level(dtcom, len(dtcom)-1)
		maxmod,bestpart = stability(part,str2,dtmod,weighted)
		if(maxmod-dtmod>0):
			Print("***************")
		modfile = open(modfilename, 'a')
		modfile.write(str2+ ":	"+ str(gtmod)+"  "+str(maxmod)+"  "+str(nmi)+"\n")
		modfile.close()

		#WRITE DETECTED COMMUNITIES TO FILE
		if(not os.path.exists(p)):
			os.makedirs(p)
		#p = "./syntheticNetworkGeneration/diggnetwork/"
		dtcom = _get_com_wise_nodes(partition_at_level(dtcom, len(dtcom)-1))

		comfile = open(p+ network,'w')
		comfile.write(str2+":\n============================================================\n")
		comfile.write(str(len(dtcom))+"\n")
		for d in dtcom:
			towrite = ' '.join([str(node) for node in dtcom[d]])
			comfile.write(towrite + "\n")
		comfile.write("=======================================================\n")
		comfile.close()

	return output

def parallelimplementation(networklist,weighted,modfilename,networkpath,p):
	#Open file to compare modularity values
	modfile = open(modfilename,'w')
	modfile.write("network								   GroundTruth	Detected-Louvain\n")
	modfile.close()
	
	#Prepare args for parallel processings
	numnetworks = len(networklist)
	cores=24
	chunksize = numnetworks/cores
	splits = []
	for i in range(cores):
		splits.append((i)*chunksize)
	splits.append(numnetworks)

	args = []
	for i in range(cores):
		arguments = networklist[splits[i]:splits[i+1]]
		tmp = [weighted,modfilename,networkpath,p]
		arguments.append(tmp)
		args.append(arguments)
	p = Pool(cores)
	modularities = p.map(runformanynetworks, args)
	#Flatten list of lists returned by different cores
	modularities = [item for items in modularities for item in items]
	return modularities

def main():
	weighted=0
	
	networkpath="./syntheticNetworkGeneration/netsForcomparingBaseline/_networks/"
	#networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/morenets/alpha0.7/"
	#networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/netsByGenerateNetsv3/alpha0.7/"

	pathtosave="./syntheticNetworkGeneration/netsForcomparingBaseline/results/"
	#pathtosave="./syntheticNetworkGeneration/results/nets121/alpha0.7/dt_db_mod_files_corrctImp_without_half/"
	#pathtosave="./syntheticNetworkGeneration/results/nets121_incWeights/alpha0.7/dt_db_mod_files_corrctImp_without_half/"

	networklist = os.listdir(networkpath)
	networklist=sorted(networklist)
	modfilename = pathtosave+"union_WithStablility"
	p= pathtosave + "DetectedCommUnionStability/"
	parallelimplementation(networklist,weighted,modfilename,networkpath,p)
	
if __name__ == '__main__':
	main()
