import networkx as nx
import random
import pickle
import math
import json
from sklearn.metrics import *
from aux.__modularity import __git_modularity
from aux.mixmod_wt_correctingImplementation import __modularity
#from aux.modularity21_optimized import __modularityoptimized
from aux.Status import Status
from aux.show_intermediate_steps import *
from aux.__inducedGraph import __git_induced_graph
from pythonlouvainmaster import best_partition
from collections import defaultdict

__PASS_MAX = -1
__MIN = 0.0000001
      
def _get_com_wise_nodes(dictionary):
    #m = max(dictionary.values())
    louvain_p = defaultdict(set)
    for l in dictionary.keys():
        louvain_p[dictionary[l]].add(l)
    return louvain_p

def is_multi_layer(e1, e2, node_c):
    if e2 in node_c and e1 in node_c[e2]:
        return True
    return False

def is_commu(e1, e2, commu):
    for c in commu:
        if e1 in commu[c] and e2 in commu[c]:
            return True
    return False

def partition_at_level(dendogram, level) :
    partition = dendogram[0].copy()
    for index in range(1, level + 1) :
        for node, community in partition.iteritems() :
            partition[node] = dendogram[index][community]
    return partition

def __renumber(dictionary) :
    count = 0
    ret = dictionary.copy()
    new_values = dict([])

    for key in dictionary.keys() :
        value = dictionary[key]
        new_value = new_values.get(value, -1)
        if new_value == -1 :
            new_values[value] = count
            new_value = count
            count = count + 1
        ret[key] = new_value

    return ret

def _get_commu_dict(node2com):

    commu={}
    count = 1
    new_values = dict([])
    for n in node2com.keys():
        v = node2com[n]
        new_value = new_values.get(v, -1)
        if new_value == -1 :
            new_values[v] = count
            new_value = count
            commu[new_value] = set()
            count = count + 1
        commu[new_value].add(n)
    return commu


def __neighcom(node, graph, status) :
    weights = []
    for neighbor in graph[node]:
        if neighbor != node :
            neighborcom = status.node2com[neighbor]
            weights.append(neighborcom)
    return weights

def __one_level(graph, status, status_list, level_count, verbose = 0) :
    modif = True
    nb_pass_done = 0
    p_temp = status.node2com
    status_list.append(p_temp)
    cur_mod = __modularity(_get_com_wise_nodes(partition_at_level(status_list, level_count)), status,graph)
    status_list.pop()
    new_mod = cur_mod

    #print "# id_node from_com to_com local_mod mod"

    while modif  and nb_pass_done != __PASS_MAX :
        cur_mod = new_mod
        modif = False
        nb_pass_done += 1
        cur_mod2 = cur_mod

        print("outside loop")
        for node in graph.nodes():
            com_node = status.node2com[node]
            neigh_communities = __neighcom(node, graph, status)
            status.node2com[node] = -1
            best_com = com_node
            best_increase = 0
            
            for com in neigh_communities:

                temp_dict = {com:_get_com_wise_nodes(status.node2com)[com]}
                base_mod_of_community = __modularity(temp_dict, status, graph)

                status.node2com[node] = com
                
                #p_temp = __renumber(status.node2com)
                status_list.append(p_temp)

                temp_dict = {com:_get_com_wise_nodes(status_list[-1])[com]}
                incr = __modularity(temp_dict, status, graph) - base_mod_of_community
                
                #incr =  __modularity(_get_com_wise_nodes(partition_at_level(status_list, level_count)), status,graph) - cur_mod2
                status_list.pop()

                if incr > best_increase :
                    best_increase = incr
                    best_com = com

                status.node2com[node] = -1
            status.node2com[node] = best_com

            if(verbose):
                pass #print("Node ", node, " moved from ", com_node, " to ", best_com, " increase: ", best_increase)
            
            #p_temp = __renumber(status.node2com)
            status_list.append(p_temp)
            cur_mod2 =  __modularity(_get_com_wise_nodes(partition_at_level(status_list, level_count)), status, graph)
            status_list.pop()
            
            if best_com != com_node :
                modif = True
            
                '''p_temp2 = __renumber(status.node2com)
                status_list.append(p_temp2)
                incr =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status,graph) - cur_mod2
                #print "m here"
                #print node, com_node, best_com, incr, best_increase, __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status), cur_mod2
                status_list.pop()'''

        p_temp = status.node2com
        status_list.append(p_temp)
        
        #Changing mod to optimized modularity definition here
        new_mod = __modularity(_get_com_wise_nodes(partition_at_level(status_list, level_count)), status,graph)
        status_list.pop()

        if new_mod - cur_mod < __MIN :
            break

def __git_one_level(graph, status, status_list, level_count) :
    modif = True
    nb_pass_done = 0
    p_temp = __renumber(status.node2com)
    status_list.append(p_temp)
    cur_mod = __git_modularity(status)
    status_list.pop()
    new_mod = cur_mod

    print "# id_node from_com to_com local_mod mod"

    while modif  and nb_pass_done != __PASS_MAX :
        cur_mod = new_mod
        modif = False
        nb_pass_done += 1
        cur_mod2 = cur_mod

        for node in graph.nodes():
            com_node = status.node2com[node]
            neigh_communities = __neighcom(node, graph, status)
            status.node2com[node] = -1
            best_com = com_node
            best_increase = 0
            
            for com in neigh_communities:
                status.node2com[node] = com
                
                p_temp = __renumber(status.node2com)
                status_list.append(p_temp)
                incr =  __git_modularity(status) - cur_mod2
                status_list.pop()

                if incr > best_increase :
                    best_increase = incr
                    best_com = com

                status.node2com[node] = -1

            status.node2com[node] = best_com
            
            p_temp = __renumber(status.node2com)
            status_list.append(p_temp)
            cur_mod2 =  __git_modularity(status)
            status_list.pop()

            if best_com != com_node :
                modif = True
                
                p_temp2 = __renumber(status.node2com)
                status_list.append(p_temp2)
                incr =  __git_modularity(status) - cur_mod2
                #print "m here"
                #print node, com_node, best_com, incr, best_increase, __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status), cur_mod2
                status_list.pop()

        p_temp = __renumber(status.node2com)
        status_list.append(p_temp)
        new_mod = __git_modularity(status)
        status_list.pop()
        if new_mod - cur_mod < __MIN :
            break


def induced_graph(partition, graph) :
    ret = nx.Graph()
    ret.add_nodes_from(partition.values())

    for node1, node2, datas in graph.edges_iter(data = True) :
        weight = datas.get("weight", 1)
        com1 = partition[node1]
        com2 = partition[node2]
        w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
        ret.add_edge(com1, com2, weight = w_prec + weight)

    return ret
def last_partition(dendrogram,communities_itr1,level):
	partition = communities_itr1
	for index in range(0, level + 1):
	    for node, community in partition.items():
	        partition[node] = dendrogram[index][community]
	return partition

def louvain(graph, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu) :
    current_graph = graph.copy()
    status = Status()
    status.init(current_graph)
    #print(layer)

    status.layer=layer
    status.node_l=node_l
    status.node_c=node_c
    status.top=top
    status.bot=bot
    status.edge_l=edge_l
    status.edge_c=edge_c
    status.couple = couple
    status.mu = mu

    mod = __modularity(_get_com_wise_nodes(status.node2com), status,current_graph)
    status_list = list()
    level_count = 0
    __one_level(current_graph, status, status_list, level_count)
    new_mod = __modularity(_get_com_wise_nodes(status.node2com), status,current_graph)
    partition = __renumber(status.node2com)
    status_list.append(partition)
    print str(mod)+" "+str(new_mod)+" OUT"

    mod = new_mod
    '''with open(str2+'_commu_benching_all_march21_louvain_step1.pickle', 'wb') as handle:
                    pickle.dump(partition_at_level(status_list, 0), handle)'''
    prev_graph = current_graph
    current_graph = induced_graph(partition, current_graph)
    status.init(current_graph)    
    printgraph(current_graph)
    #detectedcomms(status_list)
    print "######################"

    communities_itr1 = status_list[0]
    print(communities_itr1)

    current_graph1 = nx.Graph()
    weightededgelist =[(e[0],e[1],e[2]) for e in current_graph.edges_iter(data='weight')]
    current_graph1.add_weighted_edges_from(weightededgelist)
    #print(list(current_graph1.edges(data = True)))
    dendo = best_partition(current_graph1)
    communities = []
    communities.append(communities_itr1)
    communities.extend(dendo)
    finalmodularity = __modularity(_get_com_wise_nodes(partition_at_level(communities, len(communities)-1)), status,prev_graph)
    return communities, finalmodularity
    #communities = last_partition(dendo,communities_itr1,len(dendo) - 1)
    #print(set(communities.values()))
    '''final_communities= defaultdict(list)
                val_list=[]
                for key,val in communities.items():
                    if val in val_list:
                        final_communities[val].append(communities_itr1[key])
                    else:
                        final_communities[val] = communities_itr1[key]
                        val_list.append(val)
            
                return final_communities'''
    #return status_list[:-1]
    '''
    while True :
        level_count+=1
        #print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        #printgraph(current_graph)
        __git_one_level(current_graph, status, status_list, level_count)
        partition = __renumber(status.node2com)
        status_list.append(partition)
        new_mod = __git_modularity(status)
        
        #print str(mod)+" "+str(new_mod)+" IN"
        if new_mod - mod < __MIN :
            break
        mod = new_mod
        
        current_graph = induced_graph(partition, current_graph)
        #printgraph(current_graph)
        status.init(current_graph)
        print "######################"
    detectedcomms(status_list)
    return status_list[:-1], mod
    	'''


    #current_graph=current_graph1


    #return status_list[:-1], mod
    '''   
    while True:
        __git_one_level(current_graph, status, weight, resolution, randomize)
        new_mod = __git_modularity(status)
        if new_mod - mod < __MIN:
            break
        partition = __renumber(status.node2com)
        status_list.append(partition)
        mod = new_mod
        current_graph = __git_induced_graph(partition, current_graph, weight)
        status.init(current_graph, weight)
    print("final mod= ",new_mod)
    detectedcomms(status_list)
    return status_list[:-1], mod
    #   '''
def computegtmod(filename):
    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    
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
    return mod

def getSeries(filename):
    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    #fc = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)
    dendogram, mod = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)

    return mod, dendogram
    #print(fc)

def generateinfofilefortesting(str2,str21,pathtosave):   
    infofile = open(str21+'.info', 'w')
    with open(str2+'_ml_network.pickle') as handle:
           fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork


    infofile.write('Ground-truth\n')
    for k,v in commu.items():
        l =str(k)+' set('+ str(list(v))+')'
        infofile.write(l+'\n')
    true_labels = [None]*200
    for k, v in commu.items():
        for node in v:
            #print("node: ",node)
            true_labels[node-1] = k

    infofile.write('\n_commu_benching_all_march21_louvain.pickle\n')
    louvain_p = detectedcomms(partition_at_level(commus,len(commus)-1))
    for comms, nodelist in louvain_p.items():
        toprint = str(comms)+' '+str(nodelist)+'\n'
        infofile.write(toprint)
    
    predicted_labels = [None]*len(true_labels)
    for k, v in louvain_p.items():
        for node in v:
            #print("node: ",node)
            predicted_labels[node-1] = k
    
    nmival = normalized_mutual_info_score(true_labels, predicted_labels)
    print("NMI: ",nmival)
    nmifile = open(pathtosave+'nmis', 'a')
    nmifile.write(str2+": "+str(nmival))
    nmifile.close()

def comparisionOfModularityValues(network,modu,pathtosave):
    gtFile="./Raphael_27.6.17/infos/"+str(network)+".info"
    gtf=open(gtFile)
    gtmod=0
    for line in gtf:
        #print()
        if(len(line.split(" "))==1 and len(line.split("-"))==1 and len(line.split("_"))==1):
            gtmod = float(line.strip())
            print("gtmod "+str(gtmod))
            break
    modfile = open(pathtosave+"modComparisionDualModLouvain",'a')
    modfile.write(network+ ":    "+ str(gtmod)+"  "+str(modu)+"\n")
    modfile.close()
    
import os
import sys


# i = int(sys.argv[1])
# p = float(sys.argv[2])
# a = float(sys.argv[3])
# mu = float(sys.argv[4])
# d = float(sys.argv[5])

# print p, a, mu, d
# str1="./test"+str(i)+"/networks_alpha"+str(a)+"/"
# str2=str1+"networks_p"+str(p)+"/networks_mu"+str(mu)+"/networks_density"+str(d)+"/new_format"
# print str2
# if not os.path.exists(str2) :
#     print p, a, mu, d
#     print "No"
#     exit()
#getSeries(str2)

'''
str2 = "./nets/smallnetwork"
modu, commus = getSeries(str2)
print("Modularity: ",modu, commus)

sys.exit()

'''

pathtosave="./resultsDualModLouvain/"
modfile = open(pathtosave+"modComparisionDualModLouvain_correctmodimplementation",'w')
modfile.write("network                                   GroundTruth    Detected-Louvain\n")
modfile.close()

networklist = os.listdir('/home/user2/Downloads/sem2/mtp_prish/Louvain_mixmod/Raphael_27.6.17/synthetics')
for network in networklist:
    #str21 = "./nets/infos/"+str(network)
    str2 = "./nets/"+str(network)

    modu, commus = getSeries(str2)

    

    '''
    with open(str2+'_commu_benching_all_march21_louvain.pickle', 'wb') as handle:
                    pickle.dump(partition_at_level(commus, len(commus)-1), handle)
                with open(str2+'_modu_benching_all_march21_louvain.pickle', 'wb') as handle:
                    pickle.dump(modu, handle)
                with open(str2+'_commu_benching_frac_march21_louvain.pickle', 'wb') as handle:
                    pickle.dump(commus, handle)
    '''
    gtmod = computegtmod(str2)
    print "FINAL_MODULARITY*** ", modu , "gtmod ", gtmod
    modfile = open(pathtosave+"modComparisionDualModLouvain_correctmodimplementation",'a')
    modfile.write(network+ ":    "+ str(gtmod)+"  "+str(modu)+"\n")
    modfile.close()
    
    #Write info file for testing
    #generateinfofilefortesting(str2,str21,pathtosave)
    #comparisionOfModularityValues(network,modu,pathtosave)
