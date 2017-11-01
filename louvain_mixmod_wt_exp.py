import networkx as nx
import random
import pickle
import math
from mixmod_wt.Status import Status
#from mixmod_wt.new_modularity import __modularity
from mixmod_wt.mixmod_wt_correctimp_exp3 import __modularity
#from mixmod_wt.aux import  __modularity
from collections import defaultdict
import matplotlib.pyplot as plt
from mixmod_wt.aux import _get_com_wise_nodes, printsomeinformation
from multiprocessing import Pool

__PASS_MAX = -1
__MIN = 0.0000001
#modctr = 0


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
    count = 1
    #count = 0
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


#New __one_level function (from Raphael's code)
def __one_level(graph, status, status_list, level_count, verbose=0) :
    #print("graph edges: ",graph.edges(data = True))
    modif = True
    nb_pass_done = 0
    #p_temp = __renumber(status.node2com)
    p_temp = status.node2com
    
    status_list.append(p_temp)
    #cur_mod = __modularity(_get_commu_dict(status_list[-1]), status, graph)
    cur_mod = __modularity(_get_com_wise_nodes(status_list[-1]), status, graph)
    
    status_list.pop()
    new_mod = cur_mod
    #print "# id_node from_com to_com local_mod mod"
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

                temp_dict = {com:_get_com_wise_nodes(status.node2com)[com]}
                base_mod_of_community = __modularity(temp_dict, status, graph)

                status.node2com[node] = com

                #p_temp = __renumber(status.node2com)
                p_temp = status.node2com
                status_list.append(p_temp)
                
                #incr = __modularity(_get_commu_dict(status_list[-1]), status, graph) - cur_mod2

                temp_dict = {com:_get_com_wise_nodes(status_list[-1])[com]}
                incr = __modularity(temp_dict, status, graph) - base_mod_of_community
                
                status_list.pop()

                if incr > best_increase :
                    best_increase = incr
                    best_com = com

                if(verbose): printsomeinformation(node, com, best_com, incr, status.node_l, status.node_c)

                status.node2com[node] = -1

            status.node2com[node] = best_com
            
            #p_temp = __renumber(status.node2com)
            p_temp = status.node2com
            status_list.append(p_temp)
            #cur_mod2 =  __modularity(_get_commu_dict(status_list[-1]), status, graph)
            cur_mod2 =  __modularity(_get_com_wise_nodes(status_list[-1]), status, graph)
            status_list.pop()

            if best_com != com_node :
                modif = True
            
            if(verbose):
                pass
                #print("{0} {1} {2}".format(node, com_node, best_com))
                #nx.draw(graph, with_labels = True)
                #plt.show()

                '''p_temp2 = __renumber(status.node2com)
                status_list.append(p_temp2)
                incr =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status) - cur_mod2
                
                print node, com_node, best_com, incr, best_increase, __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status), cur_mod2
                status_list.pop()'''

        #p_temp = __renumber(status.node2com)
        p_temp = status.node2com
        status_list.append(p_temp)
        #new_mod = __modularity(_get_commu_dict(status_list[-1]), status, graph)
        new_mod = __modularity(_get_com_wise_nodes(status_list[-1]), status, graph)
        
        print("In __one_level new_mod: {0:.4f} cur_mod: {1:.4f}".format(new_mod,cur_mod))
        if(verbose): print("Status list[-1]: ",status_list[-1])
        
        status_list.pop()
        if new_mod - cur_mod < __MIN :
            break

#__modularity(_get_commu_dict(status_list[-1]), status)


def induced_graph_multilayer(partition, graph, status):
    new_layer =defaultdict(set)
    new_node_l=defaultdict(set)
    new_node_c=defaultdict(set)
    new_couple=defaultdict(set)
    layer = status.layer

    ret = nx.Graph()
    #id_extra_com = len(partition) + 1
    id_extra_com = max(partition.values()) + 1

    list_node_com = {}
    part = {}

    partition_rebuild = {}
    for id_node in partition:
        partition_rebuild.setdefault(partition[id_node], [])
        partition_rebuild[partition[id_node]].append(id_node)
    
    for id_com in partition_rebuild.keys():
        
        layer_node = {}
        is_top = False
        is_bot = False 
        for id_node in partition_rebuild[id_com]:
            if id_node in layer[1]:
                layer_node[id_node] = 1
                is_top = True
            else:
                layer_node[id_node] = 2
                is_bot = True

        if is_top and is_bot: # add two nodes into induced graph
            ret.add_node(id_com)
            ret.add_node(id_extra_com)
            #print("idcom: ", id_com, "idextracom: ",id_extra_com)
            #updating status
            new_layer[1].add(id_com)
            new_layer[2].add(id_extra_com)
            ###

            part[id_com] = id_com
            part[id_extra_com] = id_com # id_extra_com will be remain with id_com community

            for id_node in layer_node:
                if layer_node[id_node] == 1:
                    list_node_com[id_node] = id_com
                else:
                    list_node_com[id_node] = id_extra_com

            id_extra_com += 1

        else: # add one node into induced graph
            ret.add_node(id_com)
            #print("idcom: ", id_com)
            #updating status
            if(is_top): new_layer[1].add(id_com)
            elif(is_bot): new_layer[2].add(id_com)
            #######

            part[id_com] = id_com

            for id_node in layer_node:
                #part[id_node] = id_com
                list_node_com[id_node] = id_com

    for node1, node2, datas in graph.edges_iter(data = True):
        weight = datas.get("weight", 1)
        com1 = list_node_com[node1]
        com2 = list_node_com[node2]
        if(com1==com2):             # for selfloops add double edge weights
            weight *= 2
        w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
        ret.add_edge(com1, com2, weight = w_prec + weight)

    #updating status
    for node1,node2 in ret.edges_iter():
        if((node1 in new_layer[1] and node2 in new_layer[1]) or (node1 in new_layer[2] and node2 in new_layer[2])):
            #add to node_l
            new_node_l[node1].add(node2)
            new_node_l[node2].add(node1)
        else:
            #add to node_c
            new_node_c[node1].add(node2)
            new_node_c[node2].add(node1)

    #updating status
    new_couple[1]=set(ret.nodes())
        
    status.layer = new_layer
    status.couple = new_couple
    status.node_c = new_node_c
    status.node_l = new_node_l
    return ret, part,status

def louvain(graph, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu) :
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

    status.init(current_graph)

    mod = __modularity(_get_commu_dict(status.node2com), status, current_graph)
    status_list = list()
    level_count = 0
    
    __one_level(current_graph, status, status_list, level_count)
    #new_mod = __modularity(_get_commu_dict(status.node2com), status, current_graph)
    new_mod = __modularity(_get_com_wise_nodes(status.node2com), status, current_graph)
    partition = __renumber(status.node2com)
    
    status_list.append(partition)
    ##print str(mod)+" "+str(new_mod)+" OUT"
    mod = new_mod
    current_graph,part,status = induced_graph_multilayer(partition, current_graph,status)
    #print("Louvain, partition: ",partition)
    status.init(current_graph)

    #print("status.layer: ",status.layer)

    #sys.exit()

    while True :
        level_count+=1
        #print level_count
        __one_level(current_graph, status, status_list, level_count, 1)
        
        partition = __renumber(status.node2com)
        #partition = status.node2com
        status_list.append(partition)

        #new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)  
        
        #new_mod = __modularity(_get_commu_dict(partition), status, current_graph)
        new_mod = __modularity(_get_com_wise_nodes(partition), status, current_graph)
        
        if new_mod - mod < __MIN :
            print("In Louvain new_mod: {0:.4f} mod = {1:.4f} AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".format(new_mod, mod))
            break
        mod = new_mod
        #current_graph = induced_graph(partition, current_graph)
        current_graph, part,status = induced_graph_multilayer(partition, current_graph,status)
        #status.init(current_graph)
        status.init(current_graph, part)

        
    return status_list[:-1], mod

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
    return mod,commu

def getSeries(filename):
    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    dendogram, mod = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)
    return mod, dendogram

def write_commus_infile(network,detected_commu,gtcom):
    #writting the detected communities-------------------------------
    detected_communities_file = open(pathtowritecommu+str(network)+".pickle","w")
    pickle.dump((gtcom,detected_commu),detected_communities_file,2)
    detected_communities_file.close()
    #-----------------------------------------------------------------

    
import os
import sys
import pickle

networklist = os.listdir('./Raphael_27.6.17/synthetics')
pathtowritecommu = "./resultsmixmod/detected_communities/"
pathtosave = './resultsmixmod/modularity_comparisions/'
modfilename = pathtosave+"modComparisionMixModLouvain_wt_correctedImplementation_exp4"


'''
#Comment following four lines if you want to run for all networks
str2 = "./nets/network_0.2_1.0_0.05_1.0_0.0"
#str2 = "./nets/smallnetwork"
modu, commus = getSeries(str2)
print("Modularity: ", modu, "Communities: ",_get_com_wise_nodes(partition_at_level(commus, len(commus)-1)))
print("GT Mod: ",computegtmod(str2))

sys.exit()
'''

def runformanynetworks(networklist):
    output = []
    for network in networklist:
        str2 = "./nets/"+str(network)
        dtmod, dtcom = getSeries(str2)
        gtmod,gtcom = computegtmod(str2)
        output.append((gtmod, dtmod))
        modfile = open(modfilename, 'a')
        modfile.write(str2+ ":    "+ str(gtmod)+"  "+str(dtmod)+"\n")
        modfile.close()
        #detected_commu = _get_com_wise_nodes(partition_at_level(dtcom, len(dtcom)-1))
        #write_commus_infile(network,detected_commu,gtcom)
    return output


def parallelimplementation(networklist):
    #Open file to compare modularity values
    modfile = open(modfilename,'w')
    modfile.write("network                                   GroundTruth    Detected-Louvain\n")
    modfile.close()

    #Prepare args for parallel processings
    numnetworks = len(networklist)

    cores=4
    chunksize = numnetworks/cores

    splits = []
    for i in range(cores):
        splits.append((i)*chunksize)
    splits.append(numnetworks)

    args = []
    for i in range(cores):
        args.append(networklist[splits[i]:splits[i+1]])

    p = Pool(cores)
    modularities = p.map(runformanynetworks, args)

    #Flatten list of lists returned by different cores
    modularities = [item for items in modularities for item in items]

    print modularities
    return modularities



parallelimplementation(networklist)

sys.exit()


# i = int(sys.argv[1])
# p = float(sys.argv[2])
# a = float(sys.argv[3])
# mu = float(sys.argv[4])
# d = float(sys.argv[5])

# #print p, a, mu, d
# str1="./test"+str(i)+"/networks_alpha"+str(a)+"/"
# str2=str1+"networks_p"+str(p)+"/networks_mu"+str(mu)+"/networks_density"+str(d)+"/new_format"
# #print str2
# if not os.path.exists(str2) :
#     #print p, a, mu, d
#     #print "No"
#     exit()

#modu, commus = getSeries(str2)

##print "FINAL_MODULARITY*** ", modu

'''with open(str2+'_commu_benching_all_march21_louvain_mixmod.pickle', 'wb') as handle:
    pickle.dump(partition_at_level(commus, len(commus)-1), handle)
with open(str2+'_modu_benching_all_march21_louvain_mixmod.pickle', 'wb') as handle:
    pickle.dump(modu, handle)
with open(str2+'_commu_benching_frac_march21_louvain_mixmod.pickle', 'wb') as handle:
    pickle.dump(commus, handle)
#print ""'''
