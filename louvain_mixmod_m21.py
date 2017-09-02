import networkx as nx
import random
import pickle
import math
from mixmod.Status import Status
from mixmod.modularity import __modularity

__PASS_MAX = -1
__MIN = 0.0000001

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

def __one_level(graph, status, status_list, level_count) :
    modif = True
    nb_pass_done = 0
    p_temp = __renumber(status.node2com)
    status_list.append(p_temp)
    cur_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
    status_list.pop()
    new_mod = cur_mod

    while modif  and nb_pass_done != __PASS_MAX :
        cur_mod = new_mod
        modif = False
        nb_pass_done += 1

        for node in graph.nodes() :
            com_node = status.node2com[node]
            neigh_communities = __neighcom(node, graph, status)
            status.node2com[node] = -1
            best_com = com_node
            best_increase = 0
            for com in neigh_communities:
                status.node2com[node] = com
                
                p_temp = __renumber(status.node2com)
                status_list.append(p_temp)
                incr =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status) - cur_mod
                status_list.pop()

                if incr > best_increase :
                    best_increase = incr
                    best_com = com

                status.node2com[node] = -1

            status.node2com[node] = best_com
            
            if best_com != com_node :
                modif = True
        
        p_temp = __renumber(status.node2com)
        status_list.append(p_temp)
        new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
        status_list.pop()
        if new_mod - cur_mod < __MIN :
            break

def induced_graph_multilayer(partition, graph, layer):
    ret = nx.Graph()
    id_extra_com = len(partition) + 1

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

            part[id_com] = id_com

            for id_node in layer_node:
                part[id_node] = id_com
                list_node_com[id_node] = id_com

    for node1, node2, datas in graph.edges_iter(data = True):
        weight = datas.get("weight", 1)
        com1 = list_node_com[node1]
        com2 = list_node_com[node2]
        w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
        ret.add_edge(com1, com2, weight = w_prec + weight)

    return ret, part

def induced_graph(partition, graph) :
    ret = nx.Graph()
    ret.add_nodes_from(partition.values())

    for node1, node2, datas in graph.edges_iter(data = True):
        weight = datas.get("weight", 1)
        com1 = partition[node1]
        com2 = partition[node2]
        w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
        ret.add_edge(com1, com2, weight = w_prec + weight)

    return ret

def louvain(graph, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu) :
    current_graph = graph.copy()
    status = Status()
    status.init(current_graph)

    status.layer=layer
    status.node_l=node_l
    status.node_c=node_c
    status.top=top
    status.bot=bot
    status.edge_l=edge_l
    status.edge_c=edge_c
    status.couple = couple
    status.mu = mu

    mod = __modularity(_get_commu_dict(status.node2com), status)
    status_list = list()
    level_count = 0
    __one_level(current_graph, status, status_list, level_count)
    new_mod = __modularity(_get_commu_dict(status.node2com), status)
    partition = __renumber(status.node2com)
    status_list.append(partition)
    ##print str(mod)+" "+str(new_mod)+" OUT"
    mod = new_mod
    
    current_graph = induced_graph(partition, current_graph)
    #current_graph, part = induced_graph_multilayer(partition, current_graph, layer)
    
    status.init(current_graph)

    while True :
        level_count+=1
        #print level_count
        __one_level(current_graph, status, status_list, level_count)
        
        partition = __renumber(status.node2com)
        status_list.append(partition)
        new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
        ##print 'new_mod2',new_mod
        #print str(mod)+" "+str(new_mod)+" IN"
        if new_mod - mod < __MIN :
            break
        mod = new_mod
        #current_graph = induced_graph(partition, current_graph)
        current_graph, part = induced_graph_multilayer(partition, current_graph, layer)
        #status.init(current_graph)
        status.init(current_graph, part)

    return status_list[:-1], mod

def getSeries(filename):
    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    dendogram, mod = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)
    return mod, dendogram
    
import os
import sys

#Comment following four lines if you want to run for all networks
str2 = "./nets/network_0.9_1.0_0.05_1.0_0.0"
modu, commus = getSeries(str2)
print("Modularity: ", modu, commus)

sys.exit()


'''
with open('_commu_benching_all_march21_louvain_mixmod.pickle', 'wb') as handle:
    pickle.dump(partition_at_level(commus, len(commus)-1), handle)


#str2 = sys.argv[1]
pathtosave = './resultsmixmod/'
modfile = open(pathtosave+"modComparisionMixModLouvain_nomultilayerinstage1",'w')
modfile.write("network                                   GroundTruth    DetectedOriginalLouvain\n")
modfile.close()

networklist = os.listdir('./Raphael_27.6.17/synthetics')
for network in networklist:
    str2 = "./nets/"+str(network)
    #str21 = '/home/user/Downloads/sem2/mtp_prish/Louvain_mixmod/resultsmixmod/modComparisionMixModLouvain'
    modu, commus = getSeries(str2)
    print "FINAL_MODULARITY*** ", modu
    gtFile="./Raphael_27.6.17/infos/"+str(network)+".info"
    gtf=open(gtFile)
    gtmod=0
    for line in gtf:
        ##print()
        if(len(line.split(" "))==1 and len(line.split("-"))==1 and len(line.split("_"))==1):
            gtmod = float(line.strip())
            #print("gtmod "+str(gtmod))
            break

    modfile = open(pathtosave+"modComparisionMixModLouvain_nomultilayerinstage1",'a')
    modfile.write(str2+ ":    "+ str(gtmod)+"  "+str(modu)+"\n")
    modfile.close()

'''


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
