import networkx as nx
import random
import pickle
import math
import json
from sklearn.metrics import *
from aux.__modularity import  __git_modularity
from mixmod_wt.mixmod_wt_correctingImplementation_without_half import __modularity
from mixmod_wt.aux import  read_raw_network
from aux.Status import Status
from aux.show_intermediate_steps import *
from aux.__inducedGraph import __git_induced_graph
from multiprocessing import Pool
from pythonlouvainmaster import best_partition
from collections import defaultdict
from pprint import pprint
    
import os
import sys
import matplotlib.pyplot as plt
import matplotlib

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

def computegtmod(filename,weighted):
    fnetwork = 0
    #with open(filename+'_ml_network.pickle') as handle:
    #    fnetwork = pickle.load(handle)
    #ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu =read_raw_network(filename,weighted)

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

def generateinfofilefortesting(str2,str21):   
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
            print("node: ",node)
            true_labels[node-1] = k

    infofile.write('\n_commu_benching_all_march21_louvain.pickle\n')
    louvain_p = detectedcomms(partition_at_level(commus,len(commus)-1))
    for comms, nodelist in louvain_p.items():
        toprint = str(comms)+' '+str(nodelist)+'\n'
        infofile.write(toprint)
    
    predicted_labels = [None]*len(true_labels)
    for k, v in louvain_p.items():
        for node in v:
            print("node: ",node)
            predicted_labels[node-1] = k
    
    nmival = normalized_mutual_info_score(true_labels, predicted_labels)
    print("NMI: ",nmival)

    nmifile = open('nmis', 'a')
    nmifile.write(str2+": "+str(nmival))
    nmifile.close()

def plotGraph(gtmodlist,dtmodlist,x,path):
    
    plt.figure()
    matplotlib.style.use('ggplot')
    plt.plot(x,gtmodlist,'g')
    plt.plot(x,dtmodlist,'r')
    plt.title("modularity comparision with original_mod_always")
    plt.ylabel("modularity")
    plt.savefig(path+"modComparisionPlot")
    plt.clf()

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

    '''mod = __modularity(_get_commu_dict(status.node2com), status)
                status_list = list()
                level_count = 0
                __one_level(current_graph, status, status_list, level_count)
                new_mod = __modularity(_get_commu_dict(status.node2com), status)
                partition = __renumber(status.node2com)
                status_list.append(partition)
                print str(mod)+" "+str(new_mod)+" OUT"
            
                mod = new_mod
                with open(str2+'_commu_benching_all_march21_louvain_step1.pickle', 'wb') as handle:
                    pickle.dump(partition_at_level(status_list, 0), handle)
            
                current_graph = induced_graph(partition, current_graph)
                status.init(current_graph)    
                printgraph(current_graph)
                #detectedcomms(status_list)
                print "######################"
            
                communities_itr1 = status_list[0]
                print(communities_itr1)'''

    '''current_graph1 = nx.Graph()
                weightededgelist =[(e[0],e[1],e[2]) for e in current_graph.edges_iter(data='weight')]
                current_graph1.add_weighted_edges_from(weightededgelist)'''
    #print(list(current_graph1.edges(data = True)))
    dendo = best_partition(current_graph)
    communities = []
    #communities.append(communities_itr1)
    communities.extend(dendo)
    finalmodularity = __modularity(_get_commu_dict(partition_at_level(communities, len(communities)-1)), status,current_graph)
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

def getSeries(filename,weighted):
    fnetwork = 0
    #with open(filename+'_ml_network.pickle') as handle:
    #    fnetwork = pickle.load(handle)
    #ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    #fc = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = read_raw_network(filename,weighted)
    dendogram, mod = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)
    return mod, dendogram
    #print(fc)

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

    #networklist = os.listdir('./Raphael_27.6.17/net_networks')
    #pathtosave ="./resultsOriginalModAlways/"
    #modfilename = pathtosave+"modComparisionOriginalLouvainAlways_newnets_oldmod_with_cionly"

def runformanynetworks(args):
    output = []

    networklist  = args[:-1]
    tmp = args[-1]
    weighted = args[-1][0]
    modfilename = args[-1][1]
    networkpath = args[-1][2]
    
    for network in networklist:
        #str2 = "./nets/"+str(network)
        #str2 = "./Raphael_27.6.17/net_networks/" +str(network)
        #str2 = './syntheticNetworkGeneration/netsForDtDmDb/_networks/netsByGenerateNetsv3/alpha0.7/' + str(network)
        str2 = networkpath + str(network)
        
        dtmod, dtcom = getSeries(str2,weighted)
        gtmod,gtcom = computegtmod(str2,weighted)
        output.append((gtmod, dtmod))
        #print(str2+ ", " + str(gtmod) + ", " + str(dtmod))
        #break
        modfile = open(modfilename, 'a')
        modfile.write(str2+ ":    "+ str(gtmod)+"  "+str(dtmod)+"\n")
        modfile.close()


        #nmifile = open(nmifilename,'a')
        #dtcomdict =_get_com_wise_nodes(partition_at_level(dtcom,len(dtcom)-1))
        #print(type(dtcomdict))
        #nmifile.write(str2+":   "+generate_nmi_file(gtcom,dtcomdict))
        #nmifile.close()

        #detected_commu = _get_com_wise_nodes(partition_at_level(dtcom, len(dtcom)-1))
        #write_commus_infile(network,detected_commu,gtcom)
    return output

def parallelimplementation(networklist,weighted,modfilename,networkpath):
    #Open file to compare modularity values
    modfile = open(modfilename,'w')
    modfile.write("network                                   GroundTruth    Detected-Louvain\n")
    modfile.close()

    #nmifile = open(nmifilename,'a')
    #nmifile.write("network                                   NMI ")
    #nmifile.close()

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
        arguments = networklist[splits[i]:splits[i+1]]
        tmp = [weighted,modfilename,networkpath]
        arguments.append(tmp)
        args.append(arguments)

    p = Pool(cores)
    modularities = p.map(runformanynetworks, args)

    #Flatten list of lists returned by different cores
    modularities = [item for items in modularities for item in items]

    print modularities
    return modularities

def main():
    weighted=0
    #pathtosave = "./syntheticNetworkGeneration/results/nets121_incWeights/alpha0.7/dt_db_mod_files_corrctImp_without_half/"
    pathtosave = "./syntheticNetworkGeneration/results/nets121/alpha0.7/dt_db_mod_files_corrctImp_without_half/"
    #networkpath =  './syntheticNetworkGeneration/netsForDtDmDb/_networks/netsByGenerateNetsv3/alpha0.7/'
    networkpath = './syntheticNetworkGeneration/netsForDtDmDb/_networks/morenets/alpha0.7/'

    networklist = os.listdir(networkpath)
    modfilename = pathtosave+"all_louvain_for_corrImpl"

    gtmodlist=[]
    dtmodlist=[]
    x=[]
    i=1
    print(len(networklist))
    modfile = open(modfilename,'w')
    modfile.write("network                                   GroundTruth    Detected-Louvain\n")
    modfile.close()
    runformanynetworks(networklist)
    #parallelimplementation(networklist)

    sys.exit()


    for network in networklist:
        str21 = "./nets/infos/"+str(network)
        str2 = "./nets/"+str(network)
        print(network)
        modu, commus = getSeries(str2)

        print "FINAL_MODULARITY*** ", modu

        '''with open(str2+'_commu_benching_all_march21_louvain.pickle', 'wb') as handle:
                        pickle.dump(partition_at_level(commus, len(commus)-1), handle)
                    with open(str2+'_modu_benching_all_march21_louvain.pickle', 'wb') as handle:
                        pickle.dump(modu, handle)
                    with open(str2+'_commu_benching_frac_march21_louvain.pickle', 'wb') as handle:
                        pickle.dump(commus, handle)'''

        #Write info file for testing
        #generateinfofilefortesting(str2,str21)
        dtmodlist.append(modu)
        x.append(i)
        i=i+1
        gtFile="./Raphael_27.6.17/infos/"+str(network)+".info"
        gtf=open(gtFile)
        gtmod=0
        for line in gtf:
            #print()
            if(len(line.split(" "))==1 and len(line.split("-"))==1 and len(line.split("_"))==1):
                gtmod = float(line.strip())
                print("gtmod "+str(gtmod))
                gtmodlist.append(gtmod)
                break
        modfile = open(pathtosave+"modComparisionOriginalLouvainAlways",'a')
        modfile.write(str2+ ":    "+ str(gtmod)+"  "+str(modu)+"\n")
        modfile.close()

    plotGraph(gtmodlist,dtmodlist,x,pathtosave)


'''
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
    p_temp = __renumber(status.node2com)
    status_list.append(p_temp)
    cur_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
    status_list.pop()
    new_mod = cur_mod

    print "# id_node from_com to_com local_mod mod"

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
                status.node2com[node] = com
                
                p_temp = __renumber(status.node2com)
                status_list.append(p_temp)
                incr =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status) - cur_mod2
                status_list.pop()

                if incr > best_increase :
                    best_increase = incr
                    best_com = com

                status.node2com[node] = -1
            status.node2com[node] = best_com

            if(verbose):
                pass #print("Node ", node, " moved from ", com_node, " to ", best_com, " increase: ", best_increase)
            
            p_temp = __renumber(status.node2com)
            status_list.append(p_temp)
            cur_mod2 =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status, verbose)
            status_list.pop()
            
            if best_com != com_node :
                modif = True
            
                p_temp2 = __renumber(status.node2com)
                status_list.append(p_temp2)
                incr =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status) - cur_mod2
                #print "m here"
                #print node, com_node, best_com, incr, best_increase, __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status), cur_mod2
                status_list.pop()

        p_temp = __renumber(status.node2com)
        status_list.append(p_temp)
        
        #Changing mod to optimized modularity definition here
        new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
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
'''

if __name__ == '__main__':
    main()

    
