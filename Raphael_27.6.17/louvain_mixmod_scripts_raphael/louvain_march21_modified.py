import networkx as nx
import random
import pickle
import math

__PASS_MAX = -1
__MIN = 0.0000001

def __modularity(commu, status):

    layer=status.layer
    node_l=status.node_l
    node_c=status.node_c       
    top=status.top
    bot=status.bot
    edge_l=status.edge_l
    edge_c=status.edge_c
    couple=status.couple
    mu = status.mu

    f=0 
    intra_inter={}
    for c in commu:
        intra_inter[c]=set()
        
        for n in commu[c]:
            for l in layer:
                if n in layer[l]:
                    intra_inter[c].add(l)
    
    #print intra_inter
    
    modularity=0    
    x1={}
    x2={}   
    for c in commu:
        x1[c]=0
        x2[c]=0
        modc_layer=0
        
        if len(intra_inter[c])>1:
            for l in layer:
                d_layer=0
                I_layer=0
                m_layer=0
                n_layer=0
                n_co_com_layer=0
            
                for n in layer[l]:
                    if n in node_l:
                        m_layer+=len(node_l[n])
                    #if n in node_c:
                    #   m_layer+=len(node_c[n]) 
                    if n in commu[c]:    #if the node belongs to current community
                        n_layer+=1
                        if n in node_l:
                            d_layer+=len(node_l[n])
                    
                            for nei in node_l[n]:
                                if nei in commu[c]: #if the neighbour belongs to current community
                                    I_layer+=1
                    
                        if n in node_c:
                            for nei in node_c[n]:
                                if nei not in commu[c]: #connected to atleast one crosslayer node outside community
                                    n_co_com_layer+=1 
                                    break
                        '''
                            for nei in node_c[n]:
                                if nei in commu[c]: #if the neighbour belongs to current community
                                    I_layer+=1
                        '''
                                            
                I_layer=float(I_layer)/2.0
                m_layer=float(m_layer)/2.0
                if I_layer>-1 and m_layer>0:
                    mod=((I_layer/m_layer)-((float(d_layer)/(2*m_layer))*(float(d_layer)/(2*m_layer))))
                else:
                    mod=0   
                #if n_layer > 0:
                #   mod=mod*pow(2.718,-(float(n_co_com_layer)/float(n_layer)))
            
                #print mod  
                #print c,l,mod,mu,mu*mod,I_layer,m_layer,d_layer
                if f==1:
                    modc_layer+=mu*mod
                else:   
                    modc_layer+=mod
                    x1[c]+=mod
            
            modc_couple=0
            for co in couple:
                d_couple_top=0.0
                d_couple_bot=0.0
                I_couple=0.0
                m_couple=0.0
            
                top_tot=0.0
                top_con=0.0
                bot_tot=0.0
                bot_con=0.0
                
                m_layer_top=0.0
                m_layer_bot=0.0
                
                for n in couple[co]:
                    if n in node_c:
                        for nei in node_c[n]:
                            if (n in layer[top[co]] and nei in layer[bot[co]]) or (n in layer[bot[co]] and nei in layer[top[co]]):
                                m_couple+=1
                
                    if n in layer[top[co]]: #n belongs to the top layer of coupling
                        if n in node_l:
                            m_layer_top+=len(node_l[n])
                            
                        if n in commu[c]:    #if the node belongs to current community
                            top_tot+=1
                            if n in node_c:
                                flagg=0
                                for nei in node_c[n]:
                                    if nei in layer[bot[co]]:
                                        d_couple_top+=1
                                        if nei in commu[c]: #if the neighbour belongs to current community
                                            I_couple+=1
                                            flagg=1
                                if flagg==1: #connected to at least 1 within community node in bottom layer
                                    top_con+=1
                            if n in node_l and n not in node_c:
                                d_couple_top+=len(node_l[n])
                                    
                                            
                    if n in layer[bot[co]]: #n belongs to the bottom layer of coupling
                        if n in node_l:
                            m_layer_bot+=len(node_l[n])
                        if n in commu[c]:    #if the node belongs to current community
                            bot_tot+=1
                            if n in node_c:
                                flagg=0
                                for nei in node_c[n]:
                                    if nei in layer[top[co]]:
                                        d_couple_bot+=1 
                                        if nei in commu[c]: #if the neighbour belongs to current community
                                            flagg=1                 
                                                #break
                                if flagg==1: #connected to at least 1 within community node in bottom layer
                                    bot_con+=1
                            if n in node_l and n not in node_c:
                                d_couple_bot+=len(node_l[n])        
                                
                I_couple=float(I_couple)
                m_couple=float(m_couple)/2.0
                m_layer_top=float(m_layer_top)/2.0
                m_layer_bot=float(m_layer_bot)/2.0
                
                if I_couple>-1 and (m_couple+m_layer_top+m_layer_bot)>0:
                    
                    mod=((I_couple/(m_couple+m_layer_top+m_layer_bot))-((d_couple_top*d_couple_bot)/((m_couple+2.0*m_layer_top+2.0*m_layer_bot)*(m_couple+2.0*m_layer_top+2.0*m_layer_bot))))
                else:
                    mod=0   
                #print I_couple, m_couple, d_couple_top, d_couple_bot   
                #print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
                if f==1:
                    modc_couple+=2*(1-mu)*mod
                else:
                    modc_couple+=mod
                    x2[c]+=mod
                #print modc_couple           
                #print "ha hh"
            modularity+=modc_layer+modc_couple
            
        else:
            for l in layer:
                d_layer=0
                I_layer=0
                m_layer=0
                n_layer=0
                n_co_com_layer=0
            
                for n in layer[l]:
                    if n in node_l:
                        m_layer+=len(node_l[n])
                
                c_layer=0.0 
                for n in layer[l]:
                    if n in node_c:
                        c_layer+=len(node_c[n])
                            
                    #if n in node_c:
                    #   m_layer+=len(node_c[n]) 
                    if n in commu[c]:    #if the node belongs to current community
                        n_layer+=1
                        if n in node_l:
                            d_layer+=len(node_l[n])
                    
                            for nei in node_l[n]:
                                if nei in commu[c]: #if the neighbour belongs to current community
                                    I_layer+=1
                    
                        if n in node_c and len(node_c[n])>0:
                            n_co_com_layer+=1 
                            d_layer+=len(node_c[n])
                            
                        
                                            
                I_layer=float(I_layer)/2.0
                m_layer=float(m_layer)/2.0
                
                if I_layer>-1 and (1.0*m_layer+(c_layer/2.0))>0:# and m_layer1>0:
                    
                    mod=((I_layer/(1.0*m_layer+(c_layer/2.0)))-((float(d_layer)/(2.0*m_layer+c_layer))*(float(d_layer)/(2.0*m_layer+c_layer))))
                    #else:
                    #mod=(-((d_layer/(2*m_layer))*(d_layer/(2*m_layer))))    
                else:
                    mod=0 
                #if n_layer > 0:
                #   mod=mod*pow(2.718,-((float(n_co_com_layer))/float(n_layer)))
            
                #print mod  
                #print c,l,mod,mu,mu*mod,I_layer,m_layer,c_layer,d_layer, n_co_com_layer, n_layer
                if f==1:
                    modc_layer+=mu*mod
                else:   
                    modc_layer+=mod
                    x1[c]+=mod
            modularity+=modc_layer      
                                    
    #print x1,x2 
    return 0.333*modularity
      
def getModularityQ(commu, status):

    layer=status.layer
    node_l=status.node_l
    node_c=status.node_c       
    top=status.top
    bot=status.bot
    edge_l=status.edge_l
    edge_c=status.edge_c
    couple=status.couple
    mu = status.mu

    f=0
    modularity=0        
    for c in commu:
        modc_layer=0
        for l in layer:
            d_layer=0
            I_layer=0
            m_layer=0
            
            for n in layer[l]:
                if n in node_l:
                    m_layer+=len(node_l[n])
                if n in commu[c]:    #if the node belongs to current community
                    if n in node_l:
                        d_layer+=len(node_l[n])
                    
                        for nei in node_l[n]:
                            if nei in commu[c]: #if the neighbour belongs to current community
                                I_layer+=1
            I_layer=float(I_layer)/2.0
            m_layer=float(m_layer)/2.0
            if I_layer>-1:
                mod=((I_layer/m_layer)-((d_layer/(2*m_layer))*(d_layer/(2*m_layer))))
            else:
                mod=0   
            #print c,l,mod,mu,mu*mod,I_layer,m_layer,d_layer
            if f==1:
                modc_layer+=mu*mod
            else:   
                modc_layer+=mod
        modc_couple=0
        for co in couple:
            d_couple_top=0
            d_couple_bot=0
            I_couple=0
            m_couple=0
            
            for n in couple[co]:
                if n in node_c:
                    for nei in node_c[n]:
                        if (n in layer[top[co]] and nei in layer[bot[co]]) or (n in layer[bot[co]] and nei in layer[top[co]]):
                            m_couple+=1
                
                if n in layer[top[co]]: #n belongs to the top layer of coupling
                    if n in commu[c]:    #if the node belongs to current community
                        if n in node_c:
                            for nei in node_c[n]:
                                if nei in layer[bot[co]]:
                                    d_couple_top+=1
                                    if nei in commu[c]: #if the neighbour belongs to current community
                                        I_couple+=1
                if n in layer[bot[co]]: #n belongs to the bottom layer of coupling
                    if n in commu[c]:    #if the node belongs to current community
                        if n in node_c:
                            for nei in node_c[n]:
                                if nei in layer[top[co]]:
                                    d_couple_bot+=1                     
                                    
            I_couple=float(I_couple)
            m_couple=float(m_couple)/2.0
            if I_couple>-1 and m_couple!=0:
                mod=((I_couple/m_couple)-((d_couple_top*d_couple_bot)/((m_couple)*(m_couple))))
            else:
                mod=0   
            #print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
            if f==1:
                modc_couple+=2*(1-mu)*mod
            else:
                modc_couple+=mod        
        
        modularity+=modc_layer+modc_couple                  
    #print len(layer)
    #print len(couple)  
    return (1/float(len(layer)+len(couple)))*modularity 

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

class Status :
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}

    layer={}
    node_l={}
    node_c={}       
    top={}
    bot={}
    edge_l={}
    edge_c={}
    couple={}
    mu = 0

    def __init__(self) :
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])

        self.layer=dict([])
        self.node_l=dict([])
        self.node_c=dict([])
        self.top=dict([])
        self.bot=dict([])
        self.edge_l=dict([])
        self.edge_c=dict([])
        self.couple=dict([])
        self.mu = 0

    def __str__(self) :
        return ("node2com : " + str(self.node2com) + " degrees : "
            + str(self.degrees) + " internals : " + str(self.internals)
            + " total_weight : " + str(self.total_weight)) 

    def copy(self) :
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = self.node2com.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.total_weight = self.total_weight
        new_status.layer=self.layer.copy()
        new_status.node_l=self.node_l.copy()
        new_status.node_c=self.node_c.copy()
        new_status.top=self.top.copy()
        new_status.bot=self.bot.copy()
        new_status.edge_l=self.edge_l.copy()
        new_status.edge_c=self.edge_c.copy()
        new_status.couple=self.couple.copy()
        new_status.mu = self.mu

    def init(self, graph, part = None) :
        """Initialize the status of a graph with every node in one community"""
        count = 0
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])

        self.total_weight = graph.size(weight = 'weight')
        for node in graph.nodes() :
            self.node2com[node] = count
            deg = float(graph.degree(node, weight = 'weight'))
            if deg < 0 :
                raise ValueError("Bad graph type, use positive weights")
            self.degrees[count] = deg
            self.gdegrees[node] = deg
            self.loops[node] = float(graph.get_edge_data(node, node,
                                             {"weight":0}).get("weight", 1))
            self.internals[count] = self.loops[node]
            count = count + 1

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
            
            p_temp = __renumber(status.node2com)
            status_list.append(p_temp)
            cur_mod2 =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
            status_list.pop()

            if best_com != com_node :
                modif = True
                
                '''p_temp2 = __renumber(status.node2com)
                status_list.append(p_temp2)
                incr =  __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status) - cur_mod2
                
                print node, com_node, best_com, incr, best_increase, __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status), cur_mod2
                status_list.pop()'''

        p_temp = __renumber(status.node2com)
        status_list.append(p_temp)
        new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
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
    #print str(mod)+" "+str(new_mod)+" OUT"
    
    mod = new_mod
    with open(str2+'_commu_benching_all_march21_louvain_step1.pickle', 'wb') as handle:
        pickle.dump(partition_at_level(status_list, 0), handle)
    
    current_graph = induced_graph(partition, current_graph)
    status.init(current_graph)    
    print "######################"
    while True :
        level_count+=1
        __one_level(current_graph, status, status_list, level_count)
        partition = __renumber(status.node2com)
        status_list.append(partition)
        new_mod = __modularity(_get_commu_dict(partition_at_level(status_list, level_count)), status)
        
        #print str(mod)+" "+str(new_mod)+" IN"
        if new_mod - mod < __MIN :
            break
        mod = new_mod
        current_graph = induced_graph(partition, current_graph)
        status.init(current_graph)
        print "######################"

    return status_list[:-1], mod

def getSeries(filename):
    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork

    print "LDELDEL"
    print couple
    sys.exit(1)

    dendogram, mod = louvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu)
    return mod, dendogram
    
import os
import sys

str2 = sys.argv[1]

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
modu, commus = getSeries(str2)

#print "FINAL_MODULARITY*** ", modu

with open(str2+'_commu_benching_all_march21_louvain.pickle', 'wb') as handle:
    pickle.dump(partition_at_level(commus, len(commus)-1), handle)
with open(str2+'_modu_benching_all_march21_louvain.pickle', 'wb') as handle:
    pickle.dump(modu, handle)
with open(str2+'_commu_benching_frac_march21_louvain.pickle', 'wb') as handle:
    pickle.dump(commus, handle)
#print ""


