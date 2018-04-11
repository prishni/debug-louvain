import networkx as nx
import random
import pickle
import math

def getModularityQ_modified(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c):
    f=0 
    modularity=0    
    x1={}
    x2={}   
    for c in commu:
        x1[c]=0
        x2[c]=0
        modc_layer=0
        '''
        clus_single=0
        clus_multi=0
        
        for l in layer:
            for n in layer[l]:
                if n in commu[c]:
                    clus1=0
                    clus2=0
                    if n in node_l:
                        for nei1 in node_l[n]:
                            for nei2 in node_l[n]:
                                if nei2>nei1 and nei1 in node_l and nei2 in node_l[nei1]:
                                    clus1+=1
                                    clus2+=1
                        if len(node_l[n])>1:
                            clus_single+=float(2*clus1)/(float)((len(node_l[n]))*(len(node_l[n])-1))            
                                        
                    if n in node_c:
                        for nei1 in node_c[n]:
                            for nei2 in node_c[n]:
                                if nei2>nei1 and nei1 in node_l and nei2 in node_l[nei1]:
                                    clus2+=1
                    if n in node_c and n in node_l and (len(node_c[n])+len(node_l[n]))>1:
                        deno=(len(node_c[n])+len(node_l[n]))
                        clus_multi+=float(2*clus2)/(float)(deno*(deno-1))
                    else:
                        if n not in node_l and n in node_c and (len(node_c[n]))>1:
                            deno=len(node_c[n])             
                            clus_multi+=float(2*clus2)/(float)(deno*(deno-1))
                        else:
                            if n not in node_c and n in node_l and (len(node_l[n]))>1:
                                deno=len(node_l[n])             
                                clus_multi+=float(2*clus2)/(float)(deno*(deno-1))
        
        clus_multi=clus_multi/float(len(commu[c]))
        clus_single=clus_single/float(len(commu[c]))
        '''
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
            if I_layer>-1:
                mod=((I_layer/m_layer)-((d_layer/(2*m_layer))*(d_layer/(2*m_layer))))
            else:
                mod=0   
            if n_layer > 0:
                mod=mod*pow(2.718,-(float(n_co_com_layer)/float(n_layer)))
            
            #print mod  
            print c,l,mod,mu,mu*mod,I_layer,m_layer,d_layer
            if f==1:
                modc_layer+=mu*mod
            else:   
                modc_layer+=mod
                x1[c]+=mod
        modc_couple=0
        for co in couple:
            d_couple_top=0
            d_couple_bot=0
            I_couple=0
            m_couple=0
            
            top_tot=0
            top_con=0
            bot_tot=0
            bot_con=0
            
            for n in couple[co]:
                if n in node_c:
                    for nei in node_c[n]:
                        if (n in layer[top[co]] and nei in layer[bot[co]]) or (n in layer[bot[co]] and nei in layer[top[co]]):
                            m_couple+=1
                
                if n in layer[top[co]]: #n belongs to the top layer of coupling
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
                                
                                            
                if n in layer[bot[co]]: #n belongs to the bottom layer of coupling
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
                                
            I_couple=float(I_couple)
            m_couple=float(m_couple)/2.0
            if I_couple>-1 and m_couple>0:
                if (bot_tot*top_tot) >0:
                    mod=((I_couple/m_couple)-((d_couple_top*d_couple_bot)/((m_couple)*(m_couple))))*(float(bot_con*top_con)/float(bot_tot*top_tot))
                else:   
                    mod=((I_couple/m_couple)-((d_couple_top*d_couple_bot)/((m_couple)*(m_couple))))
            else:
                mod=0   
            #print I_couple, m_couple, d_couple_top, d_couple_bot   
            print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
            if f==1:
                modc_couple+=2*(1-mu)*mod
            else:
                modc_couple+=mod
                x2[c]+=mod
            print modc_couple           
            #print "ha hh"
        modularity+=modc_layer+modc_couple                  
    print x1,x2 
    return 0.333*modularity 

 
    
def getModularityQ(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c):
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

def getModularity_adapt(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c):
    #edge_l
    #print edge_l,edge_c
    edge_l={}
    for l in layer:
        s=0.0
        for n in layer[l]:
            if n in node_l: 
                s+=len(node_l[n])
        edge_l[l]=s/2.0 
                
    #edge_c
    edge_c={}
    for c in couple:
        s=0.0
        for n in layer[top[c]]:
            if n in node_c: 
                s+=len(node_c[n])
        for n in layer[bot[c]]:
            if n in node_c: 
                s+=len(node_c[n])       
        edge_c[c]=s/2.0
    
    #print edge_l,edge_c
    modularity=0
    for c in commu:
        for n1 in commu[c]:
            for n2 in commu[c]:
                if not n1==n2:
                    for l in layer:
                        if n1 in layer[l] and n2 in layer[l]:
                            aij=0
                            d1=0
                            d2=0
                            if n1 in node_l:
                                d1=len(node_l[n1])
                                if n2 in node_l[n1]:
                                    aij=1
                            if n2 in node_l:
                                d2=len(node_l[n2])
                            modularity+=(aij/(2*(edge_l[l])))-((d1*d2)/((2*(edge_l[l]))*(2*(edge_l[l]))))
                    for co in couple:
                        if (n1 in layer[top[co]] and n2 in layer[bot[co]]):# or (n2 in layer[top[co]] and n1 in layer[bot[co]]):
                            aij=0
                            d1=0
                            d2=0
                            if n1 in node_c:
                                d1=len(node_c[n1])
                                if n2 in node_c[n1]:
                                    aij=1
                            if n2 in node_c:
                                d2=len(node_c[n2])
                            modularity+=(aij/((edge_c[co])))-((d1*d2)/((edge_c[co])*(edge_c[co])))
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

def build_network(layer, node_l, node_c, top, bot, couple, edge_l, edge_c):
    G = nx.Graph()
    for n in node_l:
        for n2 in node_l[n]:
            G.add_edge(n, n2)
    for n in node_c:
        for n2 in node_c[n]:
            G.add_edge(n, n2)
    return G

def remove_edge(G, node_l, node_c, commu):
    ori = nx.number_connected_components(G)
    n = ori
    sl_c = 0
    sl_nc = 0
    ml_c = 0
    ml_nc = 0
    while n <= ori:
        edge_vals = nx.edge_betweenness_centrality(G, weight='weight')
        max_ = max(edge_vals.values())
        # print max_
        # print ".",
        for k, v in edge_vals.iteritems():
            if float(v) == max_:
                G.remove_edge(k[0],k[1])
                is_ml = is_multi_layer(k[0], k[1], node_c)
                is_com = is_commu(k[0], k[1], commu)
                if is_ml:
                    if is_com:
                        ml_c+=1
                    else:
                        ml_nc +=1
                else:
                    if is_com:
                        sl_c+=1
                    else:
                        sl_nc+=1

        n = nx.number_connected_components(G)
    # print sl_c, sl_nc, ml_c, ml_nc
    return sl_c, sl_nc, ml_c, ml_nc

def getModularity(G, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu):
    communities = nx.connected_components(G)
    commu={}
    com_ID=1
    for c in communities:
        commu[com_ID]=set()
        for n in c:
            commu[com_ID].add(int(n))
        com_ID+=1
    # print commu
    # return getModularity_adapt(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c)
    return getModularityQ_modified(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c)

def newmann_girvain(G, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu):
    modu_list = []
    commu_list = []
    edges_p_list = []
    maxQ = 0.0
    Q = 0.0
    while True:
        sl_c, sl_nc, ml_c, ml_nc = remove_edge(G, node_l, node_c, commu)
        Q = getModularity(G, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu);
        # print "Modularity of decomposed G: %f" % Q
        if Q > maxQ:
            maxQ = Q
            communities = list(nx.connected_components(G))
        commu_list.append(list(nx.connected_components(G)))
        modu_list.append(Q)
        edges_p_list.append([sl_c, sl_nc, ml_c, ml_nc])
        if G.number_of_edges() == 0:
            break
    if maxQ > 0.0:
        print "Max Q: %f" % maxQ
        print "Communities:", communities
    else:
        print "Max Q: %f" % maxQ

    return commu_list, modu_list, edges_p_list, communities, maxQ


def getSeries(filename):
    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    fcommu, Qf, edge_params, bcomu, bmodu = newmann_girvain(ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu)
    print "FINAL"
    return Qf, fcommu, edge_params, bcomu, bmodu


import os
import sys

i = int(sys.argv[1])
p = float(sys.argv[2])
a = float(sys.argv[3])
mu = float(sys.argv[4])
d = float(sys.argv[5])

print p, a, mu, d
str1="./test"+str(i)+"/networks_alpha"+str(a)+"/"
str2=str1+"networks_p"+str(p)+"/networks_mu"+str(mu)+"/networks_density"+str(d)+"/new_format"
if not os.path.exists(str2) :
    print p, a, mu, d
    print "No"
    exit()
modu, commu, edges_p, bcomu, bmodu = getSeries(str2) 
with open(str2+'_commu_benching_all_modi2.pickle', 'wb') as handle:
    pickle.dump(commu, handle)
with open(str2+'_modu_benching_all_modi2.pickle', 'wb') as handle:
    pickle.dump(modu, handle)
with open(str2+'_commu_benching_frac_modi2.pickle', 'wb') as handle:
    pickle.dump(bcomu, handle)
with open(str2+'_modu_benching_frac_modi2.pickle', 'wb') as handle:
    pickle.dump(bmodu, handle)
with open(str2+'_edge_benching_p_all_modi2.pickle', 'wb') as handle:
    pickle.dump(edges_p, handle)
print ""
