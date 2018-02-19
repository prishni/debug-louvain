import networkx as nx
import random
import pickle

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

def getSeries(filename):
    fp=open(filename,'r')
    line=fp.readline()
    line=line.rstrip()
    n_layer=int(line)
    layer={}
    node_l={}
    l_ID=1
    edge_l={}
    edge_c={}
    f_el = open(filename+'_edges_list_commod', 'w')
    for i in range(0,n_layer):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        layer[l_ID]=set()
        #print line
        for n in line:
            layer[l_ID].add(int(n))
        line=fp.readline()
        line=int(line.rstrip())
        n_edge=line
        #print n_edge
        edge_l[l_ID]=n_edge
        for j in range(0,n_edge):
            line=fp.readline()
            line=line.rstrip()
            line=line.split()
            n1=int(line[0])
            n2=int(line[1]) 
            if n1 not in node_l:
                node_l[n1]=set()
            node_l[n1].add(n2)      
            if n2 not in node_l:
                node_l[n2]=set()
            node_l[n2].add(n1)
            f_el.write(str(n1-1)+' '+str(n2-1)+'\n')

        l_ID+=1
        
    line=fp.readline()
    line=line.rstrip()
    n_couple=int(line)
    #print n_couple
    node_c={}       
    top={}
    bot={}
    c_ID=1
    couple={}

    for i in range(0,n_couple):
        line=fp.readline()
        #print line
        line=line.rstrip()
        line=line.split()
        top[c_ID]=int(line[0])
        bot[c_ID]=int(line[1])
        
        couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
        
        line=fp.readline()
        line=int(line.rstrip())
        n_edge=line
        #print n_edge
        edge_c[c_ID]=n_edge
        count_edge = 0
        for j in range(0,n_edge):
            line=fp.readline()
            line=line.rstrip()
            line=line.split()
            n1=int(line[0])
            n2=int(line[1])
            if n1 not in node_c:
                node_c[n1]=set()
            node_c[n1].add(n2)
            if n2 not in node_c:
                node_c[n2]=set()
            node_c[n2].add(n1)  
            count_edge += 1
            f_el.write(str(n1-1)+' '+str(n2-1)+'\n')
        edge_c[c_ID] = count_edge
        c_ID=c_ID+1

    line=fp.readline()
    line=line.rstrip()
    #print line
    n_comm=int(line)
    commu={}
    com_ID=1
    for i in range(0,n_comm):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        commu[com_ID]=set()
        for n in line:
            commu[com_ID].add(int(n))
        com_ID+=1       

    count_c = 0

    # for cn in node_c:
    #     count_c +=len(node_c)

    # print "COUNTC"+str(count_c/2),

    # for c in commu:
    #     oset = set(commu[c])
    #     cset = set()
    #     for c2 in commu:
    #         if c!=c2:
    #             cset.update(commu[c2])

    #     for on in oset:
    #         if on in node_c :
    #             for cn in list(node_c[on]):
    #                 if cn in cset and random.random() <=g:
    #                     node_c[on].discard(cn)
    #                     node_c[cn].discard(on)
    # count_c=0
    # for cn in node_c.keys():
    #     if len(node_c[cn]) == 0:
    #         del node_c[cn]
    # for cn in node_c:
    #     count_c +=len(node_c)

    # print len(layer),len(node_l),len(node_c),top,bot,len(commu),len(couple)
    # print edge_c, edge_l, count_c/2
    slc=0
    slnc=0
    mlc=0
    mlnc=0

    for n in node_l:
        for k in node_l[n]:
            if is_commu(k, n, commu):
                slc+=1
            else :
                slnc+=1
    for n in node_c:
        for k in node_c[n]:
            if is_commu(k, n, commu):
                mlc+=1
            else :
                mlnc+=1
    # print slc, slnc, mlc, mlnc
    # print len(commu)
    # cc = {}
    # for i in commu:
    #     if len(commu[i]) not in cc:
    #         cc[len(commu[i])] = 0
    #     cc[len(commu[i])] += 1
    # print cc
    mu=0
    ml_network = build_network(layer, node_l, node_c, top, bot, couple, edge_l, edge_c)
    f_met = open(filename+'_meta_commod', 'w')
    f_met.write("lr = ["+str(edge_l[1])+", "+str(edge_l[2])+", "+str(edge_c[1])+"]\n")
    f_met.write("nr = [100, 100]")

    with open(filename+'_ml_network.pickle', 'wb') as handle:
        pickle.dump([ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu], handle)

import os
#alpha,p,mu, p1,p2
list_alpha = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
list_mu = [0.05, 0.20, 0.4, 0.55, 0.6, 0.75]
list_p = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
list_p1 = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
list_p2 = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
'''
plist = [0.1, 0.25, 0.4, 0.6, 0.8]
alphalist = [0.2, 0.4, 0.6, 0.8, 1.0]
mulist = [0.05, 0.2, 0.4, 0.55, 0.75]
densitylist = [0.004, 0.01, 0.025, 0.04, 0.055, 0.07]
'''
count =0
count2=0

import sys
from os import listdir
from os.path import isfile, join

mypath = sys.argv[1]

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for f in onlyfiles:
    if "_edges_list_commod" not in f and "_meta_commod" not in f and "_ml_network.pickle" not in f:
        print f
        getSeries(mypath + '/' + f)

# for a in list_alpha:
#     count=0
#     count2=0
#     for p in list_p:
#         for mu in list_mu:
#             for p1 in list_p1:
#                 for p2 in list_p2:
#                     str1="./n100_net_new/"
#                     str2=str1+"network_"+str(a)+"_"+str(p)+"_"+str(mu)+"_"+str(p1)+"_"+str(p2)
#                     #print str2
#                     count2+=1
#                     if not os.path.exists(str2) :
#                         # print p, a, mu, d
#                         # print "No"
#                         count+=1
#                         continue
#                     getSeries(str2)

#     print count,count2
