import networkx as nx
import random
import pickle

def getSeries(filename):

    fnetwork = 0
    with open(filename+'_ml_network.pickle') as handle:
        fnetwork = pickle.load(handle)
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    c_1 = 0
    c_2 = 0
    c_12 = 0
    f_el = open(filename+'_edges_list_commod2', 'w')
    for i in node_l:
        if i <=100:
            for n in list(node_l[i]):
                node_l[n].discard(i)
                f_el.write(str(i-1)+" "+str(n-1)+"\n")
                c_1+=1
    for i in node_l:
        if i > 100:
            for n in list(node_l[i]):
                node_l[n].discard(i)
                f_el.write(str(i-1)+" "+str(n-1)+"\n")
                c_2+=1
    for i in node_c:
        for n in list(node_c[i]):
            node_c[n].discard(i)
            f_el.write(str(i-1)+" "+str(n-1)+"\n")
            c_12+=1

    # print c_1, c_2, c_12
    f_met = open(filename+'_meta_commod2', 'w')
    f_met.write("lr = ["+str(c_1)+", "+str(c_2)+", "+str(c_12)+"]\n")
    f_met.write("nr = [100, 100]")

import os
'''
plist = [0.1, 0.25, 0.4, 0.6, 0.8]
alphalist = [0.2, 0.4, 0.6, 0.8, 1.0]
mulist = [0.05, 0.2, 0.4, 0.55, 0.75]
densitylist = [0.004, 0.01, 0.025, 0.04, 0.055, 0.07]
'''
list_alpha = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
list_mu = [0.05, 0.20, 0.4, 0.55, 0.6, 0.75]
list_p = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
list_p1 = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
list_p2 = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
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
#     c=0
#     for p in list_p:
#         for mu in list_mu:
#             for p1 in list_p1:
#                 for p2 in list_p2:
#                     # print p, a, mu, d, 
#                     str1="./n100_net_new/"
#                     str2=str1+"network_"+str(a)+"_"+str(p)+"_"+str(mu)+"_"+str(p1)+"_"+str(p2)
#                     if not os.path.exists(str2) :
#                         # print p, a, mu, d
#                         # print "No"
#                         continue
#                     getSeries(str2)
#                 # print ""
#                     c+=1
#     print c

