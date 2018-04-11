import networkx as nx
import random
import pickle
import ast
import os

def getSeries(filename):
    edges = open(filename+'/hetero.net','r').readlines()
    met = open(filename+'/.meta','r')
    met = ast.literal_eval(met.readlines()[0].strip().strip("lr = "))
    f = []
    f.append(filename +'/edge_00')
    f.append(filename +'/edge_11')
    f.append(filename +'/edge_01')
    c = 0
    for (i, m) in enumerate(met):
        f_met = open(f[i], 'w')
        for j in range(m):
            f_met.write(edges[c+j])
        c+=m
'''
plist = [0.8]
alphalist = [0.4, 0.6, 0.8, 1.0]
mulist = [0.4]
densitylist = [0.04]
'''

list_alpha = [0.2,0.4,0.6,0.8]
list_mu = [0.05,0.1,0.2,0.3]
list_p = [0.2,0.4,0.6,0.8]
list_p1 = [0.2,0.4,0.6,0.8]
list_p2 = [0.0,0.1,0.2,0.3]
new_dir = "../syntheticNetworkGeneration/netsForcomparingBaseline/"

for alpha in list_alpha:
    new_dir_alpha = new_dir + "/alpha-" + str(alpha)

    for p in list_p:
        new_dir_p = new_dir_alpha + "/p-" + str(p)
        
        for mu in list_mu:
            new_dir_mu = new_dir_p + "/mu-" + str(mu)
        
            for p1 in list_p1:
                new_dir_p1 = new_dir_mu + "/p1-" + str(p1)
                
                for p2 in list_p2:
                    new_dir_p2 = new_dir_p1 + "/p2-" + str(p2)
                    getSeries(new_dir_p2)
'''                   
plist = [0.1, 0.25, 0.4, 0.6, 0.8]
alphalist = [0.2, 0.4, 0.6, 0.8, 1.0]
mulist = [0.05, 0.2, 0.4, 0.55, 0.75]
densitylist = [0.004, 0.01, 0.025, 0.04, 0.055, 0.07]

for i in range(11,21):
    c=0
    for p in plist:
        for a in alphalist:
            for mu in mulist:
                for d in densitylist:
                    # print p, a, mu, d
                    str1="./test"+str(i)+"/networks_alpha"+str(a)+"/"
                    str2=str1+"networks_p"+str(p)+"/networks_mu"+str(mu)+"/networks_density"+str(d)+"/new_format"
                    if not os.path.exists(str2) :
                        # print p, a, mu, d
                        # print "No"
                        continue
                    getSeries(str2,i, p, a, mu, d)
                    c+=1
    print c
'''
