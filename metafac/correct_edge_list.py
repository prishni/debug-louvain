import networkx as nx
import random
import pickle
import os

def read_raw_network(filename,weighted):
    print(filename)
    fp=open(filename,'r')
    line=fp.readline()
    line=line.rstrip()
    n_layer=int(float(line))
    layer={}
    node_l={}
    l_ID=1
    edge_l={}
    edge_c={}
    weightednode_l ={}
    weightednode_c ={}
    # f_el = open(filename+'_edges_list_commod'+str(g), 'w')
    for i in range(0,n_layer):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        layer[l_ID]=set()
        #print line
        for n in line:
            layer[l_ID].add(int(float(n)))
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        #print n_edge
        edge_l[l_ID]=n_edge
        for j in range(0,n_edge):
            line=fp.readline()
            line=line.rstrip()
            line=line.split()
            n1=int(float(line[0]))
            n2=int(float(line[1]) )
            if n1 not in node_l:
                node_l[n1]=set()
            node_l[n1].add(n2)    
            if n2 not in node_l:
                node_l[n2]=set()
            node_l[n2].add(n1)
            
            if(weighted==1):
                weightednode_l[n1]  = weightednode_l.get(n1,dict())
                weightednode_l[n2]  = weightednode_l.get(n2,dict())
                weightednode_l[n1][n2] = int(float(line[2]))
                weightednode_l[n2][n2] = int(float(line[2]))

            
        l_ID+=1
        
    line=fp.readline()
    line=line.rstrip()
    n_couple=int(float(line))
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
        top[c_ID]=int(float(line[0]))
        bot[c_ID]=int(float(line[1]))
        
        couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
        
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        #print n_edge
        edge_c[c_ID]=n_edge
        count_edge = 0
        for j in range(0,n_edge):
            line=fp.readline()
            line=line.rstrip()
            line=line.split()
            n1=int(float(line[0]))
            n2=int(float(line[1]))
            if n1 not in node_c:
                node_c[n1]=set()
            node_c[n1].add(n2)
            if n2 not in node_c:
                node_c[n2]=set()
            node_c[n2].add(n1)  
            count_edge += 1

            if(weighted==1):
                weightednode_c[n1]  = weightednode_c.get(n1,dict())
                weightednode_c[n2]  = weightednode_c.get(n2,dict())
                weightednode_c[n1][n2] = int(float(line[2]))
                weightednode_c[n2][n2] = int(float(line[2]))

            # f_el.write(str(n1-1)+' '+str(n2-1)+'\n')
        edge_c[c_ID] = count_edge
        c_ID=c_ID+1

    line=fp.readline()
    line=line.rstrip()
    #print line
    n_comm=int(float(line))
    commu={}
    com_ID=1
    for i in range(0,n_comm):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        commu[com_ID]=set()
        for n in line:
            commu[com_ID].add(int(float(n)))
        com_ID+=1      
    mu=0

    ml_network =build_network(layer, node_l, node_c,weightednode_l,weightednode_c, top, bot, couple, edge_l, edge_c,weighted)

    return ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu,commu
    
def getSeries(filename):

    fnetwork = 0
    '''with open(filename+'_ml_network.pickle') as handle:
                    fnetwork = pickle.load(handle)'''
    #ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = fnetwork
    ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu = read_raw_network(filename,0)
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
                    new_dir_p2 = new_dir_p1 + "/p2-" + str(p2)+"/new_format"
                    getSeries(new_dir_p2)
            
'''
for i in range(11, 21):
    c=0
    for p in plist:
        for a in alphalist:
            for mu in mulist:
                for d in densitylist:
                    # print p, a, mu, d, 
                    str1="./test"+str(i)+"/networks_alpha"+str(a)+"/"
                    str2=str1+"networks_p"+str(p)+"/networks_mu"+str(mu)+"/networks_density"+str(d)+"/new_format"
                    if not os.path.exists(str2) :
                        # print p, a, mu, d
                        # print "No"
                        continue
                    getSeries(str2)
                # print ""
                    c+=1
    print c
'''

