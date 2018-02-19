import os, sys, time
import subprocess
import networkx as nx
import random
import pickle
import math
from sklearn.metrics import *
from collections import defaultdict
import matplotlib.pyplot as plt
from multiprocessing import Pool

N=100
num_layers =2
filepath = "./netsForDtDmDb/_networks/nets/"
file  ="network_0.7_0.7_0.05_0.5_0.0_3"

def build_network(layer, node_l, node_c, top, bot, couple, edge_l, edge_c):
	G = nx.Graph()
	for n in node_l:
		for n2 in node_l[n]:
			G.add_edge(n, n2)
	for n in node_c:
		for n2 in node_c[n]:
			G.add_edge(n, n2)
	return G


def read_raw_network(filename):
    #print(filename)
    fp=open(filename,'r')
    line=fp.readline()
    line=line.rstrip()
    n_layer=int(float(line))
    layer={}
    node_l={}
    l_ID=1
    edge_l={}
    edge_c={}
    # f_el = open(filename+'_edges_list_commod'+str(g), 'w')
    for i in range(0,n_layer):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        layer[l_ID]=set()
        ##print line
        for n in line:
            layer[l_ID].add(int(float(n)))
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        ##print n_edge
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
            # f_el.write(str(n1-1)+' '+str(n2-1)+'\n')

        l_ID+=1
        
    line=fp.readline()
    line=line.rstrip()
    n_couple=int(float(line))
    ##print n_couple
    node_c={}      
    top={}
    bot={}
    c_ID=1
    couple={}

    for i in range(0,n_couple):
        line=fp.readline()
        ##print line
        line=line.rstrip()
        line=line.split()
        top[c_ID]=int(float(line[0]))
        bot[c_ID]=int(float(line[1]))
        
        couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
        
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        ##print n_edge
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
            # f_el.write(str(n1-1)+' '+str(n2-1)+'\n')
        edge_c[c_ID] = count_edge
        c_ID=c_ID+1

    line=fp.readline()
    line=line.rstrip()
    ##print line
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
    mu=0.05

    ml_network =build_network(layer, node_l, node_c, top, bot, couple, edge_l, edge_c)

    return ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu

def writenet(dt,db, dm , ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):
	filename = filepath +"network_0.7_0.7_0.05_0.5_0.0_"+ str(float(dt)) +"_" + str(float(db))
	fp = open(filename,'w')

	fp.write(str(num_layers)+"\n")

	#write all nodes in each layer
	for l in layer:
		towrite =""
		for n in layer[l]:
			towrite= towrite+ str(n)+" "
		fp.write(towrite+'\n')

		fp.write(str(int(2*E[l]))+'\n')
		count=0
		for n1 in layer[l]:
			for n2 in node_l[n1]:
				count+=1
				fp.write(str(n1) + ' ' + str(n2) + '\n')
		#print("count=" ,count)
	fp.write("1\n")
	fp.write("1 2\n")
	fp.write(str(int(2*E12))+'\n')
	for n1 in node_c:
		for n2 in node_c[n1]:
			fp.write(str(n1) + ' ' + str(n2) + '\n')

	fp.write(str(len(commu))+'\n')
	for c in commu:
		towrite =""
		for n in commu[c]:
			towrite= towrite+ str(n)+" "
		fp.write(towrite+' \n')

	fp.close()

def applydt(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):
	edgestoadd = ((N*dt)/2.0) - E[1]
	edgesoutofcommu  =0.05*edgestoadd
	edgesincommu = edgestoadd - edgesoutofcommu
	print("edgestoadd: {0} , edgesincommu: {1} ,edgesoutofcommu: {2}".format(edgestoadd,edgesincommu,edgesoutofcommu))
	print("E[1]: ",E[1])
	loop = edgesincommu
	for n1 in layer[1]:
		if(loop<=0):
			break
		for n2 in layer[1]:
			if(loop<=0):
				break
			if(nodecomm[n1]!=nodecomm[n2] or n1==n2):
				continue
			if(n1 in node_l[n2] or n2 in node_l[n1]):
				continue
			#print("adding {0} and {1}".format(n1,n2))
			node_l[n1].add(n2)
			node_l[n2].add(n1)
			loop=loop-1
	added = edgesincommu-loop
	E[1]+=added
	print("added in commu: ",added," E[1]= ",E[1])

	loop = edgesoutofcommu
	for n1 in layer[1]:
		if(loop<=0):
			break
		for n2 in layer[1]:
			if(loop<=0):
				break
			if(nodecomm[n1]==nodecomm[n2] or n1==n2):
				continue
			if(n1 in node_l[n2]):
				continue
			node_l[n1].add(n2)
			node_l[n2].add(n1)
			loop=loop-1
	added = edgesoutofcommu-loop
	E[1]+= added
	print("added out of commu: ",added," E[1]= ",E[1])

	dt_achieved = 2.0*E[1]/100.0
	print("for " +str(dt)+" dt achieved = "+str(dt_achieved)+" done")

	return node_l,E

def applydb(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):
	edgestoadd = ((N*db)/2.0) - E[2]
	edgesoutofcommu  =0.05*edgestoadd
	edgesincommu = edgestoadd - edgesoutofcommu
	print("edgestoadd: {0} , edgesincommu: {1} ,edgesoutofcommu: {2}".format(edgestoadd,edgesincommu,edgesoutofcommu))
	print("E[2]: ",E[2])
	
	loop = edgesincommu
	for n1 in layer[2]:
		if(loop<=0):
			break
		for n2 in layer[2]:
			if(loop<=0):
				break
			if(nodecomm[n1]!=nodecomm[n2] or n1==n2):
				continue
			if(n1 in node_l[n2] or n2 in node_l[n1]):
				continue
			#print("adding {0} and {1}".format(n1,n2))
			node_l[n1].add(n2)
			node_l[n2].add(n1)
			loop=loop-1
	added = edgesincommu-loop
	E[2]+= added
	print("added in commu: ",added," E[2]= ",E[2])


	loop = edgesoutofcommu
	for n1 in layer[2]:
		if(loop<=0):
			break
		for n2 in layer[2]:
			if(loop<=0):
				break
			if(nodecomm[n1]==nodecomm[n2] or n1==n2):
				continue
			if(n1 in node_l[n2]):
				continue
			node_l[n1].add(n2)
			node_l[n2].add(n1)
			loop=loop-1
	added = edgesoutofcommu-loop
	E[2]+=added
	print("added out of commu: ",added," E[2]= ",E[2])

	dt_achieved = 2.0*E[2]/100.0
	print("for " +str(db)+" db achieved = "+str(dt_achieved)+" done")

	return node_l,E


def create(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):
	node_l,E = applydt(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)
	node_l,E = applydb(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)
	
	##print("for " +str(dt)+" dt achieved = "+str(dt_achieved)+" done")
	print(dt,db)
	writenet(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)
	return node_l,E


def createNet(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):
	edgestoadd = ((N*dt)/2.0) - E[1]
	edgesoutofcommu = 0.05*edgestoadd   #so that mu remains same after addition of edges
	edgeswithincommu  =edgestoadd - edgesoutofcommu
	loop = edgeswithincommu
	
	#Adding "edgeswithincommu" number of edges within the communities in both the layers
	#for l in layer:
	loop = edgeswithincommu
	while(loop>0):
		for n1 in layer[l]:
			if(loop<=0):
				break
			for n2 in layer[l]:
				if(loop<=0):
					break
				if(n1==n2):
					continue
				if(nodecomm[n1] != nodecomm[n2]):
					continue
				if(n1 in node_l[n2]):
					continue
				node_l[n1].add(n2)
				node_l[n2].add(n1)
				loop = loop-1
		added = edgeswithincommu - loop
		break
	E[l]+=added

	#Adding "edgesoutofcommu" number of edges out of the communities in both the layers
	#for l in layer:
	loop = edgesoutofcommu
	while(loop>0):
		for n1 in layer[l]:
			if(loop<=0):
				break
			for n2 in layer[l]:
				if(loop<=0):
					break
				if(n1==n2):
					continue
				
				if(nodecomm[n1] == nodecomm[n2]):
					continue
				if(n1 in node_l[n2]):
					continue
				node_l[n1].add(n2)
				node_l[n2].add(n1)
				loop = loop-1
		
		added = edgesoutofcommu -loop
		break
	E[l]+=added

	#ml_network is not modified , need to add these new edges to the network graph
	dt_achieved = 2.0*E[1]/100.0
	#print("for " +str(dt)+" dt achieved = "+str(dt_achieved)+" done")
	writenet(dt,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)

	return ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12


def getinfo(graph,layer,commu,node_l, node_c):
	#Construct a dict = {node1: layer, node2: layer......}
    nodelayer = {}
    nodecomm ={}
    for l in layer:
        for node in layer[l]:
            nodelayer[node] = l

    #calculate Intra_inter ----------------------------------------------------------
    intra_inter={}
    for c in commu:
        intra_inter[c]=set()
        
        for n in commu[c]:
            for l in layer:
                if n in layer[l]:
                    intra_inter[c].add(l)

    #--------------------------------------------------------------------------------

    #calculate |E1| , |E2| , |E12|---------------------------------------------------
    E={}
    E12=0
    for l in layer:
        E[l]=0
        for n in layer[l]:
            for nei in node_l.get(n,set()):
                E[l]+= 1
        E[l] = E[l]/2

    for n in node_c:
        for nei in node_c[n]:
            E12+=1
    E12 = E12/2
    #--------------------------------------------------------------------------------

    for c in commu:
    	for n in commu[c]:
    		nodecomm[n] = c

    return nodecomm , nodelayer, intra_inter,E, E12

def getSeries(filename):
	ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu =read_raw_network(filename)
	nodecomm,nodelayer, intra_inter,E, E12 = getinfo(ml_network, layer,commu, node_l, node_c)
	initial_dm= int(round(E12*1.0/N))
	
	initial_dt = int(round(2.0*E[1]/N))
	initial_db = int(round(2.0*E[2]/N))
	print(initial_db,initial_dt)

	node_lcopy = node_l
	Ecopy= E

	for dt in range(initial_dt+1,101):
		node_l = node_lcopy
		E = Ecopy
		for db in range(initial_db+1,101):
			#ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12 = createNetapplydb(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)
			node_l, E = create(dt,db,initial_dm,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)
			

getSeries("./netsForDtDmDb/_networks/" + file)
