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
from copy import deepcopy

N=100
num_layers =2
mu = 0.05
#filepath = "./netsForDtDmDb/_networks/morenets/alpha0.7/"
#file  ="network_0.7_0.7_0.05_0.7_0.0_3"
#file  ="network_0.3_0.7_0.05_0.7_0.0_10"

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
    ##print(filename)
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
        ###print line
        for n in line:
            layer[l_ID].add(int(float(n)))
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        ###print n_edge
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
    ###print n_couple
    node_c={}      
    top={}
    bot={}
    c_ID=1
    couple={}

    for i in range(0,n_couple):
        line=fp.readline()
        ###print line
        line=line.rstrip()
        line=line.split()
        top[c_ID]=int(float(line[0]))
        bot[c_ID]=int(float(line[1]))
        
        couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
        
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        ###print n_edge
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
    ###print line
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

def writenet(dt,db, dm , ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12,filepath):
	filename = filepath + str(float(dm)) +"_"+ str(float(dt)) +"_" + str(float(db))
	fp = open(filename,'w')

	fp.write(str(num_layers)+"\n")

	#write all nodes in each layer
	for l in layer:
		towrite =""
		for n in layer[l]:
			towrite= towrite+ str(n)+" "
		fp.write(towrite+'\n')

		fp.write(str(int(2*E[l]))+'\n')
		#print(2*E[l])
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

#def create(dt,db,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):


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

def applyremoval(d,layer,node_l,commu,nodelayer,nodecomm,E,l):
	toremove = E[l] - (N*d/2.0)
	removeoutofcommu = int(round(mu*toremove))
	removewithincommu  =int(round(toremove - removeoutofcommu))
	#print("toremove: {0} , removewithincommu: {1} ,removeoutofcommu: {2}".format(toremove,removewithincommu,removeoutofcommu))
	#print("E[1]: ",E[l])

	loop = removewithincommu
	while(loop>0):
		modif=0
		for n1 in layer[l]:
			if(loop<=0):
				break;
			for n2 in layer[l]:
				if(loop<=0):
					break;
				if(nodecomm[n1]!=nodecomm[n2] or n1==n2):
					continue
				if(n1 not in node_l[n2] or n2 not in node_l[n1]):
					continue
				##print("adding {0} and {1}".format(n1,n2))
				node_l[n1].remove(n2)
				node_l[n2].remove(n1)
				modif =1
				loop=loop-1
		if(modif==0 or loop<=0):
			removed = removewithincommu-loop
			E[l]-=removed
			#print("removed in commu: ",removed," E["+str(l)+"]= ",E[l])
			break

	loop = removeoutofcommu
	while(loop>0):
		modif=0
		for n1 in layer[l]:
			if(loop<=0):
				break;
			for n2 in layer[l]:
				if(loop<=0):
					break;
				if(nodecomm[n1]==nodecomm[n2] or n1==n2):
					continue
				if(n1 not in node_l[n2] or n2 not in node_l[n1]):
					continue
				##print("adding {0} and {1}".format(n1,n2))
				node_l[n1].remove(n2)
				node_l[n2].remove(n1)
				modif =1
				loop=loop-1
		if(modif==0 or loop<=0):
			removed = removeoutofcommu-loop
			E[l]-=removed
			#print("removed out of commu: ",removed," E["+str(l)+"]= ",E[l])
			break
	#print("for layer",l," density achoeved: ",(2.0*E[l]/N)," desired: ",d)
	return E,node_l

def applyadditional(d,layer,node_l,commu,nodelayer,nodecomm,E,l):
	edgestoadd = ((N*d)/2.0) - E[l]
	edgesoutofcommu  =int(round(mu*edgestoadd))
	edgesincommu = int(round(edgestoadd - edgesoutofcommu))
	#print("edgestoadd: {0} , edgesincommu: {1} ,edgesoutofcommu: {2}".format(edgestoadd,edgesincommu,edgesoutofcommu))
	#print("E[1]: ",E[l])

	loop = edgesincommu
	while(loop>0):
		modif=0
		for n1 in layer[l]:
			if(loop<=0):
				break
			for n2 in layer[l]:
				if(loop<=0):
					break
				if(nodecomm[n1]!=nodecomm[n2] or n1==n2):
					continue
				if(n1 in node_l[n2] or n2 in node_l[n1]):
					continue
				##print("adding {0} and {1}".format(n1,n2))
				node_l[n1].add(n2)
				node_l[n2].add(n1)
				modif=1
				loop=loop-1
		if(modif==0 or loop<=0):
			added = edgesincommu-loop
			E[l]+=added
			#print("added in commu: ",added," E["+str(l)+"]= ",E[l])
			break

	loop = edgesoutofcommu
	while(loop>0):
		modif=0
		for n1 in layer[l]:
			if(loop<=0):
				break
			for n2 in layer[l]:
				if(loop<=0):
					break
				if(nodecomm[n1]==nodecomm[n2] or n1==n2):
					continue
				if(n1 in node_l[n2] or n2 in node_l[n1]):
					continue
				##print("adding {0} and {1}".format(n1,n2))
				node_l[n1].add(n2)
				node_l[n2].add(n1)
				modif=1
				loop=loop-1
		if(modif==0 or loop<=0):
			added = edgesoutofcommu-loop
			E[l]+=added
			#print("added out of commu: ",added," E["+str(l)+"]= ",E[l])
			break
	#print("for layer",l," density achoeved: ",(2.0*E[l]/N)," desired: ",d)
	return E,node_l


def getSeries(filename,newfilepath):
	ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu =read_raw_network(filename)
	nodecomm,nodelayer, intra_inter,E, E12 = getinfo(ml_network, layer,commu, node_l, node_c)
	
	#initial_dm= int(round(E12*1.0/N))
	initial_dm= (E12*1.0/N)
	initial_dt = int(round(2.0*E[1]/N))
	initial_db = int(round(2.0*E[2]/N))
	#print(initial_db,initial_dt,initial_dm)
	required_fractions= [0.1,0.3,0.5,0.7,0.9,1.0,3.0,5.0,7.0,9.0,10.0]
	required_dts = [v*initial_dm for  v in reversed(required_fractions)]
	
	print(len(set(required_dts)))
	
	node_lcopy = deepcopy(node_l)
	Ecopy= deepcopy(E)

	for dt in required_dts:
		#dt=5.0
		node_l = deepcopy(node_lcopy)
		E = deepcopy(Ecopy)

		if(dt < initial_dt):
			E,node_l = applyremoval(dt,layer,node_l,commu,nodelayer,nodecomm,E,1)
		else:
			E,node_l = applyadditional(dt,layer,node_l,commu,nodelayer,nodecomm,E,1)
		
		node_lcopydt = deepcopy(node_l)
		Ecopydt= deepcopy(E)
		for db in required_dts:
			#db=10.0
			#print("for bd: ",db," dt:",dt)
			node_l = deepcopy(node_lcopydt)
			E = deepcopy(Ecopydt)
			if(db < initial_db):
				E,node_l = applyremoval(db,layer,node_l,commu,nodelayer,nodecomm,E,2)
			else:
				E,node_l = applyadditional(db,layer,node_l,commu,nodelayer,nodecomm,E,2)
			writenet(dt,db, initial_dm , ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12,newfilepath)
			#print("******************************************************************")
			#break
			#create(dt,db,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12)
		#break

#getSeries("./netsForDtDmDb/_networks/baseNetworks/" + file , filepath)

