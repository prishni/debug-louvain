###########################################################################################
#
#  Code will run louvain_originalmodAlways and mixmodFinal0 on all the nets for generating 
#  modularity/nmi heatmaps. (For both weighted and nonweighted)
#
###########################################################################################


import os, sys, time
#from louvain_mixmod_m21_wt_without_global_JL_Final0 import parallelimplementation as _mixmodFinal0
#from louvain_OriginalModAlways import parallelimplementation as _AllLouvain
from union import parallelimplementation as _mixmod_union
#from louvain_mixmod_m21_wt_without_global_JL_Final1 import parallelimplementation as _mixmodFinal1

print("starts..")
alphalist = [0.3,0.5,0.8]
plist =[0.3,0.5,0.8]
p1list =[0.3,0.5,0.8]

#alphalist = [0.8]
#plist =[0.8]
#p1list =[0.8]
'''
def mixmod_final1(weighted):
	if(weighted == 0):
		pathtosave = "./syntheticNetworkGeneration/results/params/nonweighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/nonweighted/"
	else:
		pathtosave = "./syntheticNetworkGeneration/results/params/weighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted/"
		#pathtosave = "./syntheticNetworkGeneration/results/params/weighted_dm1/"
		#networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted_dm1/"
		
	for alpha in alphalist:
		readnetsfrom1 = networkpath+ "aplha" + str(alpha) 
		directory1 = pathtosave+ "aplha" + str(alpha) 
		if not os.path.exists(directory1):
			os.makedirs(directory1)
		for p in plist:
			readnetsfrom2 = readnetsfrom1 + "/P" + str(p)
			directory2 = directory1 + "/P" + str(p)
			if not os.path.exists(directory2):
				os.makedirs(directory2)
			for p1 in p1list:
				readnetsfrom3 = readnetsfrom2 + "/P1_" + str(p1) +"/"
				directory3 = directory2 + "/P1_" + str(p1)
				if not os.path.exists(directory3):
					os.makedirs(directory3)
				networklist = os.listdir(readnetsfrom3)
				networklist=sorted(networklist)
				modfilename = directory3+"/mixmodFinal1"
				_mixmodFinal0(networklist,weighted,modfilename,readnetsfrom3)

def mixmod_final0(weighted):
	if(weighted == 0):
		pathtosave = "./syntheticNetworkGeneration/results/params/nonweighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/nonweighted/"
	else:
		pathtosave = "./syntheticNetworkGeneration/results/params/weighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted/"
		#pathtosave = "./syntheticNetworkGeneration/results/params/weighted_dm1/"
		#networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted_dm1/"
		
	for alpha in alphalist:
		readnetsfrom1 = networkpath+ "aplha" + str(alpha) 
		directory1 = pathtosave+ "aplha" + str(alpha) 
		if not os.path.exists(directory1):
			os.makedirs(directory1)
		for p in plist:
			readnetsfrom2 = readnetsfrom1 + "/P" + str(p)
			directory2 = directory1 + "/P" + str(p)
			if not os.path.exists(directory2):
				os.makedirs(directory2)
			for p1 in p1list:
				readnetsfrom3 = readnetsfrom2 + "/P1_" + str(p1) +"/"
				directory3 = directory2 + "/P1_" + str(p1)
				if not os.path.exists(directory3):
					os.makedirs(directory3)
				networklist = os.listdir(readnetsfrom3)
				networklist=sorted(networklist)
				modfilename = directory3+"/mixmod_for_corrImpl"
				_mixmodFinal0(networklist,weighted,modfilename,readnetsfrom3)

def all_louvain(weighted):
	if(weighted ==0):
		pathtosave = "./syntheticNetworkGeneration/results/params/nonweighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/nonweighted/"
	else : 
		pathtosave = "./syntheticNetworkGeneration/results/params/weighted/"
		#pathtosave = "./syntheticNetworkGeneration/results/params/weighted_dm1/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted/"
		#networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted_dm1/"

	for alpha in alphalist:
		readnetsfrom1 = networkpath+ "aplha" + str(alpha) 
		directory1 = pathtosave+ "aplha" + str(alpha) 
		if not os.path.exists(directory1):
			os.makedirs(directory1)
		for p in plist:
			readnetsfrom2 = readnetsfrom1 + "/P" + str(p)
			directory2 = directory1 + "/P" + str(p)
			if not os.path.exists(directory2):
				os.makedirs(directory2)
			for p1 in p1list:
				readnetsfrom3 = readnetsfrom2 + "/P1_" + str(p1) +"/"
				directory3 = directory2 + "/P1_" + str(p1)
				if not os.path.exists(directory3):
					os.makedirs(directory3)
				networklist = os.listdir(readnetsfrom3)
				networklist=sorted(networklist)
				modfilename = directory3+"/all_louvain_for_corrImpl"
				_AllLouvain(networklist,weighted,modfilename,readnetsfrom3)

'''
def mixmod_union(weighted):
	if(weighted == 0):
		pathtosave = "./syntheticNetworkGeneration/results/params/nonweighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/nonweighted/"
		
	else:
		pathtosave = "./syntheticNetworkGeneration/results/params/weighted/"
		networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted/"
		#pathtosave = "./syntheticNetworkGeneration/results/params/weighted_dm1/"
		#networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted_dm1/"
		
	for alpha in alphalist:
		readnetsfrom1 = networkpath+ "aplha" + str(alpha) 
		directory1 = pathtosave+ "aplha" + str(alpha) 

		if not os.path.exists(directory1):
			os.makedirs(directory1)
		for p in plist:
			readnetsfrom2 = readnetsfrom1 + "/P" + str(p)
			directory2 = directory1 + "/P" + str(p)
			if not os.path.exists(directory2):
				os.makedirs(directory2)
			for p1 in p1list:
				readnetsfrom3 = readnetsfrom2 + "/P1_" + str(p1) +"/"
				directory3 = directory2 + "/P1_" + str(p1)
				if not os.path.exists(directory3):
					os.makedirs(directory3)
				networklist = os.listdir(readnetsfrom3)
				networklist=sorted(networklist)
				modfilename = directory3+"/mixmod_union"
				detectedCommuFile = directory3+"/communityDetected"
				_mixmod_union(networklist,weighted,modfilename,readnetsfrom3,detectedCommuFile)

weighted = 0
#all_louvain(weighted)
#mixmod_final0(weighted)
mixmod_union(weighted)
#mixmod_final1(weighted)

weighted  =1
#all_louvain(weighted)
#mixmod_final0(weighted)
mixmod_union(weighted)
#mixmod_final1(weighted)



