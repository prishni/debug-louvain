import os, sys, time
from louvain_mixmod_m21_wt_without_global_JL_Final0 import parallelimplementation
from louvain_OriginalModAlways import parallelimplementation as parallelimplementation_AllLouvain

print("starts..")
alphalist = [0.3,0.5,0.8]
plist =[0.3,0.5,0.8]
p1list =[0.3,0.5,0.8]

def mixmod_final0():
	weighted=0
	pathtosave = "./syntheticNetworkGeneration/results/params/nonweighted/"
	networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/nonweighted/"

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
				parallelimplementation(networklist,weighted,modfilename,readnetsfrom3)

def all_louvain():
	weighted=0
	pathtosave = "./syntheticNetworkGeneration/results/params/weighted/"
	networkpath="./syntheticNetworkGeneration/netsForDtDmDb/_networks/params/weighted/"

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
				parallelimplementation_AllLouvain(networklist,weighted,modfilename,readnetsfrom3)

#all_louvain()
mixmod_final0()


