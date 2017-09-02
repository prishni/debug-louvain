import pickle
import os
import sys
from collections import defaultdict
import json

def write_GT(file,fout):
	f=open(file)
	for line in f:
		if(line=='\n'):
			break
		fout.write(line)

f = open("nets/network_0.9_1.0_0.05_1.0_0.0_commu_benching_all_march21_louvain.pickle", 'rb')
p = pickle.load(f)
m = max(p.values())
louvain_p = defaultdict(list)
for l in p.keys():
	louvain_p[p[l]].append(l)

fout = open("deletethisshit.txt", 'w')
file = "./Raphael_27.6.17/infos/network_0.9_1.0_0.05_1.0_0.0.info"
write_GT(file,fout)
fout.write("\n")
for k, v in louvain_p.items():
	strtowrite = str(k) + " " + str(v) + '\n'
	fout.write(strtowrite)