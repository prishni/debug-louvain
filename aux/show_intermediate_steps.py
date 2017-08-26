from collections import defaultdict

def detectedcomms(p):
	print("printing fron detectedcomms\n")
	#for p in plist:
	louvain_p = defaultdict(list)
	for l in p.keys():
		louvain_p[p[l]].append(l)
	for comms, nodelist in louvain_p.items():
		print(comms ,":", nodelist)
	return louvain_p

def printgraph(ngraph): 
	print(type(ngraph), ngraph)
	print "Printing graph: "
	print(list(ngraph.nodes()))
	print(list(ngraph.edges(data = True)))
	print(type(ngraph.edges(data = True)))
	print "Done printing graphs"