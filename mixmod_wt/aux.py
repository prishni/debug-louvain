from collections import defaultdict

def _get_com_wise_nodes(dictionary):
	#m = max(dictionary.values())
	louvain_p = defaultdict(set)
	for l in dictionary.keys():
	    louvain_p[dictionary[l]].add(l)
	return louvain_p

def printsomeinformation(node, com_node, best_node, incr,node_l,node_c):
	pass; return;
	print("Printing Information: {0} {1} {2} {3} {4} {5}".format(node, com_node, best_node, incr,node_l,node_c))


def __modularity(commu, status, graph):
	layer = status.layer 
	node_l = status.node_l
	node_c = status.node_c

	nodelayer = {}
	for l in layer:
		for node in layer[l]:
			nodelayer[node] = l

	x1 = {}	#stores within layer modularity values
	x2 = {}	#store coupling modularity values

	#graph[n][nei].get('weight',1)
	layerdegreesum = {}
	for l in layer:
		layerdegreesum[l] = 0
		for n1 in layer[l]:
			for n2 in node_l.get(n1, set()):
				layerdegreesum[l] += graph[n1][n2].get('weight', 1)

	couplingdegreesum = 0

	for n1 in node_c:
		for n2 in node_c[n1]:
			   	couplingdegreesum += graph[n1][n2].get('weight',1)

	#Compute modularity value for each community separately
	modularity = 0
	for c in commu:
		x1[c] = 0
		x2[c] = 0
		if(len(commu[c])==1): modl = 0 
		else:
			for l in layer:
				modl = 0
				for n1 in commu[c]:
					d1 = 0
					for n3 in node_l[n1]:
						d1 += graph[n1][n3].get('weight',1)

					for n2 in commu[c]:
						aij = 0
						d2 = 0
						if(nodelayer[n2]!=nodelayer[n1]): continue

						if((n1==n2 and n2 in node_l.get(n1, set())) or n1!=n2):
						    if(n2 in node_l[n1]):
						        if(n1!=n2): aij = graph[n1][n2].get('weight',1)
						        else: aij = graph[n1][n2].get('weight',1)

						    for n3 in node_l[n2]:
						        d2 += graph[n2][n3].get('weight',1)

						modl += (aij*1.0-(d1*d2*1.0)/layerdegreesum[l])/layerdegreesum[l]
			x1[c] += modl*1.0/2


		modc = 0
		for n1 in commu[c]:
			d1 = 0
			for n3 in node_c.get(n1, set()):
				d1 += graph[n1][n3].get('weight',1)
			
			for n2 in commu[c]:
				if(n1==n2): continue
				if(n2 not in node_c.get(n1, set())): continue

				d2 = 0
				for n3 in node_c.get(n2, set()):
					d2 += graph[n2][n3].get('weight',1)

				aij = graph[n1][n2].get('weight',1)
					
				if(couplingdegreesum == 0): 
					modc = 0
				else:
					modc += (aij*1.0-((d1*d2*1.0)/couplingdegreesum))/couplingdegreesum

		x2[c] = modc/2;
		modularity += x1[c] + x2[c]

	return (1.0/3)*modularity

