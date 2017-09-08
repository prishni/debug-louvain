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