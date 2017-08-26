import networkx as nx
from community_louvain import best_partition

g = nx.Graph()
g.add_nodes_from(range(0, 17))



edgestr = '''
[(0, 0, 47), (0, 1, 2), (0, 5, 1), (0, 7, 1), (0, 13, 1), (0, 16, 1), (1, 1, 31), 
(1, 11, 1), (1, 14, 2), (2, 2, 53), (2, 11, 1), (2, 10, 1), (2, 15, 1), (3, 3, 32),
 (4, 16, 12), (4, 8, 1), (4, 10, 2), (4, 11, 1), (4, 4, 29), (5, 5, 29), (5, 8, 1), (5, 10, 12), 
 (5, 13, 2), (5, 14, 1), (6, 6, 15), (7, 7, 77), (7, 9, 1), (7, 12, 1), (7, 13, 1), (7, 14, 1), 
 (8, 8, 41), (8, 14, 1), (8, 13, 1), (9, 9, 22), (9, 12, 17), (9, 14, 1), (10, 10, 9), (11, 11, 52), 
 (11, 14, 2), (12, 12, 38), (12, 13, 1), (12, 15, 1), (13, 13, 73), (13, 16, 1), (14, 14, 49),
  (15, 15, 10), (16, 16, 17)]
'''

weightededgelist = eval(edgestr)
unweightededgelist = [(e[0], e[1]) for e in weightededgelist]

g.add_weighted_edges_from(weightededgelist)

#g.add_edges_from(unweightededgelist)

print(list(g.edges(data = True)))
communities = best_partition(g)
print(communities)