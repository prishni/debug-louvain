from common import *
from itertools import combinations
from copy import deepcopy

a1 = [[0,1,2,3],[4,5,6]]
a2 = [[0],[1,2],[3],[4,5],[6]]

def merge(*clists):
    layers = len(clists)
    original_clist = deepcopy(flatten_list(clists))
    nc_dict = {}
    for i, c in enumerate(original_clist):
        for n in c:
            try:
                nc_dict[n].append(i)
            except:
                nc_dict[n] = [i]

    cn_dict = {}
    for k, v in nc_dict.items():
        try:
            cn_dict[tuple(v)].append(k)
        except:
            cn_dict[tuple(v)] = [k]

    merged = list(cn_dict.values())
    for c in merged:
        c.sort()
    merged.sort()
    return merged
##
##def merge(*result_list):
##    '''merge(resulta, resultb, ...) -> merged_result.'''
##    n_of_result = len(result_list)
##
##    # check results for: repeat node / consistency
##    number_of_nodes = len(flatten_list(result_list[0]))
##    set_of_nodes = set(range(number_of_nodes))
##    for r in result_list:
##        fr = flatten_list(r)
##        assert len(fr) == len(set(fr)) # no repeated node in r
##        assert set(fr) == set_of_nodes # all r has same nodes and starts from 0
##
##    # sort c in r
##        result_list = [[sorted(c) for c in r] for r in result_list]
##
##
##    total = {}
##    for i, r in enumerate(result_list):
##        for c in r:
##            for lnk in combinations(c, 2):
##                try:
##                    total[lnk] += 1
##                except:
##                    total[lnk] = 1
##
##    overlapping_lnks = [k for k, v in total.items() if v == n_of_result]
##
##
##    # because all connected components are cliques
##    # therefore all nodes in a connected components are connecting to the smallest nodes
##    # e[1] > e[0]
##    c_of_n = list(range(number_of_nodes))
##    for e in overlapping_lnks:
##        if c_of_n[e[1]] > e[0]:
##            c_of_n[e[1]] = e[0]
##
##    communities = {}
##    for nid, cid in enumerate(c_of_n):
##        try:
##            communities[cid].append(nid)
##        except:
##            communities[cid] = [nid]
##
##    for c in communities.values():
##        for lnk in combinations(sorted(c), 2):
##            assert lnk in overlapping_lnks
##                
##            
##    return list(communities.values())

##
##resulta = [list(range(100)), list(range(100, 200))]
##resultb = [list(range(50)), list(range(50, 100)), list(range(100, 150)), list(range(150, 200))]
##print(merge(resulta, resultb))
##
