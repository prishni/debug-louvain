from common import *
from sys import argv

def ret_c(C, CC):
    retrC = [flatten_list((C[x] for x in c)) for c in CC]
    return retrC

def rec_ret_c(clists):
    '''recursively retrieve clist into one c'''
    if len(clists) == 2:
        return ret_c(clists[0], clists[1])
    else:
        return ret_c(clists[0], rec_ret_c(clists[1:]))


def treec_to_cmu(E):
    '''same_func as div_to_cmu'''
    maxCID = max(E)
    d = [[] for x in range(maxCID+1)]
    for i, x in enumerate(E):
        d[x].append(i)
    return d

def tree_to_cmu(tree):
    splitpoints = [i for i, x in enumerate(tree) if x[0] == 0]
    assert(len(splitpoints) >= 1)

    tree_c = [x[1] for x in tree] # only keep c information
    nc_lists = [tree_c[start:end] for start, end in zip(splitpoints, splitpoints[1:]+[len(tree)])]

    clists = list(map(treec_to_cmu, nc_lists))
    cmu = rec_ret_c(clists)
    for c in cmu:
        c.sort()
    return cmu

def main(tree_fn, cmu_fn):
    tree = readnet(tree_fn)
    cmu = tree_to_cmu(tree)
    savenet(cmu, cmu_fn)

if __name__ == '__main__':
    if len(argv) != 2:
        scriptname = argv[0].split('\\')[-1]
        print('usage: {} tree_fn cmu_fn'.format(scriptname))
        main(tree_fn, cmu_fn)
