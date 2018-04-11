from sys import argv
from common import *

def div_to_com(fin, fout):
    d = read_div(fin)
    savenet(d, fout)

def read_div(fin):
    with open(fin) as f:
        E = [int(line) for line in f]
    return divlist_to_cmu(E)

def divlist_to_cmu(E):
    maxCID = max(E)
    d = [[] for x in range(maxCID+1)]
    for i, x in enumerate(E):
        d[x].append(i)
    return d
    

if __name__ == '__main__':
    if len(argv) != 3:
        print('usage: divtocmu.py fin fout')
    else:
        temp, fin, fout = argv
        div_to_com(fin, fout)
    
    
