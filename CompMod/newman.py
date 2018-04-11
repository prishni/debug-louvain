from itertools import chain, groupby
from collections import Counter
from random import choice, shuffle, randint
from common import *

def gennodeseq(bypass = 0, *ns):
    #multiple layer version
    nseq = list(range(sum(ns)))
    shuffle(nseq)
    pnseq = set()
        
    def nodeseq(rst = False):
        if rst: # if rst is True, just reset list and return
            nseq.extend(list(pnseq))
            pnseq.clear()
            return

        # randomly pop a num from nseq
        # add it to pnseq, then return it
        # if no nodes left, return None
        # (optional) if already |bypass| nodes passed, return None
        if bypass != 0 and len(pnseq) > bypass:
            return None
        if not nseq:
            return None
        pos = randint(0, len(nseq)-1)
        nseq[pos], nseq[-1] = nseq[-1], nseq[pos]
        rand_id = nseq.pop(-1)
        pnseq.add(rand_id)

        base_id = 0
        for layer, noderange in enumerate(ns):
            if base_id <= rand_id < base_id + noderange:
                return (layer, rand_id - base_id)
            base_id += noderange
    return nodeseq


class Node:
    def __init__(self, NID):
        self.NID = NID
        self.layer = 0
        self.degree = 0
        self.neilist = Counter()
        self.comm = None

    def __repr__(self):
        return "Node {}".format(self.NID)

    def addnei(self, neinode, weight=1):
        self.degree += weight
        self.neilist[neinode] += weight

class Community:
    def __init__(self, CID, nodes = set()):
        self.CID = CID
        self.layer = 0
        self.nodes = set()
        self.aii = 0
        self.eii = 0
        for node in nodes:
            self.nodes.add(node)
            node.comm = self

    def __len__(self):
        return len(self.nodes)

    def isempty(self):
        return len(self.nodes) == 0

    def __repr__(self):
        return "Comm {}".format(self.CID)

    
    def update_ae(self):
        self.aii = 0
        self.eii = 0

        
        for member in self.nodes:
            self.eii += member.degree

            for nei, cnt in member.neilist.items():
                if nei.comm == member.comm:
                    self.aii += cnt

                    
            
class Graph:
    def __init__(self):
        self.E = None
        self.M = 0
        self.N = 0
        self.nlist = []
        self.clist = []

    def readfile(self, filename):
        E = readnet(filename)
        self.updateE(E)

    def updateE(self, E):
        if len(E[0]) == 3:
            print('weighted E, init Newman as Weighted')
            self.updateWE(E)
        else:
            self.E = E
            self.M = len(E)
            self.N = max(flatten_list(E)) + 1

            self.nlist = [Node(i) for i in range(self.N)]
            for e in E:
                n0, n1 = self.nlist[e[0]], self.nlist[e[1]]
                n0.addnei(n1)
                n1.addnei(n0)

            self.initC()

    def updateWE(self, WE):
        self.E = None
        n0, n1, w = zip(*WE)
        self.M = sum(w)
        self.N = max(max(n0), max(n1)) + 1

        self.nlist = [Node(i) for i in range(self.N)]
        for e0, e1, w in WE:
            n0, n1 = self.nlist[e0], self.nlist[e1]
            n0.addnei(n1, w)
            n1.addnei(n0, w)

        self.initC()

    def initC(self):
        self.clist = [Community(n.NID, [n]) for n in self.nlist]

        for c in self.clist:
            c.update_ae()

    def updateC(self, cl):
        self.clist = [Community(cid, [self.nlist[n] for n in c]) for cid, c in enumerate(cl)]

        for c in self.clist:
            c.update_ae()

    def Modularity(self):
        m = self.M
        cQ_list = [c.aii/(2*m) - (c.eii/(2*m))**2 for c in self.clist]
        return sum(cQ_list)

    def move_node(self, node, dst_c):
        src_c = node.comm
        # hai shi sui shou xie de
        node.comm = dst_c
        src_c.nodes.remove(node)
        dst_c.nodes.add(node)

        for nei, cnt in node.neilist.items():
            if nei == node:
                src_c.aii -= cnt
                dst_c.aii += cnt
            else:
                if nei.comm == dst_c:
                    dst_c.aii += cnt * 2
                if nei.comm == src_c:
                    src_c.aii -= cnt * 2
                

        src_c.eii -= node.degree
        dst_c.eii += node.degree

    def mv_nd(self, nid, dst_cid):
        self.move_node(self.nlist[nid], self.clist[dst_cid])
   
    def calc_argmaxC(self, nid):
        node = self.nlist[nid]
        src_c = node.comm
        
        nvlist = node.neilist.keys()
        nclist = set([n.comm for n in nvlist])
        if src_c in nclist:
            nclist.remove(src_c)

        oriQ = self.Modularity()
        dQlist = [] #[[neic1, dQ1], [neic2, dQ2]] pair list

        for neic in nclist:
##            self.move_node(node, neic)
##            dQ = self.Modularity() - oriQ
##
##            self.move_node(node, src_c)
            dQ = self.calc_dQ(nid, neic)
##            if dQ != dQ2:
##                print(dQ, dQ2)



            
            if dQ > 0:
                dQlist.append([neic, dQ])
        self.move_node(node, src_c)

        if not dQlist:
            return None
        else:
            return max(dQlist, key = lambda x:x[1])[0]

    def calc_dQ2(self, nid, dst_cid):
        return self.calc_dQ(nid, self.clist[dst_cid])

    def calc_dQ(self, nid, dst_c):
        m = self.M
        node = self.nlist[nid]
        src_c = node.comm

        src_old_e, src_old_a = src_c.eii, src_c.aii
        dst_old_e, dst_old_a = dst_c.eii, dst_c.aii
        src_new_e = src_old_e - node.degree
        dst_new_e = dst_old_e + node.degree

        src_new_a, dst_new_a = src_old_a, dst_old_a
        for nei, cnt in node.neilist.items():
            if nei == node:
                src_new_a -= cnt
                dst_new_a += cnt
            else:
                if nei.comm == src_c:
                    src_new_a -= cnt * 2
                if nei.comm == dst_c:
                    dst_new_a += cnt * 2

        dQ = (src_new_a - src_old_a + dst_new_a - dst_old_a) / (2 * m) - (src_new_e**2 - src_old_e**2 + dst_new_e**2 - dst_old_e**2) / (4 * m * m)

        return dQ
        
    def all_pos_dQ(self, layer, nid):
        node = self.nlist[nid]
        src_c = node.comm
        
        nvlist = node.neilist.keys()
        nclist = set([n.comm for n in nvlist])
        if src_c in nclist:
            nclist.remove(src_c)

        oriQ = self.Modularity()
        dQlist = [] #[[neic1, dQ1], [neic2, dQ2]] pair list

        for neic in nclist:
            dQ = self.calc_dQ(nid, neic)
            
            if dQ > 0:
                dQlist.append([neic, dQ])

        return dQlist


        node = self.nlist[layer][nid]
        src_c = node.comm
        
        nc_list = set([n.comm for n in node.nei_list])
        if src_c in nc_list:
            nc_list.remove(src_c)
        dQ_list = self.calc_dQlist(layer, nid, nc_list)
        return [x for x in dQ_list if x[1] > 0]

            
        


    def reach_minimal(self):
        node_picker = gennodeseq(1000, self.N)
        moved = False

        while 1:
            next_n = node_picker()
            if next_n == None:
                newc = [c for c in self.clist if not c.isempty()]
                for i, c in enumerate(newc):
                    c.CID = i
                self.clist = newc
                return moved

            layer, nid = next_n
            argmax_c = self.calc_argmaxC(nid)
            if argmax_c != None:
                moved = True
                self.move_node(self.nlist[nid], argmax_c)
                node_picker(rst = True)

    def gen_cinfo(self):
        cinfo = [[n.NID for n in c.nodes] for c in self.clist]
        return cinfo
            
def gen_E_from_C(E, grph):
    nl = grph.nlist
    E2 = [[nl[e0].comm.CID, nl[e1].comm.CID] for e0, e1 in E]
    return E2

def retrieve_c(C, CC):
    retrC = [sum([C[x] for x in c], []) for c in CC]
    return retrC
        
        
def FastUnfolding(E):
    print('enter')
    g = Graph()
    g.updateE(E)
    moved = g.reach_minimal()

    if moved == True:
        cc = FastUnfolding(gen_E_from_C(E, g))
        c = retrieve_c(g.gen_cinfo(), cc)
        return c
    else:
        return g.gen_cinfo()
        


    


    



##gg = Graph()
##gg.readfile('./Research/artificial network/newman/sample_networks/karate.txt')
##print(gg.Modularity())
##
##n = [Node(x) for x in range(3)]
##n[0].addnei(n[1]), n[1].addnei(n[2]), n[2].addnei(n[0])
##n[1].addnei(n[0]), n[0].addnei(n[2]), n[2].addnei(n[1])
##c = [Community(0, n[0:2]), Community(1, [n[2]])]
##[_c.update_neiCCount() for _c in c]
