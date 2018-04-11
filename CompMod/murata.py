from itertools import chain
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
    def __init__(self, NID, layer):
        self.NID = NID
        self.layer = layer
        self.adj_list = Counter()
        self.nei_list = set()
        self.degree = 0
        self.comm = None

    def __repr__(self):
        return "N {}-{}".format(self.layer, self.NID)

    def addlink(self, adjnode):
        self.degree += 1
        self.adj_list[adjnode] += 1

    def gen_nei_list(self):
        self.nei_list.clear()
        for adjn in self.adj_list:
            self.nei_list.update(adjn.adj_list.keys())

        
        self.nei_list.remove(self)


class Community:
    def __init__(self, CID, layer, nodes = set()):
        self.CID = CID
        self.layer = layer
        self.nodes = set()
        self.aCCount = Counter()

        self.degree = 0
        self.partner = None
        
        for node in nodes:
            self.nodes.add(node)
            node.comm = self

        

    def __len__(self):
        return len(self.nodes)

    def isempty(self):
        return len(self.nodes) == 0

    def __repr__(self):
        return "C {}-{}".format(self.layer, self.CID)

    def gen_aCCount(self):
        self.aCCount.clear()
        for node in self.nodes:
            for an, cnt in node.adj_list.items():
                self.aCCount[an.comm] += cnt

    def gen_degree(self):
        self.degree = sum(self.aCCount.values())

    def gen_partner(self):
        if self.isempty():
            self.partner = None
            return

        maxcid = max(self.aCCount, key = self.aCCount.get)
        maxcount = self.aCCount[maxcid]
        partners = [c for c, cnt in self.aCCount.items() if cnt == maxcount]
        mindegree = min([c.degree for c in partners])

        for c in partners:
            if c.degree == mindegree:
                self.partner = c
                return

    def elm(self):
        return self.aCCount[self.partner]


class Graph:
    def __init__(self):
        self.E = None
        self.M = 0
        self.NN = [None] * 2
        self.nlist = [None] * 2
        self.clist = [None] * 2

    def readfile(self, fn):
        E = readnet(fn)
        self.updateE(E)

    def updateE(self, E):
        self.E = E
        self.M = len(E)
        self.NN = [max(x) + 1 for x in zip(*E)]

        self.nlist = [[Node(i, layer) for i in range(n)] for layer, n in enumerate(self.NN)]
        for e in E:
            n0 = self.nlist[0][e[0]]
            n1 = self.nlist[1][e[1]]
            n0.addlink(n1)
            n1.addlink(n0)

        for nl in self.nlist:
            for n in nl:
                n.gen_nei_list()

        self.initC()

    def initC(self):
        self.clist = [[Community(n.NID, n.layer, nodes = [n]) for n in nl] for nl in self.nlist]

        for cl in self.clist:
            for c in cl:
                c.gen_aCCount()
                c.gen_degree()
                
        for cl in self.clist:
            for c in cl:                
                c.gen_partner()

    def updateC(self, cx, cy):
        self.clist = []
        for layer, cl in enumerate([cx, cy]):
            self.clist.append([Community(cid, layer, nodes = [self.nlist[layer][n] for n in c]) for cid, c in enumerate(cl)])
                
        for cl in self.clist:
            for c in cl:
                c.gen_aCCount()
                c.gen_degree()
                
        for cl in self.clist:
            for c in cl:                
                c.gen_partner()        
                
    def MurataQ(self):
        m = self.M
        elm = 0
        alam = 0
        for cl in self.clist:
            for c in cl:
                if not c.isempty():
                    elm += c.elm()
                    alam += c.degree * c.partner.degree

        return (elm / m - alam / (m**2)) / 2

    def move_node(self, node, dst_c):
        src_c = node.comm

        node.comm = dst_c
        src_c.nodes.remove(node)
        dst_c.nodes.add(node)

        src_c.degree -= node.degree
        dst_c.degree += node.degree
        
        account = Counter()
        for n, cnt in node.adj_list.items():
            account[n.comm] += cnt

        src_c.aCCount -= account
        dst_c.aCCount += account

        for ac, cnt in account.items():
            ac.aCCount[src_c] -= cnt
            if ac.aCCount[src_c] == 0:
                del ac.aCCount[src_c]
            ac.aCCount[dst_c] += cnt
            
        src_c.gen_partner()
        dst_c.gen_partner()

        assc_list = set(src_c.aCCount.keys()).union(dst_c.aCCount.keys())
        for assc in assc_list:
            assc.gen_partner()

    def gen_cinfo(self):
        cinfo = []
        for cl in self.clist:
            cinfo.append([[n.NID for n in c.nodes] for c in cl])

        return cinfo

    def mv_nd(self, layer, NID, dst_CID):
        assert layer <= len(self.nlist) - 1
        self.move_node(self.nlist[layer][NID], self.clist[layer][dst_CID])

    def calc_dQ(self, layer, nid, dst_cid):
        node = self.nlist[layer][nid]
        src_c = node.comm
        dst_c = self.clist[layer][dst_cid]

        oriQ = self.MurataQ()
        self.move_node(node, dst_c)
        dQ = self.MurataQ() - oriQ
        self.move_node(node, src_c)

        return dQ

    def Modularity(self):
        return self.MurataQ()


    def calc_dQlist(self, layer, nid, neiC_list):
        node = self.nlist[layer][nid]
        src_c = node.comm

        oriQ = self.MurataQ()
        dQ_list = []
        for nei_c in neiC_list:
            self.move_node(node, nei_c)
            dQ = self.MurataQ() - oriQ
            if dQ > 0:
                dQ_list.append([nei_c, dQ])
        self.move_node(node, src_c)

        return dQ_list

    def all_pos_dQ(self, layer, nid):
        node = self.nlist[layer][nid]
        src_c = node.comm
        
        nc_list = set([n.comm for n in node.nei_list])
        if src_c in nc_list:
            nc_list.remove(src_c)
        dQ_list = self.calc_dQlist(layer, nid, nc_list)
        return [x for x in dQ_list if x[1] > 0]


    def calc_argmaxC(self, layer, nid):
        node = self.nlist[layer][nid]
        src_c = node.comm
        
        nc_list = set([n.comm for n in node.nei_list])
        if src_c in nc_list:
            nc_list.remove(src_c)

        dQ_list = self.calc_dQlist(layer, nid, nc_list)
        if not dQ_list:
            return None
        else:
            return max(dQ_list, key = lambda x:x[1])[0]




    

    def reach_minimal(self):
        node_picker = gennodeseq(1000, *self.NN)
        moved = False

        while 1:
            nextnode = node_picker()
            if not nextnode:
                cx = [c for c in self.clist[0] if not c.isempty()]
                for i, c in enumerate(cx):
                    c.CID = i
                cy = [c for c in self.clist[1] if not c.isempty()]
                for i, c in enumerate(cy):
                    c.CID = i
                self.clist = [cx, cy]
                return moved

            layer, nid = nextnode
            argmax_c = self.calc_argmaxC(layer, nid)
            if argmax_c != None:
                moved = True
                self.mv_nd(layer, nid, argmax_c.CID)
                node_picker(rst = True)
        


def gen_E_from_C(E, grph):
    nl0, nl1 = grph.nlist
    E2 = [[nl0[e0].comm.CID, nl1[e1].comm.CID] for e0, e1 in E]
    return E2

def retrieve_c(Cxy, CCxy):
    CCx, CCy = CCxy
    Cx = [list(flatten_list([Cxy[0][x] for x in c])) for c in CCx]
    Cy = [list(flatten_list([Cxy[1][x] for x in c])) for c in CCy]
    return (Cx, Cy)

    
def FastUnfolding(E):
    g = Graph()
    g.updateE(E)
    moved = g.reach_minimal()

    if moved == True:
        ccxy = FastUnfolding(gen_E_from_C(E, g))
        cxy = retrieve_c(g.gen_cinfo(), ccxy)
        return cxy
    else:
        return g.gen_cinfo()


def MBFU(netfn):
    E = readnet(netfn)
    C = FastUnfolding(E)
    newC = [[sorted(c) for c in sideC]for sideC in C]
    return newC
