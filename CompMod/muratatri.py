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

    def addlink(self, adjnodes):
        self.degree += 1
        self.adj_list[adjnodes] += 1

    def gen_nei_list(self):
        self.nei_list.clear()

        adjset = set()

        # put all adjn of 2 layers in 1 set
        for adjns in self.adj_list:
            adjset.update(adjns)

        neiset = set()
        for adjn in adjset:
            for neins in adjn.adj_list:
                neiset.update(neins)

        self.nei_list = set(n for n in neiset if n.layer == self.layer)
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
                self.aCCount[(an[0].comm, an[1].comm)] += cnt

    def gen_degree(self):
        self.degree = sum(self.aCCount.values())

    def gen_partner(self):
        if self.isempty():
            self.partner = None
            return

        maxcount = max(self.aCCount.values())
        partners = [c for c, cnt in self.aCCount.items() if cnt == maxcount]
        mindegree = min([c[0].degree * c[1].degree for c in partners])

        for c in partners:
            if c[0].degree * c[1].degree == mindegree:
                self.partner = c
                return

    def elm(self):
        return self.aCCount[self.partner]


class Graph:
    def __init__(self):
        self.E = None
        self.M = 0
        self.NN = [None] * 3
        self.nlist = [None] * 3
        self.clist = [None] * 3

    def readfile(self, fn):
        E = readnet(fn)
        self.updateE(E)

    def updateE(self, E):
        self.E = E
        self.M = len(E)
        self.NN = [max(x) + 1 for x in zip(*E)]

        self.nlist = [[Node(i, layer) for i in range(n)] for layer, n in enumerate(self.NN)]
        for e in E:
            n0, n1, n2 = [self.nlist[i][ei] for i, ei in enumerate(e)]
            n0.addlink(tuple([n1, n2]))
            n1.addlink(tuple([n0, n2]))
            n2.addlink(tuple([n0, n1]))

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

    def updateC(self, *c):
        self.clist = []
        for layer, cl in enumerate(c):
            self.clist.append([Community(cid, layer, nodes = [self.nlist[layer][n] for n in c]) for cid, c in enumerate(cl)])
                
        for cl in self.clist:
            for c in cl:
                c.gen_aCCount()
                c.gen_degree()
                
        for cl in self.clist:
            for c in cl:                
                c.gen_partner()        
                
    def MuratatriQ(self):
        m = self.M
        elm = 0
        alam = 0
        for cl in self.clist:
            for c in cl:
                if not c.isempty():
                    elm += c.elm()
                    alam += c.degree * c.partner[0].degree * c.partner[1].degree

        return (elm / m - alam / (m**3)) / 3

    def Modularity(self):
        return self.MuratatriQ()

    def move_node(self, node, dst_c):
        assert node.layer == dst_c.layer
        src_c = node.comm

        node.comm = dst_c
        src_c.nodes.remove(node)
        dst_c.nodes.add(node)

        src_c.degree -= node.degree
        dst_c.degree += node.degree
        
        account = Counter()
        for nl, cnt in node.adj_list.items():
            account[tuple(n.comm for n in nl)] += cnt

        # assertion
        for cl, cnt in account.items():
            if not src_c.aCCount[cl] >= cnt:
                print(cl, src_c.aCCount, account)
                raise
        # assertion over
        
        src_c.aCCount -= account
        dst_c.aCCount += account

        # assertion
        tmp = Counter(src_c.aCCount)
        src_c.gen_aCCount()
        assert tmp == src_c.aCCount
        tmp = Counter(dst_c.aCCount)
        dst_c.gen_aCCount()
        assert tmp == dst_c.aCCount        
        # assertion over

        for acs, cnt in account.items():
            for ac in acs:
                another_ac = acs[1] if ac == acs[0] else acs[0]
                old_lnk = (src_c, another_ac) if src_c.layer < another_ac.layer else (another_ac, src_c)
                new_lnk = (dst_c, another_ac) if dst_c.layer < another_ac.layer else (another_ac, dst_c)
                ac.aCCount[old_lnk] -= cnt
                if ac.aCCount[old_lnk] == 0:
                    del ac.aCCount[old_lnk]

                ac.aCCount[new_lnk] += cnt

        for acs, cnt in account.items():
            for ac in acs:
                # assertion
                tmp = Counter(ac.aCCount)
                ac.gen_aCCount()
                try:
                    assert tmp == ac.aCCount
                except:
                    print(ac)
                    print(tmp)
                    print(ac.aCCount)
                    raise


                
            
        src_c.gen_partner()
        dst_c.gen_partner()

        assc_list = set()
        for c in src_c, dst_c:
            for x in c.aCCount.keys():
                assc_list.update(x)
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


    def calc_dQlist(self, layer, nid, neiC_list):
        node = self.nlist[layer][nid]
        src_c = node.comm

        oriQ = self.MuratatriQ()
        dQ_list = []
        for nei_c in neiC_list:
            self.move_node(node, nei_c)
            dQ = self.MuratatriQ() - oriQ
            dQ_list.append([nei_c, dQ])
        self.move_node(node, src_c)

        return dQ_list

    def calc_dQ(self, layer, nid, dst_cid):
        node = self.nlist[layer][nid]
        src_c = node.comm
        dst_c = self.clist[layer][dst_cid]

        oriQ = self.MuratatriQ()
        self.move_node(node, dst_c)
        dQ = self.MuratatriQ() - oriQ
        self.move_node(node, src_c)

        return dQ

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

##
##            
##        nc_list = set(x for x in self.clist[layer] if not x.isempty())
##        nc_list.remove(src_c)
##
##
        
        dQ_list = self.calc_dQlist(layer, nid, nc_list)
        if not dQ_list:
            return None
        elif max([x[1] for x in dQ_list]) < 0:
            return None
        else:
            return max(dQ_list, key = lambda x:x[1])[0]




    

    def reach_minimal(self):
        node_picker = gennodeseq(1000, *self.NN)
        moved = False

        while 1:
            nextnode = node_picker()
            if not nextnode:
                newclist = []
                for layer, cl in enumerate(self.clist):
                    newcl = [c for c in cl if not c.isempty()]
                    for i, c in enumerate(newcl):
                        c.CID = i
                    newclist.append(newcl)
                self.clist = newclist
                return moved

            layer, nid = nextnode
            argmax_c = self.calc_argmaxC(layer, nid)
            if argmax_c != None:
                moved = True
                self.mv_nd(layer, nid, argmax_c.CID)
                node_picker(rst = True)



def gen_E_from_C(E, grph):
    nl = grph.nlist
    E2 = [[nl[i][n].comm.CID for i, n in enumerate(e)] for e in E]
    return E2

def retrieve_c(Cxyz, CCxyz):
    # using sum seems to be more clear
    CCx, CCy, CCz= CCxyz

    Cx = [list(flatten_list([Cxyz[0][x] for x in c])) for c in CCx]
    Cy = [list(flatten_list([Cxyz[1][x] for x in c])) for c in CCy]
    Cz = [list(flatten_list([Cxyz[2][x] for x in c])) for c in CCz]
    return (Cx, Cy, Cz)

    
def FastUnfolding(E):
    g = Graph()
    g.updateE(E)
    moved = g.reach_minimal()
##    print('this layer, modularity = ', g.MuratatriQ())

    if moved == True:
        ccxy = FastUnfolding(gen_E_from_C(E, g))
        cxy = retrieve_c(g.gen_cinfo(), ccxy)
        return cxy
    else:
        return g.gen_cinfo()


def MTFU(netfn):
    E = readnet(netfn)
    C = FastUnfolding(E)
    newC = [[sorted(c) for c in sideC]for sideC in C]
    return newC
    


##EE = [[11, 5, 10], [1, 2, 3], [11, 6, 9], [11, 4, 5], [7, 4, 3], [1, 10, 2], [4, 1, 6], [5, 1, 7], [11, 9, 5], [9, 7, 0], [11, 0, 9], [6, 1, 4], [3, 0, 11], [3, 5, 4], [0, 5, 11], [10, 0, 6], [3, 2, 11], [5, 10, 10], [7, 2, 0], [4, 10, 9], [5, 10, 0], [6, 11, 10], [3, 11, 10], [9, 5, 6], [1, 11, 2], [7, 3, 4], [9, 3, 3], [11, 7, 1], [6, 7, 8], [11, 1, 8], [2, 5, 6], [4, 11, 11], [0, 5, 1], [5, 4, 9], [3, 10, 3], [10, 2, 6], [1, 7, 4], [4, 4, 0], [5, 10, 6], [2, 4, 11], [7, 8, 8], [2, 0, 8], [3, 0, 1], [4, 4, 9], [8, 2, 1], [3, 10, 11], [10, 5, 9], [4, 0, 10], [8, 7, 3], [15, 16, 25], [19, 21, 34], [22, 15, 29], [22, 23, 35], [12, 12, 30], [18, 18, 30], [13, 20, 28], [14, 15, 26], [20, 12, 26], [14, 20, 35], [19, 18, 32], [18, 18, 33], [13, 17, 29], [22, 20, 30], [13, 14, 29], [17, 18, 31], [23, 19, 25], [18, 16, 25], [16, 23, 28], [22, 12, 29], [15, 20, 25], [15, 22, 32], [15, 12, 26], [23, 15, 28], [18, 19, 32], [19, 20, 33], [12, 14, 25], [22, 16, 34], [19, 18, 27], [23, 17, 35], [21, 23, 27], [23, 15, 30], [16, 21, 24], [21, 23, 35], [23, 15, 27], [21, 23, 31], [14, 18, 31], [14, 16, 28], [18, 12, 34], [10, 4, 23], [0, 5, 17], [7, 8, 22], [9, 9, 15], [9, 9, 13], [8, 1, 19], [10, 2, 18], [6, 4, 12], [11, 6, 12], [7, 10, 20], [8, 11, 21], [2, 4, 12], [1, 3, 19], [2, 3, 17], [3, 1, 22], [3, 4, 14], [3, 10, 18], [17, 17, 18], [19, 22, 13], [23, 15, 19], [19, 13, 17], [16, 23, 13], [19, 19, 20], [15, 19, 13], [19, 15, 16], [18, 14, 15], [15, 21, 20], [19, 16, 19], [12, 13, 17], [14, 23, 14]]
##ggg = Graph() ; ggg.updateE(EE)
