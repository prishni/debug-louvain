from itertools import chain, groupby
from collections import Counter
from random import choice, shuffle, randint
from copy import deepcopy
import newman, murata, muratatri
from sys import argv
from common import *
graph_dict = {'uni':newman,'bi':murata,'tri':muratatri}


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




def subgraphtypefinder(E, noderange_list):
    '''input E and layerinfo, determine the type of g'''

    
    if len(E) == 0:
        raise
    node_to_layer = belong_judger(noderange_list)
    e = E[0]
    c_of_e = [node_to_layer(n) for n in e]
    assert (None not in c_of_e)
    if len(e) == 3 and len(set(c_of_e)) == 3:
        gtype = 'tri'
        layerinfo = tuple(c_of_e)
    elif len(e) == 2 and len(set(c_of_e)) == 2:
        gtype = 'bi'
        layerinfo = tuple(c_of_e)
    elif len(e) == 2 and len(set(c_of_e)) == 1:
        gtype = 'uni'
        layerinfo = tuple(c_of_e)
    else:
        raise

    for e in E:
        try:
            assert [node_to_layer(n) for n in e] == list(layerinfo)
        except:
            print([node_to_layer(n) for n in e], list(layerinfo))
            raise

    return gtype, layerinfo


        

class Heterograph:
    def __init__(self, E, linkrange_list, noderange_list):
        self.linkrange_list = linkrange_list[:]
        self.noderange_list = noderange_list[:]
        # seperate E into different subgraphs
        self.E_list = [E[base:upper] for base, upper in range_to_pair(linkrange_list)]

        # use size of graph as weight
        self.w_list = linkrange_list[:]

        # create subgraph based on subE and layer information(noderange_list)
        self.subg_list = [SubGraph(subE, noderange_list) for subE in self.E_list]

        # init global_clist
        self.c_list = [[x] for x in range(sum(noderange_list))]

    def reach_minimal(self):
        node_picker = gennodeseq(500, sum(self.linkrange_list))
        moved = False
        while 1:
            next_n = node_picker()
            if next_n == None:
                newc = [c for c in self.c_list if c]
                self.c_result = newc

                return moved

            layer, g_nid = next_n #layer is meaningless
            argmax_cid = self.calc_amc(g_nid)
            if argmax_cid != None:
                moved = True


                for subg in self.subg_list:
                    subg.mvnd(g_nid, argmax_cid[0])
                self.pretend_mvnd(g_nid, argmax_cid[0])

                node_picker(rst = True)

    def calc_amc(self, g_nid): #arg max c
        posc = []
        for subg in self.subg_list:
            posc.extend(subg.all_pos_dQc(g_nid))

        if not posc:
            return None

        dQ_list = []

        for g_cid in list(zip(*posc))[0]:
            
            dQ = sum([w*subg.caldq(g_nid, g_cid) for w, subg in zip(self.w_list, self.subg_list)])
            dQ_list.append(tuple([g_cid, dQ/sum(self.w_list)]))


        max_c_dq = max(dQ_list, key = lambda x:x[1])
        if max_c_dq[1] > 0:
            return max_c_dq#[0]
        else:
            return None
            
    def Modu(self):
        return sum([w * subg.g.Modularity() for w, subg in zip(self.w_list, self.subg_list)]) / sum(self.w_list)
        

    def pretend_mvnd(self, g_nid, g_cid):
        for c in self.c_list:
            if g_nid in c:
                c.remove(g_nid)
        self.c_list[g_cid].append(g_nid)
        
            
            

        
        

        

        


class SubGraph:
    def __init__(self, E, noderange_list):
        self.global_E = deepcopy(E)
        self.gtype, self.layerinfo = subgraphtypefinder(E, noderange_list)

        ## select module and instantiate
        self.g = graph_dict[self.gtype].Graph()

        self.glt_n = dict() # global_local_table

        ## generate glt_n
        # case uni
        if self.gtype == 'uni':
            layer = self.layerinfo[0]
            for i, n in enumerate(sorted(set(sum(E, [])))):
                self.glt_n[n] = (0, i)
        # case bi / tri
        else:
            nodebag = [sorted(set(ns)) for ns in zip(*E)]
            for locallayer, ns in enumerate(nodebag):
                globallayer = self.layerinfo[locallayer]
                for i, n in enumerate(ns):
                    self.glt_n[n] = (locallayer, i)


        ## generate glt_c
        self.glt_c = dict()
        
        if self.gtype == 'uni':
            g_layer = self.layerinfo[0]
            pair = range_to_pair(noderange_list)[g_layer]
            for l_cid, g_cid in enumerate(range(*pair)):
                self.glt_c[g_cid] = (0, l_cid)
        else:
            for l_layer, g_layer in enumerate(self.layerinfo):
                pair = range_to_pair(noderange_list)[g_layer]
                for l_cid, g_cid in enumerate(range(*pair)):
                    self.glt_c[g_cid] = (l_layer, l_cid)

        # generate lgt_c / lgt_n
        self.lgt_n = {v:k for k, v in self.glt_n.items()} # local_global_table = reverse global_local_table
        self.lgt_c = {v:k for k, v in self.glt_c.items()}

        # generate real_c, use for update
        if self.gtype == 'uni':
            g_layer = self.layerinfo[0]
            pair = range_to_pair(noderange_list)[g_layer]
            real_c = [[] for x in range(*pair)]
            for cid, c in enumerate(real_c):
                g_cid = self.lgt_c[(0, cid)]
                g_nid = g_cid
                if g_nid in self.glt_n:
                    c.append(self.glt_n[g_nid][1])
            real_cl = [real_c]

        else:
            real_cl = []
            for l_layer, g_layer in enumerate(self.layerinfo):
                pair = range_to_pair(noderange_list)[g_layer]
                real_c = [[] for x in range(*pair)]
                for cid, c in enumerate(real_c):
                    g_cid = self.lgt_c[(l_layer, cid)]
                    g_nid = g_cid
                    if g_nid in self.glt_n:
                        c.append(self.glt_n[g_nid][1])
                real_cl.append(real_c)

            
            


                    


        self.local_E = [list(map(lambda x:self.glt_n.get(x)[1], e)) for e in self.global_E]
        self.g.updateE(self.local_E)
        self.g.updateC(*real_cl)
            
            

    def __repr__(self):
        return "{} SubGraph in layer {}".format(self.gtype, self.layerinfo)

    def all_pos_dQc(self, g_nid):
        if g_nid not in self.glt_n:
            return []
        layer, nid = self.glt_n[g_nid]
        dqlist = self.g.all_pos_dQ(layer, nid)

        #return list of g_cid where dq is positive
        return [(self.lgt_c[(c.layer, c.CID)], dq) for c, dq in dqlist] 

    def caldq(self, g_nid, g_cid):
        if g_nid not in self.glt_n:
            return 0
        layern, nid = self.glt_n[g_nid]
        try:
            layerc, cid = self.glt_c[g_cid]
        except:
            print(g_nid, g_cid, self.glt_n, glt_c)
            raise

        assert layern == layerc
        layer = layern

        try: # compliant with newman
            return self.g.calc_dQ(layer, nid, cid)
        except:
            return self.g.calc_dQ2(nid, cid)

    def mvnd(self, g_nid, g_cid):
        if g_nid not in self.glt_n:
            return
        layer, nid = self.glt_n[g_nid]
        layer, cid = self.glt_c[g_cid]
        if self.gtype == 'uni':
            self.g.mv_nd(nid, cid)
        elif self.gtype in ('bi', 'tri'):
            self.g.mv_nd(layer, nid, cid)
        else:
            raise


def gen_E_from_C(E, C):
    c_of_n = gen_cofn_from_c(C)
    E2 = [[c_of_n[n] for n in e] for e in E]
    return E2

def retrieve_c(C, CC):
    retrC = [sum([C[x] for x in c], []) for c in CC]
    return retrC

def gen_cofn_from_c(C):
    c_of_n = dict()
    for cid, c in enumerate(C):
        for n in c:
            c_of_n[n] = cid

    return c_of_n

def gen_nr_from_c(nr, C):
    c_of_n = gen_cofn_from_c(C)
    l_of_n = belong_judger(nr)

    new_nr = [0] * len(nr)
    for c in C:
        lset = set([l_of_n(n) for n in c])
        assert len(lset) == 1
        new_nr[lset.pop()] += 1
    return new_nr

            
def FastUnfolding(E, lr, nr):
    print('enter FU')
    g = Heterograph(E, lr, nr)
    moved = g.reach_minimal()

    if moved == True:
        new_nr = gen_nr_from_c(nr, g.c_result)
        new_E = gen_E_from_C(E, g.c_result)
        cc = FastUnfolding(new_E, lr, new_nr)
        c = retrieve_c(g.c_result, cc)
        return c
    else:
        return g.c_result
        

        
        
        
##E = [[0,1],[1,2],[7,3,0],[8,5,2],[7,8],[4,6],[5,6],[3,6]]
##lr = [2,2,1,3]
##nr = [3,3,3]
##        
##hg = Heterograph(E, lr, nr)



def main(hnet, metafn, resfn):
    E = readnet(hnet)
    temp = open(metafn).read()
    exec(temp, globals())
    res = FastUnfolding(E, lr, nr)
##    printnet([sorted(x) for x in res])
    savenet([sorted(x) for x in res], resfn)
    return [sorted(x) for x in res]

def combi_cd(dirfn):
    C = main(dirfn+'/hetero.net', dirfn+'/.meta', dirfn+'/combi.cmu')
    return C


if __name__ == '__main__':
    if len(argv) < 4:
        print('argv: combiner.py hnet meta res.cmu')
    else:
        main(argv[1], argv[2], argv[3])

##
##E = readnet('./hetero.net')
##lr = [124, 188-124]
##nr = [12, 12, 18]
##hg = Heterograph(E, lr, nr)
##E1 = E[:lr[0]] # 84 nodes, 118 line, tripartite
##E2 = E[lr[0]:] # 28 nodes, ? line, unipartite
##
##d1 = [{e: i for i, e in enumerate(set(nr))} for nr in zip(*E1)]
##d2 = {e:i for i, e in enumerate(sorted(set(sum(E2, []))))}
##
##newE1 = [[d1[i][n] for i, n in enumerate(e)] for e in E1]
##newE2 = [list(map(d2.get, e)) for e in E2]
##
##import newman, muratatri, murata
##
##g1 = muratatri.Graph()
##g2 = newman.Graph()
##
##g1.updateE(newE1)
##g2.updateE(newE2)
##
##NN = 24,24,24
##node_picker = gennodeseq(1000, *NN)
