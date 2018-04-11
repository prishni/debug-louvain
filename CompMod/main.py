from common import *
from mergecom import merge
from copy import deepcopy
from operator import sub
from cnewman import LPM, BNF, TNF
from murata import MBFU
from muratatri import MTFU
import combiner
import os
from sys import argv


def subgraphtypefinder(SubE, noderange_list):
    '''input E and layerinfo, determine the type of g'''


    if len(SubE) == 0:
        raise
    node_to_layer = belong_judger(noderange_list)
    e = SubE[0]
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

    for e in SubE:
        try:
            assert [node_to_layer(n) for n in e] == list(layerinfo)
        except:
            print([node_to_layer(n) for n in e], list(layerinfo))
            raise

    return gtype, layerinfo

def gen_oldnew_dict(nodeset):
    # new is solid new id
    old_new_dict = {old:new for new, old in enumerate(sorted(set(nodeset)))}
    return old_new_dict

def gen_newold_dict(nodeset):
    new_old_dict = {new:old for new, old in enumerate(sorted(set(nodeset)))}
    return new_old_dict

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

class HGraph:
    def __init__(self, folder):
        'folder as input. folder/hetero.net, folder/.meta as default input file.'
        self.E = readnet(folder+'/hetero.net')
        E = self.E
        with open(folder+'/.meta') as metaf:
            x = metaf.readline().strip().split(' = ')
            assert(x[0] == 'lr')
            lr = eval(x[1])
            x = metaf.readline().strip().split(' = ')
            assert(x[0] == 'nr')
            nr = eval(x[1])
        self.nr, self.lr = nr, lr
        self.SubE_list = (deepcopy(E[a:b]) for a, b in range_to_pair(lr))
        self.SubG_list = [SubGraph(SubE, nr, lr) for SubE in self.SubE_list]
        for i, SubG in enumerate(self.SubG_list):
            SubG.subid = i
            SubG.folder = folder

    def community_detection(self):
        result_layers = {}
        for SubG in self.SubG_list:
            for l, r in SubG.comu_detec().items():
                try:
                    result_layers[l].append(r)
                except:
                    result_layers[l] = [r]

		##        for l in result_layers:
		##            print('layer: ', l)
		##            print(merge(*result_layers[l]))
		##            print()

        result_merged_falldown = [merge(*result_layers[l]) for l in range(len(self.nr))]

        layer_start_id = [x[0] for x in range_to_pair(self.nr)]


        # get nid_falldown into nid, then merge all layers.
        result_merged = []
        for l, clist in enumerate(result_merged_falldown):
            for c in clist:
                result_merged.append([n+layer_start_id[l] for n in c])

		##        print(result_merged)
        new_nr = gen_nr_from_c(self.nr, result_merged)
        new_E = gen_E_from_C(self.E, result_merged)
        cc = combiner.FastUnfolding(new_E, self.lr, new_nr)
        c = retrieve_c(result_merged, cc)
        for x in c: x.sort()
        return c

class SubGraph:
    def __init__(self, subE, nr, lr):
        self.nr, self.lr = nr, lr
        self.subE = deepcopy(subE)
        self.gtype, self.layerinfo = subgraphtypefinder(subE, nr)
        self.falldownE = self.falldown() # using subE, nr, lr, gtype, layerinfo
        self.compressE, self.ond_layers = self.compress() # using falldownE, nr, lr, gtype, layerinfo
        self.nod_layers = {l:{new:old for old, new in d.items()} for l, d in self.ond_layers.items()}
        self.missed_node_layers = self.gen_missed_node()


    def falldown(self):
        SubE, nr, lr = self.subE, self.nr, self.lr
        gtype, layerinfo = self.gtype, self.layerinfo
        # the start id for each layer
        layer_start_id = [x[0] for x in range_to_pair(nr)]

        # offset that each e should falldown
        offsets = [layer_start_id[layer] for layer in layerinfo]

        falldownSubE = [list(map(sub, e, offsets)) for e in SubE]

        return falldownSubE


    def compress(self):
        falldown_SubE = self.falldownE
        nr, lr = self.nr, self.lr
        gtype, layerinfo = self.gtype, self.layerinfo

        nodeset_in_layer = {l:set() for l in set(layerinfo)}

        for e in falldown_SubE:
            for l, n in zip(layerinfo, e):
                nodeset_in_layer[l].add(n)

        on_dict_of_layer = {}
        for l, nodeset in nodeset_in_layer.items():
            on_dict_of_layer[l] = gen_oldnew_dict(nodeset)


        compressed_SubE = [[on_dict_of_layer[l][n] for l, n in zip(layerinfo, e)] for e in falldown_SubE]

        return compressed_SubE, on_dict_of_layer

    def gen_missed_node(self):
        '''missed node, using falldown nid (start from 0 each layer)'''
        ond = self.ond_layers
        nr = self.nr
        missed_node_layers = {l:set(range(nr[l])) - set(d.keys()) for l, d in ond.items()}

        return missed_node_layers

    def comu_detec(self):
		##        using binetfinder in the first step
        cd_method = {'uni':LPM, 'bi':BNF, 'tri':TNF}[self.gtype]

        subnetfn = self.folder+'/sub{}.net'.format(self.subid)
        savenet(self.compressE, subnetfn)
        result = cd_method(subnetfn)
        result_layers = {l:None for l in self.ond_layers}


        # subtle differences between processing uni and multipartite network
        if self.gtype == 'uni':
            l = self.layerinfo[0]
            nod = self.nod_layers[l]
            real_result = [list(map(nod.get, c)) for c in result]
            real_result.append(list(self.missed_node_layers[l])) # append all missed nodes in one community
            return {l:real_result}

        else:
            real_result_layers = {}
            for i, l in enumerate(self.layerinfo):
                nod = self.nod_layers[l]
                real_result = [[nod[n] for n in c] for c in result[i]]
                real_result.append(list(self.missed_node_layers[l])) # append all missed nodes in one community
                real_result_layers[l] = real_result
            return real_result_layers

##exec(open('artnet/.meta').read(), globals())
##
##HE = readnet('artnet/hetero.net')
##SubE0, SubE1, SubE2 = (deepcopy(HE[a:b]) for a, b in range_to_pair(lr))

def main(folderfn, net, resfn=None):
    if resfn == None:
        resfn = folderfn +net+ '_merge.cmu'
    else:
        resfn = folderfn +resfn
    print('read {}hetero.net and {}.meta'.format(folderfn, folderfn))
    print('save result to '+resfn)
    H = HGraph(folderfn)
    C = H.community_detection()
    savenet(C, resfn)
    print('read {}/hetero.net and {}/.meta'.format(folderfn, folderfn))
    print('save result to '+resfn)
    return C

def convertInput(file,tmpfolder):
    fin = open(file)
    fout1 = open(tmpfolder+'hetero.net','w')
    fout2 = open(tmpfolder+'.meta','w')
    fout3 = open(tmpfolder+'ground_truth_communities','w')
    lr =[]
    nr =[100,100]
    fin.readline()  
    fin.readline()
    lr.append(int(fin.readline().strip()))
    i = 0
    while(i<lr[0]):
        line = fin.readline().split()
        e1 = int(float(line[0]))-1
        e2 = int(float(line[1]))-1
        fout1.write(str(e1)+' '+str(e2)+'\n')
        i+=1
    
    fin.readline()
    lr.append(int(fin.readline().strip()))
    i=0
    while(i<lr[1]):
        line = fin.readline().split()
        e1 = int(float(line[0]))-1
        e2 = int(float(line[1]))-1
        fout1.write(str(e1)+' '+str(e2)+'\n')
        i+=1

    fin.readline()  
    fin.readline()
    
    lr.append(int(fin.readline().strip()))
    i=0
    while(i<lr[2]):
        line = fin.readline().split()
        e1 = int(float(line[0]))-1
        e2 = int(float(line[1]))-1
        fout1.write(str(e1)+' '+str(e2)+'\n')
        i+=1

    
    commu = int(float(fin.readline().strip()))
    i=0
    while(i<commu):
        line = fin.readline().strip().split()
        fout3.write(str(int(float(line[0]))-1))    
        for l in line[1:]:
            fout3.write(' '+str(int(float(l))-1))  
        fout3.write('\n')
        i+=1

    fout2.write('lr = ['+str(lr[0])+', '+str(lr[1])+', '+str(lr[2])+']\n' )
    fout2.write('nr = [100, 100]' )

    fout1.close()
    fout2.close()
    fout3.close()

'''
if __name__ == '__main__':
    if not (len(argv) in (2,3)):
        print('usage: main.py folder_name [result_fn]\n"hetero.net" and ".meta" in folder\ndefault result_filename = folder/merge.cmu')
    elif len(argv) == 2:
        temp, folderfn = argv
        main(folderfn)
    else: #len(argv) == 2:
        temp, folderfn, resfn = argv
        resfn = folderfn + '/' + resfn
        main(folderfn, resfn)
'''
if __name__ == '__main__':
    folder_name = "../_networks/"
    tmpfolder = "temphetero/"
    #temp, folderfn = argv    
    networks = os.listdir(folder_name)
    for net in networks:
        convertInput(folder_name+net,tmpfolder)
        main(tmpfolder,net)


    
    