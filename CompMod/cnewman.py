import subprocess
from subprocess import Popen, PIPE
from divtocmu import div_to_com, read_div, divlist_to_cmu
from sys import argv
import os
from common import *
from hierarchy import tree_to_cmu
from newman import FastUnfolding
startupinfo = None
if os.name == 'nt':
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
import datetime
import os

##class FrenchError(Exception):
##    def __init__(self, value):
##        self.value = value
##    def __str__(self):
##        return self.value

##def NuMn(uninet):
##    'newman on subuni, generage resfn(su.)[bin tree div l1]'
##    t = datetime.datetime.now()
##    resfn = uninet.rstrip('.net')+'-res-'+t.strftime('%y%m%d%H%M%S%f')+'.xxx'
##
##    binfn = resfn.replace('xxx', 'bin')
##    treefn = resfn.replace('xxx', 'tree')
##    p = Popen(['newman/convert.exe', '-i', uninet, '-o', binfn], startupinfo = startupinfo)
##    p.wait()
##    p = Popen(['newman/community.exe', '-l', '-1', binfn], startupinfo = startupinfo, \
##              stdout=open(treefn,'w'), stderr=open(os.devnull,'w'))
##    p.wait()
##    
##    tree = readnet(treefn)
##    cmu = tree_to_cmu(tree)
##    os.remove(binfn)
##    os.remove(treefn)
##    return cmu
##    p = Popen(['newman/hierarchy.exe', treefn], startupinfo = startupinfo, stdout = PIPE)
##    p.wait()
##    lvinfo = p.stdout.read().decode()
##    lvstr = lvinfo.split('\n')[0]
##    assert lvstr.startswith('Number of levels: ')
##    lv = int(lvstr.split(' ')[-1]) - 1
##    divfn = resfn.replace('xxx', 'div')
##    p = Popen(['newman/hierarchy.exe', '-l', str(lv), treefn], startupinfo = startupinfo, stdout = open(divfn, 'w'))
##    p.wait()
##    hierinfo = open(divfn).read()
##    divlist = [x.split(' ')[1] for x in hierinfo.strip().split('\n')]
##    div = list(map(int, divlist))
##    os.remove(binfn)
##    os.remove(treefn)
##    os.remove(divfn)
##    
##    #div_to_com
##    maxCID = max(div)
##    d = [[] for x in range(maxCID+1)]
##    for i, x in enumerate(div):
##        d[x].append(i)
##
##    return d

##def NuMnEW(WE,  dirfn='.'):
##    'dirfn is used to save tmpfile. weighted only'
##    wlist = []
##    E = []
##    for n0, n1, w in WE:
##        E.append([n0, n1])
##        wlist.append(w)
##    
##    
##    nodeset = set(flatten_list(E))
##    old_nodes = sorted(nodeset)
##    ondict = {old_nid:new_nid for new_nid, old_nid in enumerate(old_nodes)}
##    nodict = {v:k for k, v in ondict.items()}
##
##    compressedE = [[ondict[n] for n in e]for e in E]
##    compressedWE = [e+[w] for e, w in zip(compressedE, wlist)]
##
##
##
##    
##    t = datetime.datetime.now()
##    tempfn = dirfn+'/temp_numn_subnet_'+t.strftime('%y%m%d%H%M%S%f')+'.net'
##    savenet(compressedWE, tempfn)
##    try:
##        compressedcmu = NuMnW(tempfn)
##    except FrenchError:
##        print('FE')
##        compressedcmu = FastUnfolding(compressedWE)
##    
##    cmu = [[nodict[n] for n in c] for c in compressedcmu]
##
##    os.remove(tempfn)
##    return cmu
##
##
##def NuMnW(uninet):
##    'weighted version of NuMn'
##    t = datetime.datetime.now()
##    resfn = uninet.rstrip('.net')+'-res-'+t.strftime('%y%m%d%H%M%S%f')+'.xxx'
##
##    binfn = resfn.replace('xxx', 'bin')
##    p = Popen(['newman/convert.exe', '-i', uninet, '-o', binfn, '-w'], startupinfo = startupinfo)
##    p.wait()
##    p = Popen(['newman/community.exe', '-l', '-1', binfn, '-w'], startupinfo = startupinfo, \
##              stdout=PIPE, stderr=PIPE)
##
##    treestr = ''
##
##    for perr in p.stderr:
##        perr = perr.decode()
##        if perr.startswith('modularity'):
##            treestr += p.stdout.read().decode()
##        if perr.startswith('pass number') and \
##           (float(perr.split()[-1]) > 1 or int(perr.split()[2]) > 500):
##            p.kill()
##            os.remove(binfn)
##            raise FrenchError('The fucking stupid French programming is burning')
##    
##    treestr += p.stdout.read().decode()
##    tree = [[int(n) for n in e.split()]for e in treestr.split('\n') if e]
##    try:
##        cmu = tree_to_cmu(tree)
##    except IndexError as e:
##        raise FrenchError(str(e))
##    os.remove(binfn)
##    return cmu

def BNF(binet):
    "BiNetFinder on subbi, generate resfn(sb.)[div0, div1, l0 l1]"
    t = datetime.datetime.now()
    resfn = binet.rstrip('.net')+'-res-'+t.strftime('%y%m%d%H%M%S%f')+'.xxx'
    print("resfn: {0}".format(resfn))
    assert resfn.count('.') == 1
    div0 = resfn.replace('xxx', 'div0')
    div1 = resfn.replace('xxx', 'div1')

    p = Popen(['BiNetFinder.exe', binet, div0, div1], startupinfo = startupinfo)
    p.wait()

    c0 = read_div(div0)
    c1 = read_div(div1)

    os.remove(div0)
    os.remove(div1)
    return [c0, c1]

def TNF(trinet):
    "TriNetFinder on subtri"
    t = datetime.datetime.now()
    resfn = trinet.rstrip('.net')+'-res-'+t.strftime('%y%m%d%H%M%S%f')+'.xxx'

    assert resfn.count('.') == 1
    div0 = resfn.replace('xxx', 'div0')
    div1 = resfn.replace('xxx', 'div1')
    div2 = resfn.replace('xxx', 'div2')

    p = Popen(['TriNetFinder.exe', trinet, div0, div1, div2], startupinfo = startupinfo)
    p.wait()

    c0 = read_div(div0)
    c1 = read_div(div1)
    c2 = read_div(div2)


    os.remove(div0)
    os.remove(div1)
    os.remove(div2)
    return [c0, c1, c2]

def LPMEW(WE,  dirfn='.'):
    'dirfn is used to save tmpfile. weighted only'
    wlist = []
    E = []
    for n0, n1, w in WE:
        E.append([n0, n1])
        wlist.append(w)
    
    
    nodeset = set(flatten_list(E))
    old_nodes = sorted(nodeset)
    ondict = {old_nid:new_nid for new_nid, old_nid in enumerate(old_nodes)}
    nodict = {v:k for k, v in ondict.items()}

    compressedE = [[ondict[n] for n in e]for e in E]
    compressedWE = [e+[w] for e, w in zip(compressedE, wlist)]



    
    t = datetime.datetime.now()
    tempfn = dirfn+'/temp_numn_subnet_'+t.strftime('%y%m%d%H%M%S%f')+'.net'
    savenet(compressedWE, tempfn)
    compressedcmu = LPMW(tempfn)

    #decompress
    cmu = [[nodict[n] for n in c] for c in compressedcmu]

    os.remove(tempfn)
    return cmu

def LPM(uninet):
    E = readnet(uninet)
    dirfn = os.path.dirname(uninet)
    WE = [[e0, e1, 1] for e0, e1 in E]
    return LPMEW(WE, dirfn)
    
    
def LPMW(weighted_uninet):
    "LPM, on weighted uninet"
    p = Popen(['LPMp.exe', weighted_uninet], stdout = PIPE, startupinfo = startupinfo)
    started = False
    div = []
    for pout in p.stdout:
        pout = pout.decode()
        if started:
            div.append(int(pout))
        if pout.startswith('--START--'):
            started = True
    return divlist_to_cmu(div)
