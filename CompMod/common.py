from itertools import chain

def flatten_list(lst):
    return list(chain.from_iterable(lst))

def readnet(fn):
    with open(fn) as f:
        E = [list(map(int, line.split())) for line in f]
    return E

def savenet(net, fn):
    with open(fn, 'w') as f:
        for edge in net:
            temp = f.write(' '.join(map(str, edge))+'\n')
        f.flush()

def printnet(net):
    for edge in net:
        print(' '.join(map(str, edge)))

def readdict(fn):
    return eval(open(fn).read())

def writedict(d, fn):
    with open(fn, 'w') as f:
        f.write('{\n')
        for r in sorted(d.keys()):
            f.write("'"+r+"'"+':'+str(d[r])+',\n')
        f.write('}\n')

def range_to_pair(rnglist):
    '''[10, 5, 10] --> [(0, 10), (10, 15), (15, 25)]'''
    base = 0
    r_pair = []
    for rng in rnglist:
        upper = base + rng
        r_pair.append(tuple([base, upper]))
        base = upper
    return r_pair

def belong_judger(rnglist):
    '''rnglist = [5,5,5]; return a func b b(4) = 0; b(7) = 1'''
    pairlist = range_to_pair(rnglist)
    def belongto(t):
        '''if t within a pair, return the position of pair'''
        for i, pair in enumerate(pairlist):
            base, upper = pair
            if  base <= t < upper:
                return i
        print("{} not in {}".format(t, pairlist))
        raise 
    return belongto
