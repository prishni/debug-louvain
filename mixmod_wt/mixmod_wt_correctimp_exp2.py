modctr = 0

def __modularity(commu, status, graph):
    global modctr
    modctr += 1
    #print("modularity called", modctr, "edgewt: ", [graph[1][nbr].get('weight',1) for nbr in graph[1]])
    #print("From modularity, node_c: ", status.node_c)
    #print("From modularity, node_l: ", status.node_l)

    layer=status.layer
    node_l=status.node_l
    node_c=status.node_c       
    top=status.top
    bot=status.bot
    edge_l=status.edge_l
    edge_c=status.edge_c
    couple=status.couple
    mu = status.mu

    #Construct a dict = {node1: layer, node2: layer......}
    nodelayer = {}
    for l in layer:
        for node in layer[l]:
            nodelayer[node] = l

    #calculate Intra_inter ----------------------------------------------------------
    f=0    
    intra_inter={}
    for c in commu:
        intra_inter[c]=set()
        
        for n in commu[c]:
            for l in layer:
                if n in layer[l]:
                    intra_inter[c].add(l)

    #--------------------------------------------------------------------------------

    #calculate |E1| , |E2| , |E12|---------------------------------------------------
    E={}
    Eloop={}
    E12=0
    for l in layer:
        E[l]=0
        Eloop[l]=0
        for n in layer[l]:
            for nei in node_l.get(n,set()):

                E[l]+= graph[n][nei].get('weight',1)
                '''
                if (n==nei): #for loops add weight twice
                    E[l]+= graph[n][nei].get('weight',1)
                '''
        E[l] = E[l]/2

    for n in node_c:
        for nei in node_c[n]:
            E12+=graph[n][nei].get('weight',1)
    E12 = E12/2
    #--------------------------------------------------------------------------------

   
    modularity=0    
    x1={}
    x2={}    
    
    for c in commu:
        x1[c]=0
        x2[c]=0
        modc_layer=0
        if len(intra_inter[c])>1:
            for l in layer:
                Aij=0.0
                hihj=0.0
                
                #compute summation Aij------------------------
                for n in commu[c]:
                    if n in layer[l]:
                        for nei in node_l.get(n,set()):
                            if nei in commu[c]:
                                Aij+=graph[n][nei].get('weight',1)
                                '''
                                if (nei==n):  #for loops add weight twice
                                    Aij+=graph[n][nei].get('weight',1)
                                '''
                Aij = Aij/2

                #compute summation hihj--------------------------
                for n1 in commu[c]:
                    if n1 in layer[l]:
                        for n2 in commu[c]:
                            if n2 in layer[l]:
                                if(n1==n2):
                                    if(n1 not in node_l[n1]): 
                                        continue
                                '''
                                    else:
                                        hi = 2*sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                                        hj = 2* sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                                else:
                                    hi = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                                    hj =sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                                '''
                                hi = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                                hj =sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                                hihj+= (hi*hj)
                hihj = hihj/2
                #-------------------------------------------------
                try:
                    mod = (1.0/(2*E[l]))*(Aij - (hihj*1.0/(2*E[l])))
                except:
                    mod=0
                
                if f==1:# f is always 0
                    modc_layer+=mu*mod
                else:    
                    modc_layer+=mod
                    x1[c]+=mod
        #----------------------------------------------------------------------
            modc_couple=0
            Aij=0
            cicj=0
            
            #compute Aij------------------------------
            for n in commu[c]:
                if n in node_c:
                    for nei in node_c[n]:
                        if nei in commu[c]:
                            '''
                            if (n==nei):
                                Aij+= 2*graph[n][nei].get('weight',1)
                            else:
                                Aij+= graph[n][nei].get('weight',1)
                            '''
                            Aij+= graph[n][nei].get('weight',1)
            Aij = Aij/2

            #compute cicj-----------------------------
            for n1 in commu[c]:
                for n2 in commu[c]:
                    if(n1==n2 or nodelayer[n1]==nodelayer[n2] ) :
                        #If nodes in same layer, continue
                        continue
                    if n1 in node_c:
                        ci = sum([graph[n1][nbr].get('weight',1) for nbr in node_c[n1]])
                    else:
                        ci = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                    '''
                    if n1 in node_l: exp1
                        hi = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                        if n1 not in node_c:
                            ci=hi
                    '''
                    if n2 in node_c:
                        cj = sum([graph[n2][nbr].get('weight',1) for nbr in node_c[n2]])
                    else:
                        cj = sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                        
                    '''
                    if n2 in node_l: exp1
                        hj = sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                        if n2 not in node_c:
                            cj=hj
                    '''
                    cicj +=  (ci*cj)
                #cicj += ((ci+hi)*(cj+hj)) exp1
            cicj = cicj/2

            try:
                norm = 1.0/(2*sum([E[l] for l in layer]) + E12) #norm = 1/(2*|E1| + 2*|E2| + |E12|)
                mod = norm*(Aij - (norm*cicj))
            except:
                mod=0
            
            if f==1:
                modc_couple+=2*(1-mu)*mod
            else:
                modc_couple+=mod
                x2[c]+=mod
            #print modc_couple            
            ##print "ha hh"
            #-----------------------------------------------------------------------
            modularity+=modc_layer+modc_couple
            
        else:
            l=list(intra_inter[c])
            l=l[0];
            Aij=0.0
            hihj=0.0
            
            #compute summation AIJ------------------------
            for n in commu[c]:
                for nei in node_l.get(n,set()):
                    if nei in commu[c]:
                        '''
                        if(nei==n):
                            Aij+= 2*graph[n][nei].get('weight',1)
                        else:
                            Aij+=graph[n][nei].get('weight',1)
                        '''
                        Aij+=graph[n][nei].get('weight',1)

            Aij = Aij/2

            #compute summation hihj--------------------------
            for n1 in commu[c]:
                for n2 in commu[c]:
                    if(n1==n2):
                        if(n1 not in node_l[n1]):
                            continue                            
                    hi = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                    ci = sum([graph[n1][nbr].get('weight',1) for nbr in node_c.get(n1,set())])
                    hj =sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                    cj = sum([graph[n2][nbr].get('weight',1) for nbr in node_c.get(n2,set())])
                    '''
                    if(n1==n2):
                        hihj+= 2*(hi*hj)
                    else:
                        hihj+= (hi*hj)
                    '''
                    hihj+= (hi*hj)
                    #hihj+= ((hi+ci)*(hj+cj))

            hihj = hihj/2
            #-------------------------------------------------
            try:
                norm = 1.0/(2*E[l] + E12)
                mod = norm*(Aij*1.0 - (norm*hihj))
            except:
                mod=0
            
            if f==1: #f is always 0
                modc_layer+=mu*mod
            else:    
                modc_layer+=mod
                x1[c]+=mod

            modularity+=modc_layer        
                                    
    ##print x1,x2    
    return 0.333*modularity    

