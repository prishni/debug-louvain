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
    in_layer_in_comm = status.in_layer_in_comm
    in_layer_out_comm = status.in_layer_out_comm
    out_layer_in_comm = status.out_layer_in_comm
    out_layer_out_comm = status.out_layer_out_comm
   
    modularity=0    

    #compute total edges-------------------------------
    E=0
    Eloop=0
    for n in node_l:
            for nei in node_l.get(n,set()):
                #print(n," - ",nei," weight: ",graph[n][nei].get('weight',1))
                if n== nei:
                    Eloop+= graph[n][nei].get('weight',1)
                else:
                    E += graph[n][nei].get('weight',1)
    E=E/2
    E= E+Eloop

    #print("node_l: ",node_l)
    #print("node_c: ",node_c)
    #print("E: ",E)
    for c in commu:
        Aij=0.0
        aijloop=0.0
        didj=0.0

        for n in commu[c]:
            for nei in node_l.get(n,set()):
                if nei==n:
                    aijloop+=graph[n][nei].get('weight',1)
                else:
                    if nei in commu[c]:
                        Aij+=graph[n][nei].get('weight',1)
            Aij = Aij/2
            Aij+=aijloop

        #compute summation didj--------------------------
        for n1 in commu[c]:
            for n2 in commu[c]:
                if(n1==n2):
                    if(n1 not in node_l[n1]):
                        continue
                di = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                dj =sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                didj+= ((di)*(dj))

        didj = didj/2

        modularity+=(1.0/(2*E))*(Aij-(didj/(2*E)))
                                    
    ##print x1,x2    
    return modularity    

