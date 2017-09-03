def __modularity(commu, status, graph):

    layer=status.layer
    node_l=status.node_l
    node_c=status.node_c       
    top=status.top
    bot=status.bot
    edge_l=status.edge_l
    edge_c=status.edge_c
    couple=status.couple
    mu = status.mu

    f=0    
    intra_inter={}
    for c in commu:
        intra_inter[c]=set()
        
        for n in commu[c]:
            for l in layer:
                if n in layer[l]:
                    intra_inter[c].add(l)
    
    #print intra_inter
    
    modularity=0    
    x1={}
    x2={}    
    for c in commu:
        x1[c]=0
        x2[c]=0
        modc_layer=0
        
        if len(intra_inter[c])>1:
            for l in layer:
                d_layer=0.0
                I_layer=0.0
                m_layer=0.0
                n_layer=0.0
                n_co_com_layer=0.0
            
                for n in layer[l]:
                    if n in node_l:
                        m_layer+=len(node_l[n])
                    #if n in node_c:
                    #    m_layer+=len(node_c[n])    
                    if n in commu[c]:    #if the node belongs to current community
                        n_layer+=1
                        if n in node_l:
                            d_layer+=len(node_l[n])
                    
                            for nei in node_l[n]:
                                if nei in commu[c]: #if the neighbour belongs to current community
                                    I_layer+=1
                    
                        if n in node_c:
                            for nei in node_c[n]:
                                if nei not in commu[c]: #connected to atleast one crosslayer node outside community
                                    n_co_com_layer+=1 
                                    break
                        '''
                            for nei in node_c[n]:
                                if nei in commu[c]: #if the neighbour belongs to current community
                                    I_layer+=1
                        '''
                                            
                I_layer=float(I_layer)/2.0
                m_layer=float(m_layer)/2.0
                if I_layer>-1 and m_layer>0:
                    mod=((I_layer/m_layer)-((float(d_layer)/(2*m_layer))*(float(d_layer)/(2*m_layer))))
                else:
                    mod=0    
                #if n_layer > 0:
                #    mod=mod*pow(2.718,-(float(n_co_com_layer)/float(n_layer)))
            
                ##print mod    
                #print c,l,mod,mu,mu*mod,I_layer,m_layer,d_layer
                if f==1:
                    modc_layer+=mu*mod
                else:    
                    modc_layer+=mod
                    x1[c]+=mod
            
            modc_couple=0
            for co in couple:
                d_couple_top=0.0
                d_couple_bot=0.0
                I_couple=0.0
                m_couple=0.0
            
                top_tot=0.0
                top_con=0.0
                bot_tot=0.0
                bot_con=0.0
                
                m_layer_top=0.0
                m_layer_bot=0.0
                
                for n in couple[co]:
                    if n in node_c:
                        for nei in node_c[n]:
                            if (n in layer[top[co]] and nei in layer[bot[co]]) or (n in layer[bot[co]] and nei in layer[top[co]]):
                                m_couple+=1
                
                    if n in layer[top[co]]: #n belongs to the top layer of coupling
                        if n in node_l:
                            m_layer_top+=len(node_l[n])
                            
                        if n in commu[c]:    #if the node belongs to current community
                            top_tot+=1
                            if n in node_c:
                                flagg=0
                                for nei in node_c[n]:
                                    if nei in layer[bot[co]]:
                                        d_couple_top+=1
                                        if nei in commu[c]: #if the neighbour belongs to current community
                                            I_couple+=1
                                            flagg=1
                                if flagg==1: #connected to at least 1 within community node in bottom layer
                                    top_con+=1
                            if n in node_l and n not in node_c:
                                d_couple_top+=len(node_l[n])
                                    
                                            
                    if n in layer[bot[co]]: #n belongs to the bottom layer of coupling
                        if n in node_l:
                            m_layer_bot+=len(node_l[n])
                        if n in commu[c]:    #if the node belongs to current community
                            bot_tot+=1
                            if n in node_c:
                                flagg=0
                                for nei in node_c[n]:
                                    if nei in layer[top[co]]:
                                        d_couple_bot+=1    
                                        if nei in commu[c]: #if the neighbour belongs to current community
                                            flagg=1                    
                                                #break
                                if flagg==1: #connected to at least 1 within community node in bottom layer
                                    bot_con+=1
                            if n in node_l and n not in node_c:
                                d_couple_bot+=len(node_l[n])        
                                
                I_couple=float(I_couple)
                m_couple=float(m_couple)/2.0
                m_layer_top=float(m_layer_top)/2.0
                m_layer_bot=float(m_layer_bot)/2.0
                
                if I_couple>-1 and (m_couple+m_layer_top+m_layer_bot)>0:
                    
                    mod=((I_couple/(m_couple+m_layer_top+m_layer_bot))-((d_couple_top*d_couple_bot)/((m_couple+2.0*m_layer_top+2.0*m_layer_bot)*(m_couple+2.0*m_layer_top+2.0*m_layer_bot))))
                else:
                    mod=0    
                ##print I_couple, m_couple, d_couple_top, d_couple_bot    
                #print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
                if f==1:
                    modc_couple+=2*(1-mu)*mod
                else:
                    modc_couple+=mod
                    x2[c]+=mod
                #print modc_couple            
                ##print "ha hh"
            modularity+=modc_layer+modc_couple
            
        else:
            for l in layer:
                d_layer=0.0
                I_layer=0.0
                m_layer=0.0
                n_layer=0.0
                n_co_com_layer=0.0
            
                for n in layer[l]:
                    if n in node_l:
                        m_layer+=len(node_l[n])
                
                c_layer=0.0    
                for n in layer[l]:
                    if n in node_c:
                        c_layer+=len(node_c[n])
                            
                    #if n in node_c:
                    #    m_layer+=len(node_c[n])    
                    if n in commu[c]:    #if the node belongs to current community
                        n_layer+=1
                        if n in node_l:
                            d_layer+=len(node_l[n])
                    
                            for nei in node_l[n]:
                                if nei in commu[c]: #if the neighbour belongs to current community
                                    I_layer+=1
                    
                        if n in node_c and len(node_c[n])>0:
                            n_co_com_layer+=1 
                            #d_layer+=len(node_c[n])
                            
                        
                                            
                I_layer=float(I_layer)/2.0
                m_layer=float(m_layer)/2.0
                
                if I_layer>-1 and m_layer>0:#(1.0*m_layer+(c_layer/2.0))>0:# and m_layer1>0:
                    
                    #mod=((I_layer/(1.0*m_layer+(c_layer/2.0)))-((float(d_layer)/(2.0*m_layer+c_layer))*(float(d_layer)/(2.0*m_layer+c_layer))))
                    mod=((I_layer/(1.0*m_layer))-((float(d_layer)/(1.0*m_layer))*(float(d_layer)/(1.0*m_layer))))
                    #else:
                    #mod=(-((d_layer/(2*m_layer))*(d_layer/(2*m_layer))))    
                else:
                    mod=0 
                #if n_layer > 0:
                #    mod=mod*pow(2.718,-((float(n_co_com_layer))/float(n_layer)))
            
                ##print mod    
                #print c,l,mod,mu,mu*mod,I_layer,m_layer,c_layer,d_layer, n_co_com_layer, n_layer
                if f==1:
                    modc_layer+=mu*mod
                else:    
                    modc_layer+=mod
                    x1[c]+=mod
            modularity+=modc_layer        
                                    
    #print x1,x2    
    return 0.333*modularity    


def getModularityQ(commu, status):
    '''
    Not getting called from anywhere
    '''
    layer=status.layer
    node_l=status.node_l
    node_c=status.node_c       
    top=status.top
    bot=status.bot
    edge_l=status.edge_l
    edge_c=status.edge_c
    couple=status.couple
    mu = status.mu

    f=0
    modularity=0        
    for c in commu:
        modc_layer=0
        for l in layer:
            d_layer=0
            I_layer=0
            m_layer=0
            
            for n in layer[l]:
                if n in node_l:
                    m_layer+=len(node_l[n])
                if n in commu[c]:    #if the node belongs to current community
                    if n in node_l:
                        d_layer+=len(node_l[n])
                    
                        for nei in node_l[n]:
                            if nei in commu[c]: #if the neighbour belongs to current community
                                I_layer+=1
            I_layer=float(I_layer)/2.0
            m_layer=float(m_layer)/2.0
            if I_layer>-1:
                mod=((I_layer/m_layer)-((d_layer/(2*m_layer))*(d_layer/(2*m_layer))))
            else:
                mod=0   
            ##print c,l,mod,mu,mu*mod,I_layer,m_layer,d_layer
            if f==1:
                modc_layer+=mu*mod
            else:   
                modc_layer+=mod
        modc_couple=0
        for co in couple:
            d_couple_top=0
            d_couple_bot=0
            I_couple=0
            m_couple=0
            
            for n in couple[co]:
                if n in node_c:
                    for nei in node_c[n]:
                        if (n in layer[top[co]] and nei in layer[bot[co]]) or (n in layer[bot[co]] and nei in layer[top[co]]):
                            m_couple+=1
                
                if n in layer[top[co]]: #n belongs to the top layer of coupling
                    if n in commu[c]:    #if the node belongs to current community
                        if n in node_c:
                            for nei in node_c[n]:
                                if nei in layer[bot[co]]:
                                    d_couple_top+=1
                                    if nei in commu[c]: #if the neighbour belongs to current community
                                        I_couple+=1
                if n in layer[bot[co]]: #n belongs to the bottom layer of coupling
                    if n in commu[c]:    #if the node belongs to current community
                        if n in node_c:
                            for nei in node_c[n]:
                                if nei in layer[top[co]]:
                                    d_couple_bot+=1                     
                                    
            I_couple=float(I_couple)
            m_couple=float(m_couple)/2.0
            if I_couple>-1 and m_couple!=0:
                mod=((I_couple/m_couple)-((d_couple_top*d_couple_bot)/((m_couple)*(m_couple))))
            else:
                mod=0   
            ##print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
            if f==1:
                modc_couple+=2*(1-mu)*mod
            else:
                modc_couple+=mod        
        
        modularity+=modc_layer+modc_couple                  
    ##print len(layer)
    ##print len(couple)  
    return (1/float(len(layer)+len(couple)))*modularity 
