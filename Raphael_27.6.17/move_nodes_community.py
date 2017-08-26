import math
import copy

def getModularityQ(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c):
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
				d_layer=0
				I_layer=0
				m_layer=0
				n_layer=0
				n_co_com_layer=0
			
				for n in layer[l]:
					if n in node_l:
						m_layer+=len(node_l[n])
					#if n in node_c:
					#	m_layer+=len(node_c[n])	
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
				#	mod=mod*pow(2.718,-(float(n_co_com_layer)/float(n_layer)))
			
				#print mod	
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
				#print I_couple, m_couple, d_couple_top, d_couple_bot	
				#print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
				if f==1:
					modc_couple+=2*(1-mu)*mod
				else:
					modc_couple+=mod
					x2[c]+=mod
				#print modc_couple			
				#print "ha hh"
			modularity+=modc_layer+modc_couple
			
		else:
			for l in layer:
				d_layer=0
				I_layer=0
				m_layer=0
				n_layer=0
				n_co_com_layer=0
			
				for n in layer[l]:
					if n in node_l:
						m_layer+=len(node_l[n])
				
				c_layer=0.0	
				for n in layer[l]:
					if n in node_c:
						c_layer+=len(node_c[n])
							
					#if n in node_c:
					#	m_layer+=len(node_c[n])	
					if n in commu[c]:    #if the node belongs to current community
						n_layer+=1
						if n in node_l:
							d_layer+=len(node_l[n])
					
							for nei in node_l[n]:
								if nei in commu[c]: #if the neighbour belongs to current community
									I_layer+=1
					
						if n in node_c and len(node_c[n])>0:
							n_co_com_layer+=1 
							d_layer+=len(node_c[n])
							
						
											
				I_layer=float(I_layer)/2.0
				m_layer=float(m_layer)/2.0
				
				if I_layer>-1 and (1.0*m_layer+(c_layer/2.0))>0:# and m_layer1>0:
					
					mod=((I_layer/(1.0*m_layer+(c_layer/2.0)))-((float(d_layer)/(2.0*m_layer+c_layer))*(float(d_layer)/(2.0*m_layer+c_layer))))
					#else:
					#mod=(-((d_layer/(2*m_layer))*(d_layer/(2*m_layer))))    
				else:
					mod=0 
				#if n_layer > 0:
				#	mod=mod*pow(2.718,-((float(n_co_com_layer))/float(n_layer)))
			
				#print mod	
				#print c,l,mod,mu,mu*mod,I_layer,m_layer,c_layer,d_layer, n_co_com_layer, n_layer
				if f==1:
					modc_layer+=mu*mod
				else:	
					modc_layer+=mod
					x1[c]+=mod
			modularity+=modc_layer		
									
	#print x1,x2	
	return 0.333*modularity	



def getSeries(filename,commu1):
	fp=open(filename,'r')
	#fp1=open(filename+"mu_modu.txt",'w')
	line=fp.readline()
	line=line.rstrip()
	n_layer=int(line)
	layer={}
	node_l={}
	l_ID=1
	edge_l={}
	edge_c={}
	for i in range(0,n_layer):
		line=fp.readline()
		line=line.rstrip()
		line=line.split()
		layer[l_ID]=set()
		#print line
		for n in line:
			layer[l_ID].add(int(n))
		line=fp.readline()
		line=int(line.rstrip())
		n_edge=line
		edge_l[l_ID]=n_edge
		for j in range(0,n_edge):
			line=fp.readline()
			line=line.rstrip()
			line=line.split()
			n1=int(line[0])
			n2=int(line[1])	
			if n1 not in node_l:
				node_l[n1]=set()
			node_l[n1].add(n2)		
			if n2 not in node_l:
				node_l[n2]=set()
			node_l[n2].add(n1)
		l_ID+=1
		
	line=fp.readline()
	line=line.rstrip()
	n_couple=int(line)
	node_c={}		
	top={}
	bot={}
	c_ID=1
	couple={}

	for i in range(0,n_couple):
		line=fp.readline()
		line=line.rstrip()
		line=line.split()
	
		top[c_ID]=int(line[0])
		bot[c_ID]=int(line[1])
		couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
		
		line=fp.readline()
		line=int(line.rstrip())
		n_edge=line
		edge_c[c_ID]=n_edge
		for j in range(0,n_edge):
			line=fp.readline()
			line=line.rstrip()
			line=line.split()
			n1=int(line[0])
			n2=int(line[1])	
			if n1 not in node_c:
				node_c[n1]=set()
			node_c[n1].add(n2)		
			if n2 not in node_c:
				node_c[n2]=set()
			node_c[n2].add(n1)	
		c_ID=c_ID+1

	line=fp.readline()
	line=line.rstrip()
	n_comm=int(line)
	commu={}
	com_ID=1
	for i in range(0,n_comm):
		line=fp.readline()
		line=line.rstrip()
		line=line.split()
		commu[com_ID]=set()
		for n in line:
			commu[com_ID].add(int(n))
		com_ID+=1		

	#print len(layer),len(node_l),len(node_c),top,bot,len(commu),len(couple)
	mu=0
	modu=getModularityQ(layer,couple,node_l,node_c,top,bot,commu1,mu,edge_l,edge_c)
	return modu



'''
FORMAT OF INPUT FILE
***************************

no_layers
layer1 vertices
no_layer1_edges
layer1_edge1
layer1_edge2
...
layer2 vertices
no_layer2_edges
layer2_edge1
layer2_edge2
...
no_couplings
coupling1_top_layer
coupling1_bot_layer
no_coupling1_edges
coupling1_edge1
coupling1_edge2
...
coupling2_top_layer
coupling2_bot_layer
no_coupling2_edges
coupling2_edge1
coupling2_edge2
...
no_communities
community1_vertices
community2_vertices
'''

a_list=[0.1]
p_list=[0.9]
#a_list=[0.1,0.2,0.4,0.5,0.6,0.8,0.9]
#p_list=[0.1,0.3,0.5,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
mu=0.05
p1=1.0
p2=0.0

for p in p_list:
	for a in a_list:
		ss1='./infos/network_'+str(a)+'_'+str(p)+'_'+str(mu)+'_'+str(p1)+'_'+str(p2)+'.info'
		ss1 = "network_0.9_1.0_0.05_1.0_0.0.info"
		print ss1
		fp1=open(ss1,'r')
		commu1={}
		commu2={}
		commu3={}
		flag=-1
		for line in fp1:
			if flag==0:
				line1=(line.strip()).split()
			
				if len(line1)==1:
					flag=1
				
				else:
					if len(line1)>=1:
						line2=line.strip()
						line2=((line2.lstrip(str(id1)+" set(")).lstrip('[')).rstrip("])").split(', ')
						#print line2
						commu1[id1]=set()
						for node in line2:
							node1=int(node)
							commu1[id1].add(node1)
						id1+=1
						
			if flag==1:
				line1=(line.strip()).split()
			
				if len(line1)==1:
					flag=2
				
				else:
					if len(line1)>=1:
						line2=line.strip()
						#print (line2.lstrip(str(id1-1)+" ")).lstrip('[')
						line2=((line2.lstrip(str(id1-1)+" ")).lstrip('[')).rstrip("]").split(', ')
						
						commu2[id1]=set()
						for node in line2:
							node1=int(node)
							commu2[id1].add(node1)
						id1+=1			
						
			if flag==2:
				line1=(line.strip()).split()
			
				if len(line1)==1:
					flag=3
				
				else:
					if len(line1)>=1:
						line2=line.strip()
						line2=((line2.lstrip(str(id1-1)+" ")).lstrip('[')).rstrip("]").split(', ')
						#print line2
						commu3[id1]=set()
						for node in line2:
							node1=int(node)
							commu3[id1].add(node1)
						id1+=1				
										
			if line.strip()=="Ground-truth":
				flag=0
				id1=1
			if line.strip()=="_commu_benching_all_march21_louvain.pickle":
				flag=1
				id1=1
			if line.strip()=="_commu_benching_all_march21_louvain_mixmod.pickle":
				flag=2
				id1=1
						
			
		fp1.close()
		
		
		
		ss='./synthetics/network_'+str(a)+'_'+str(p)+'_'+str(mu)+'_'+str(p1)+'_'+str(p2)
		ss = './synthetics/network_0.9_1.0_0.05_1.0_0.0'
		modu1=getSeries(ss,commu1)
		modu2=getSeries(ss,commu2)
		modu3=getSeries(ss,commu3)
		print modu1,modu2,modu3
		
		commu={}
		choice=input('Press 0 for Normal Louvain, 1 for Mixmod Louvain: ')
		while True:			
			if choice==0:
				commu=copy.deepcopy(commu2)
				modu=modu2
				break
			else:
			 	if choice==1:
			 		commu=copy.deepcopy(commu3)
			 		modu=modu3
			 		break
			 	else:
			 		choice=input('Please enter correctly;\nPress 0 for Normal Louvain, 1 for Mixmod Louvain: ')		
		if choice==0:
			print "Comparing Normal Louvain with Ground Truth===================>"	 		
		else:	
			print "Comparing MixMod Louvain with Ground Truth===================>"
		
		print "\n=====> Ground Truth ====> No. Communities, modularity ",len(commu1),modu1
		for c in commu1:
			str1="GT"+str(c)+"->  Length "+str(len(commu1[c]))+"->	"+str(commu1[c])
			print str1	
		
		print "\n=====> Detected ====> No. Communities, modularity ",len(commu),modu
		for c in commu:	
			max1=-1
			cmax=-1
			for c1 in commu1:
				intersec = len(commu1[c1].intersection(commu[c]))
				if intersec > max1:
					max1=intersec
					cmax=c1
			str1="DT"+str(c)+"->  Length "+str(len(commu[c]))+"->	"+str(commu[c])+" ---->Max Overlap GT"+str(cmax)+", "+str(max1)+"/"+str(len(commu[c])) 
			print str1		
			
		while True:			
			node=input('Please enter the node to move: ')
			tar=input('Pleas enter the target community: ')
			if node<1 or node>200:
				print "Please Enter correct node" 
			else:
				if tar<1 or tar>len(commu):
					print "Please Enter correct community"
				else:
					for c in commu:
						if node in commu[c]:
							commu[c].remove(node)
					commu[tar].add(node)
					modu11=getSeries(ss,commu)
					print 'Previous Modu, Current Modu', modu,modu11
					
					
					print "\n=====> Detected ====> No. Communities, modularity ",len(commu),modu11
					for c in commu:	
						max1=-1
						cmax=-1
						for c1 in commu1:
							intersec = len(commu1[c1].intersection(commu[c]))
							if intersec > max1:
								max1=intersec
								cmax=c1
						str1="DT"+str(c)+"->  Length "+str(len(commu[c]))+"->	"+str(commu[c])+" ---->Max Overlap GT"+str(cmax)+", "+str(max1)+"/"+str(len(commu[c])) 
						print str1		
					
					choice2=input('Do you want to stop (0) or continue (1)?: ')
					if choice2==0:
						break
							
					choice1=input('Do you want to keep the last movement (0) or revert to original detected (1)?: ')
					if choice1==0:
						modu=modu11
					else:
						commu={}
						
						if choice==0:
							commu=copy.deepcopy(commu2)
							modu=modu2
						else:
							commu=copy.deepcopy(commu3)
							modu=modu3
							
						print "\n=====> Detected ====> No. Communities, modularity ",len(commu),modu
						for c in commu:	
							max1=-1
							cmax=-1
							for c1 in commu1:
								intersec = len(commu1[c1].intersection(commu[c]))
								if intersec > max1:
									max1=intersec
									cmax=c1
							str1="DT"+str(c)+"->  Length "+str(len(commu[c]))+"->	"+str(commu[c])+" ---->Max Overlap GT"+str(cmax)+", "+str(max1)+"/"+str(len(commu[c])) 
							print str1			 
