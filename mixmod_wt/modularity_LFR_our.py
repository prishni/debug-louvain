def getModularityQ(commu,status,graph):
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
	x1={}
	x2={}	
	for c in commu:
		x1[c]=0
		x2[c]=0
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
			if I_layer>-1 and m_layer>0:
				mod=((I_layer/m_layer)-((float(d_layer)/(2*m_layer))*(float(d_layer)/(2*m_layer))))
			else:
				mod=0	
			#print mod	
			print c,l,mod,mu,mu*mod,I_layer,m_layer,d_layer
			if f==1:
				modc_layer+=mu*mod
			else:	
				modc_layer+=mod
				x1[c]+=mod
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
			if I_couple>-1 and m_couple>0:
				mod=((I_couple/m_couple)-((float(d_couple_top)*float(d_couple_bot))/((m_couple)*(m_couple))))
			else:
				mod=0	
			#print I_couple, m_couple, d_couple_top, d_couple_bot	
			print c,co,mod,mu,2*(1-mu)*mod,I_couple,m_couple,d_couple_top,d_couple_bot
			if f==1:
				modc_couple+=2*(1-mu)*mod
			else:
				modc_couple+=mod
				x2[c]+=mod
			print modc_couple			
			#print "ha hh"
		modularity+=modc_layer+modc_couple					
	print x1,x2	
	return 0.333*modularity	

def getModularity_adapt(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c):
	#edge_l
	#print edge_l,edge_c
	edge_l={}
	for l in layer:
		s=0.0
		for n in layer[l]:
			if n in node_l: 
				s+=len(node_l[n])
		edge_l[l]=s/2.0	
				
	#edge_c
	edge_c={}
	for c in couple:
		s=0.0
		for n in layer[top[c]]:
			if n in node_c: 
				s+=len(node_c[n])
		for n in layer[bot[c]]:
			if n in node_c: 
				s+=len(node_c[n])		
		edge_c[c]=s/2.0
	
	print edge_l,edge_c
	modularity=0
	x1={}
	x2={}
	x3={}
	x4={}
	x5={}
	for c in commu:
		x1[c]=0
		x2[c]=0
		x3[c]=0
		x4[c]=0
		x5[c]=0
		for n1 in commu[c]:
			for n2 in commu[c]:
				if n1>n2:
					for l in layer:
						if n1 in layer[l] and n2 in layer[l]:
							aij=0
							d1=0
							d2=0
							if n1 in node_l:
								d1=len(node_l[n1])
								if n2 in node_l[n1]:
									aij=1
									
							if n2 in node_l:
								d2=len(node_l[n2])
								
							ss=(aij/((edge_l[l])))-((d1*d2)/((2*(edge_l[l]))*(2*(edge_l[l]))))
							x4[c]+=ss
							modularity+=ss
							
				 	for co in couple:
				 		if (n1 in layer[top[co]] and n2 in layer[bot[co]]) or (n2 in layer[top[co]] and n1 in layer[bot[co]]):
				 			aij=0
							d1=0
							d2=0
							if n1 in node_c:
								d1=len(node_c[n1])
								if n2 in node_c[n1]:
									aij=1
								
							if n2 in node_c:
								d2=len(node_c[n2])
							if edge_c[co]>0:	
								ss=(aij/((edge_c[co])))-((d1*d2)/((edge_c[co])*(edge_c[co])))
								x1[c]+=(aij/((edge_c[co])))
								
								x2[c]+=((d1*d2)/((edge_c[co])*(edge_c[co])))
									
							else:
								ss=0	
							modularity+=ss
							x5[c]+=ss
							#x2[c]+=((d1*d2)/((edge_c[co])*(edge_c[co])))
	print x1
	print x2
	print x4
	print x5					
	return 0.333*modularity	


def getSeries(filename):
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
	#modu=getModularity_adapt(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c)
	modu=getModularityQ(layer,couple,node_l,node_c,top,bot,commu,mu,edge_l,edge_c)
	'''
	max_modu=-2
	max_mu=-1
	count=0
	CC={}
	while mu<=1:		
		modu=getModularity(layer,couple,node_l,node_c,top,bot,commu,mu)
		fp1.write(str(mu)+" "+str(modu)+"\n")
		CC[count]=modu
		if modu>max_modu:
			max_modu=modu
			max_mu=mu
		mu=mu+0.01
		count+=1	
	#print max_modu,max_mu

	#return CC
	fp1.close()
	return [max_modu,max_mu]
	'''
	return modu,n_edge



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
...
'''
def main():
	fp1=open('./mQ/mQ_100_5_config5_no_ext.txt','w')
	for perturb in range(0,101):
		frac=float(perturb)/100.0
		str1='./config5_no_ext/config5_'+str(frac)+'.txt'
		print str1
		modu,couple=getSeries(str1)
		
		print frac,couple,modu
		fp1.write(str(frac)+' '+str(couple)+' '+str(modu)+'\n')
	fp1.close()
'''
p_list=[40,50,60,70,80,90]
#p_list=[1]
for p in p_list:
	ss='./mQ/single/mQ_single_single_1_'+str(p)+'.txt'	
	#ss='./mQ_modified_new_march1_corrected/mQ_modified_march1_100_5_config5_1_no_ext.txt'
	fp1=open(ss,'w')
	for perturb in range(0,101):
		frac=float(perturb)/100.0
		#str1='./config5_1_no_ext/config5_1_'+str(frac)+'.txt'
		str1='./new_single_layer/new_config_single/with_1_coupling/Test_config('+str(p)+'%)/test_config'+str(frac)+'.txt'
		#str1='./new_cross_layer/new_config_single/with_0.002_clique/Test_config('+str(p)+'%)/test_config'+str(frac)+'.txt'
		print str1
		modu,couple=getSeries(str1)
	
		print frac,couple,modu
		fp1.write(str(frac)+' '+str(couple)+' '+str(modu)+'\n')
	fp1.close()
	
p_list=[40,50,60,70,80,90]
#p_list=[1]
for p in p_list:
	ss='./mQ/single/mQ_single_multi_1_'+str(p)+'.txt'	
	#ss='./mQ_modified_new_march1_corrected/mQ_modified_march1_100_5_config5_1_no_ext.txt'
	fp1=open(ss,'w')
	for perturb in range(0,101):
		frac=float(perturb)/100.0
		#str1='./config5_1_no_ext/config5_1_'+str(frac)+'.txt'
		str1='./new_single_layer/new_config_multi/with_1_coupling/Test_config('+str(p)+'%)/test_config'+str(frac)+'.txt'
		#str1='./new_cross_layer/new_config_multi/with_0.002_clique/Test_config('+str(p)+'%)/test_config'+str(frac)+'.txt'
		print str1
		modu,couple=getSeries(str1)
	
		print frac,couple,modu
		fp1.write(str(frac)+' '+str(couple)+' '+str(modu)+'\n')
	fp1.close()	
'''	
if __name__ == '__main__':
	main()