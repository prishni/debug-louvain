import os
import sys



test_dirs = [d for d in os.listdir('./params/mu') if os.path.isdir(os.path.join('./params/mu', d))]
param_dirs = [d for d in os.listdir('./params') if os.path.isdir(os.path.join('./params', d))]


for param_dir in param_dirs:
	str1 = './params/' + param_dir +'/'
	for test_dir in test_dirs:
		str2 = str1 + test_dir +'/'
		dirs = [d for d in os.listdir(str2) if os.path.isdir(os.path.join(str2, d))]
		for dir1 in dirs:
			str3 = str2+dir1+'/'
			print(str3 + '\n')
			fin = open(str3+'new_format')
			fout1 = open(str3+'hetero.net','w')
			fout2 = open(str3+'.meta','w')
			fout3 = open(str3+'ground_truth_communities','w')
			lr =[]
			nr =[100,100]
			fin.readline()  
			fin.readline()
			lr.append(int(fin.readline().strip()))
			i = 0
			while(i<lr[0]):
				line = fin.readline().split()
				e1 = int(line[0])-1
				e2 = int(line[1])-1
				fout1.write(str(e1)+' '+str(e2)+'\n')
				i+=1
			
			fin.readline()
			lr.append(int(fin.readline().strip()))
			i=0
			while(i<lr[1]):
				line = fin.readline().split()
				e1 = int(line[0])-1
				e2 = int(line[1])-1
				fout1.write(str(e1)+' '+str(e2)+'\n')
				i+=1

			fin.readline()  
			fin.readline()
			
			lr.append(int(fin.readline().strip()))
			i=0
			while(i<lr[2]):
				line = fin.readline().split()
				e1 = int(line[0])-1
				e2 = int(line[1])-1
				fout1.write(str(e1)+' '+str(e2)+'\n')
				i+=1

			
			commu = int(fin.readline().strip())
			i=0
			while(i<commu):
				line = fin.readline().strip().split()
				fout3.write(str(int(line[0])-1))	
				for l in line[1:]:
					fout3.write(' '+str(int(l)-1))	
				fout3.write('\n')
				i+=1

			fout2.write('lr = ['+str(lr[0])+', '+str(lr[1])+', '+str(lr[2])+']\n' )
			fout2.write('nr = [100, 100]' )

			fout1.close()
			fout2.close()
			fout3.close()
