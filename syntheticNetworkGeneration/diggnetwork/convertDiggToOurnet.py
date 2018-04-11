import sys,os

def main():
	edgelist1 =[]
	edgelist2 =[]
	coupling =[]
	filename = "DiggNewFormat"
	hetero = open("hetero.net",'r')
	for line in hetero:
		line = line.strip().split()
		line[0] = int(line[0])
		line[1] = int(line[1])
		if(line[0]>=60):
			line[0]+=40
		if(line[1]>=60):
			line[1]+=40
		if(line[0]<60 and line[1]<60):
			edgelist1.append((line[0],line[1]))
		elif(line[0]>=60 and line[1]>=60):
			edgelist2.append((line[0],line[1]))
		else:
			coupling.append((line[0],line[1]))
	hetero.close()

	f = open(filename,'w')
	f.write('2\n')
	for i in range(1,61):
		f.write(str(i) + " ")
	f.write("\n")
	leng = len(edgelist1)
	f.write(str(leng)+"\n")
	for each in edgelist1:
		f.write(str(int(each[0])+1) +" "+ str(int(each[1])+1)+ "\n")
	for i in range(101,161):
		f.write(str(i) + " ")
	f.write("\n")
	leng = len(edgelist2)
	f.write(str(leng)+"\n")
	for each in edgelist2:
		f.write(str(int(each[0])+1) +" "+ str(int(each[1])+1)+ "\n")		
	f.write("1\n1 2\n")
	f.write(str(len(coupling))+"\n")
	for each in coupling:
		f.write(str(int(each[0])+1) +" "+ str(int(each[1])+1)+ "\n")
	f.write('7\n')
	truth = open('truth','r')
	for line in truth:
		line = line.split()
		for each in line:
			each = int(each)
			if(each>=60):
				each+=40
			f.write(str(each+1)+" ")
		f.write('\n')
	truth.close()
	f.close()

if __name__ == '__main__':
	main()