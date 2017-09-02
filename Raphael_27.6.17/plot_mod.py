import matplotlib.pyplot as plt
import os
import sys
import matplotlib

files= os.listdir("./infos")
#x=range(1,len(files))
yval_gt=[]
yval_mm=[]
yval_l =[]
x=[]
i=1
for file in files:
	flag=0;
	file = "./infos/"+file
	fp =open(file)
	x.append(i)
	i=i+1
	for line in fp:
		if(line=="\n"): continue
		try:
			if(len(line.split(" "))==1 and len(line.split("-"))==1 and len(line.split("_"))==1):
				if(flag==0):
					flag=1
					print(file)
					print("GT: "+line.strip())
					yval_gt.append(float(line.strip()))
				elif(flag==1):
					flag=2
					#print(line)
					yval_l.append(float(line.strip()))
				elif(flag==2):
					flag=0
					print("MM: "+line.strip())
					print("\n")
					yval_mm.append(float(line.strip()))
					break
		except:
			pass
	fp.close()
#print(yval_gt)
#print("\n")
#print(yval_mm)

plt.figure()
matplotlib.style.use('ggplot')
#plt.subplot(211)
#plt.plot(x, yval_gt, x, yval_mm)
#plt.axis([0,1,0,1])
plt.plot(x,yval_gt,'g')
#plt.subplot(212)
plt.plot(x,yval_mm,'r')
plt.ylabel('Modularity')
plt.savefig("mod_plot.png")
plt.clf()

#plt.show()