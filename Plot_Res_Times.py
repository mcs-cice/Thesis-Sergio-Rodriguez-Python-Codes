'''
    Author: Sergio Rodriguez Peña
'''



#####################################################################
#                                                                   #
#   Code that plots the results of Residende_Time.py output file.   #
#                                                                   #
#####################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt



''' DEFINED FUNCTIONS '''

# Function that reads a residence times file.
def lect_distribution(file):
	dat=open(file+'.txt')
	datos=[]
	for lin in dat:
		l=lin.split()
		datos.append(l)
	ind=0
	index=[]
	times=[]
	distributions=[]
	for i in np.arange(len(datos)):
		if len(datos[i])==0:
			continue
		if datos[i][0]=='Time':
			ind+=1
			index.append(datos[i][1])
			times.append([])
			distributions.append([])
	for i in np.arange(len(datos)):
		if len(datos[i])==0:
			continue
		if datos[i][0]=='Time':
			for j in np.arange(4,len(datos[i])):
				times[index.index(datos[i][1])].append(float(datos[i][j]))
		if datos[i][0]=='#':
			for j in np.arange(3,len(datos[i])):
				distributions[index.index(datos[i][2])].append(float(datos[i][j]))		
	return(index,times,distributions)



''' CODE '''

# Read the files.
index=[]
time=[]
distribution=[]
sigue='y'
print('\n')
while sigue=='y':
	file=input('Name of the file: ')
	a,b,c=lect_distribution(file)
	for i in np.arange(len(a)):
		index.append(a[i])
		time.append(b[i])
		distribution.append(c[i])
	sigue=input('Read another file (y/n): ')
	print('\n')

# Compute the average residence time.
medias=[]
error=[]
for i in np.arange(len(distribution)):
	a=0
	for j in np.arange(len(distribution[i])):
	    a=a+distribution[i][j]*time[i][j]
	medias.append(a/np.sum(distribution[i]))
	b=0
	for j in np.arange(len(distribution[i])):
	    b=b+((time[i][j]-medias[i])**2*distribution[i][j])
	error.append(np.sqrt(b/np.sum(distribution[i])))

# Show the reuslts on screen.
print('\n')
for i in np.arange(len(medias)):
	print('Avg. Residence Time %s = %.3f +- %.3f ns' %(index[i],medias[i],error[i]))
print('\n')

# Print composition of the residence time files.
print('Composition of Systems: \n')
for i in np.arange(len(index)):
	print('%s  --  %i' %(index[i],i))
print('\n')
print('END  --  n')
print('\n')

# Ask for which to plot, with the corresponding label.
new=-1
plo=[]
labels=[]
while new!='n':
	new=input('Which one plot?: ')
	if new!='n':
		plo.append(int(new))
		l=input('How to label?: ')
		labels.append(l)

# Minimum residence time probability.
m=10
for i in np.arange(len(plo)):
	if min(distribution[plo[i]]/np.sum(distribution[plo[i]]))*100<m:
		m=min(distribution[plo[i]]/np.sum(distribution[plo[i]]))*100

# Maximum residence time.
n=0
for i in np.arange(len(plo)):
	if max(time[plo[i]])>n:
		n=max(time[plo[i]])

# Adjust of the decay region and width of the minimum probability.
slope=[]
error_slope=[]
oo=[]
error_oo=[]
width=[]
for i in np.arange(len(plo)):
	x_fit=[]
	y_fit=[]
	width.append([0,0])
	for j in np.arange(len(distribution[plo[i]])):
		x_fit.append(np.log(time[plo[i]][j]))
		y_fit.append(np.log(distribution[plo[i]][j]/np.sum(distribution[plo[i]])*100))
		if distribution[plo[i]][j]==1:
			width[i][0]=time[plo[i]][j]
			break
	for j in np.arange(1,len(distribution[plo[i]])):
		if distribution[plo[i]][-j]==1:
			width[i][1]=time[plo[i]][-j]
			break
	popt,pcov=np.polyfit(x_fit,y_fit,1,cov=True)
	slope.append(popt[0])
	oo.append(popt[1])
	error_slope.append(np.sqrt(pcov[0,0]))
	error_oo.append(np.sqrt(pcov[1,1]))
	print('\n')
	print('Slope of %s = %.3f +- %.3f' %(index[plo[i]],slope[i],error_slope[i]))
	print('Width of %s = %.1f - %.1f ns' %(index[plo[i]],width[i][0],width[i][-1]))

# Plot the results.
colors=['k','violet','b','m','gold']
markers=['o','s','X','P','s']
for i in np.arange(len(plo)):
	plt.figure()
	for j in np.arange(len(plo)):
		plt.scatter(1e5,1e5,color=colors[j],s=500,label=labels[j],marker=markers[j])
		plt.scatter(time[plo[j]],distribution[plo[j]]/np.sum(distribution[plo[j]])*100,color=colors[j],s=100,alpha=0.1,marker=markers[j])
	plt.scatter(time[plo[i]],distribution[plo[i]]/np.sum(distribution[plo[i]])*100,color=colors[i],s=100,marker=markers[i])
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Time / ns',fontsize=30)
	plt.ylabel('Probability %',fontsize=30)
	plt.xticks([1,10,100,200],['1','10','100','200'],fontsize=30)
	plt.yticks([1e-4,1e-3,1e-2,1e-1,1,10,100],['0.0001','0.001','0.01','0.1','1','10','100'],fontsize=30)
	if i==(len(plo)-1):
		plt.legend(loc='upper left',bbox_to_anchor=(0.0,1.18),ncol=len(plo),edgecolor='0.0',fontsize=36)
	plt.xlim(0.1,n+50)
	plt.ylim(m/2,101)
plt.show(block=False)

