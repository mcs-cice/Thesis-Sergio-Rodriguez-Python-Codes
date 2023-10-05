'''
    Author: Sergio Rodriguez Peña
'''



###################################################################################
#                                                                                 #
#   Code that plots the clusters created with the code Cluster_File_Creation.py   #
#                                                                                 #
###################################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors



''' DEFINED FUNCTIONS '''

# Function that checks the number of isolated anions and cations in each timestep.
def check_isolated(cl):
	iso_cation=0
	for i in np.arange(len(cl)):
		if len(cl[i])==1:
			iso_cation+=1
	a=0
	for i in np.arange(len(cl)):
		a=a+len(cl[i])
	iso_anion=len(index_cation)+len(index_anion)-a
	return(iso_anion,iso_cation)


# Function that creates the progress bar shown in the terminal.
def Progress_Bar(progress,total):
    percent=100*(progress/float(total))
    bar='|'+'#'*int(percent/2)+'-'*(50-(int(percent/2)))+'|'
    print(f'%s %.1f' %(bar,percent) +' %', end='\r')

# Function that plots the cluster map.
def Plot_Cluster_Map(C):
	if not big:
		X,Y=np.meshgrid(np.arange(len(C)),np.arange(len(C)))
		Z=np.zeros((len(C),len(C)))
		for i in np.arange(len(C)):
			for j in np.arange(len(C)):
				Z[i][j]=100*C[i][j]/np.sum(C)
		minimo=10
		for i in np.arange(len(Z)):
			for j in np.arange(len(Z)):
				if Z[i][j]==0:
					continue
				else:
					if Z[i][j]<minimo:
						minimo=Z[i][j]
		fig,ax=plt.subplots()
		cl=ax.pcolormesh(X,Y,Z,norm=colors.LogNorm(vmin=1e-2,vmax=100),shading='auto')
		cbar = plt.colorbar(cl,extend='max')
		for t in cbar.ax.get_yticklabels():
		     t.set_fontsize(22)
		plt.plot(np.arange(len(C)),np.arange(len(C)),linewidth=6,color='k')
		plt.xlabel('# Anions',fontsize=30)
		plt.ylabel('# Cations',fontsize=30)
		plt.xticks(ticks=np.arange(len(C)+1),fontsize=16)
		plt.yticks(ticks=np.arange(len(C)+1),fontsize=16)
		plt.grid(visible=True,color='gray')
		plt.xlim(-0.5,len(C)+0.5)
		plt.ylim(-0.5,len(C)+0.5)
	if big:
		X,Y=np.meshgrid(np.arange(50),np.arange(50))
		Z=np.zeros((50,50))
		for i in np.arange(len(C)-50,len(C)):
			for j in np.arange(len(C)-50,len(C)):
				Z[i-(len(C)-50)][j-(len(C)-50)]=100*C[i][j]/np.sum(C)
		fig,ax=plt.subplots()
		cl=ax.pcolormesh(X,Y,Z,norm=colors.LogNorm(vmin=1e-2,vmax=100),shading='auto')
		cbar = plt.colorbar(cl,extend='max')
		for t in cbar.ax.get_yticklabels():
		     t.set_fontsize(22)
		plt.plot(np.arange(len(C)),np.arange(len(C)),linewidth=6,color='k')
		plt.xlabel('# Anions',fontsize=30)
		plt.ylabel('# Cations',fontsize=30)
		plt.grid(visible=True,color='gray')
		plt.xlim(-0.5,50.5)
		plt.ylim(-0.5,50.5)	
		xt=[]
		for i in np.arange(51):
			xt.append('%s' %(len(C)-50+i))
		plt.xticks(np.arange(len(Z)+1),xt,fontsize=16)
		plt.yticks(np.arange(len(Z)+1),xt,fontsize=16)
	plt.show(block=False)	


# Function that plots the probability of the clusters.
def Plot_Cluster_Probability(C):
	# [positive,neutral,negative]
	isolated=[0,0,0]
	pairs=[0,0,0]
	small=[0,0,0]  # 3-8
	medium=[0,0,0]  # 9-16
	large=[0,0,0]  # >16
	supercluster=[0,0,0]
	if not big:
		for i in np.arange(len(C)):
			for j in np.arange(len(C[i])):
				if C[i][j]>0:
					q=Check_Charge(i,j)
					if (i+j)==1:
						isolated[q]=isolated[q]+C[i][j]
					if (i+j)==2:
						pairs[q]=pairs[q]+C[i][j]
					if (i+j)>3 and (i+j)<=8:
						small[q]=small[q]+C[i][j]
					if (i+j)>8 and (i+j)<=16:
						medium[q]=medium[q]+C[i][j]
					if (i+j)>16:
						large[q]=large[q]+C[i][j]
		total=np.sum(isolated)+np.sum(pairs)+np.sum(small)+np.sum(medium)+np.sum(large)
		plt.figure()
		positive=[100*isolated[0]/total,100*pairs[0]/total,100*small[0]/total,100*medium[0]/total,100*large[0]/total]
		plt.bar(np.arange(5)-0.2,positive,color='r',width=0.15,label='positive')
		neutral=[100*isolated[1]/total,100*pairs[1]/total,100*small[1]/total,100*medium[1]/total,100*large[1]/total]
		plt.bar(np.arange(5),neutral,color='k',width=0.15,label='neutral')
		negative=[100*isolated[2]/total,100*pairs[2]/total,100*small[2]/total,100*medium[2]/total,100*large[2]/total]
		plt.bar(np.arange(5)+0.2,negative,color='b',width=0.15,label='negative')
		plt.ylabel('Probability %',fontsize=20)
		plt.xticks(np.arange(5),['Iso','Pairs','S-Clust \n (3-8 ions)','M-Clust \n (9-16 ions)','L-Clust \n (>16 ions)'],fontsize=20)
		plt.yticks(fontsize=20)
		plt.legend(fontsize=30,loc='upper left',bbox_to_anchor=(0.0,1.15),ncol=3,edgecolor='0.0')
		plt.xlim(-0.5,4.5)
	isolated=[0,0,0]
	pairs=[0,0,0]
	small=[0,0,0]  # 3-8
	medium=[0,0,0]  # 9-16
	large=[0,0,0]  # >16
	supercluster=[0,0,0]
	if not big:
		for i in np.arange(len(C)):
			for j in np.arange(len(C[i])):
				if C[i][j]>0:
					q=Check_Charge(i,j)
					if (i+j)==1:
						isolated[q]=isolated[q]+C[i][j]*(i+j)
					if (i+j)==2:
						pairs[q]=pairs[q]+C[i][j]*(i+j)
					if (i+j)>3 and (i+j)<=8:
						small[q]=small[q]+C[i][j]*(i+j)
					if (i+j)>8 and (i+j)<=16:
						medium[q]=medium[q]+C[i][j]*(i+j)
					if (i+j)>16:
						large[q]=large[q]+C[i][j]*(i+j)
		total=np.sum(isolated)+np.sum(pairs)+np.sum(small)+np.sum(medium)+np.sum(large)
		plt.figure()
		positive=[100*isolated[0]/total,100*pairs[0]/total,100*small[0]/total,100*medium[0]/total,100*large[0]/total]
		plt.bar(np.arange(5)-0.2,positive,color='r',width=0.15,label='positive')
		neutral=[100*isolated[1]/total,100*pairs[1]/total,100*small[1]/total,100*medium[1]/total,100*large[1]/total]
		plt.bar(np.arange(5),neutral,color='k',width=0.15,label='neutral')
		negative=[100*isolated[2]/total,100*pairs[2]/total,100*small[2]/total,100*medium[2]/total,100*large[2]/total]
		plt.bar(np.arange(5)+0.2,negative,color='b',width=0.15,label='negative')
		plt.ylabel('Probability %',fontsize=20)
		plt.xticks(np.arange(5),['Iso','Pairs','S-Clust \n (3-8 ions)','M-Clust \n (9-16 ions)','L-Clust \n (>16 ions)'],fontsize=20)
		plt.yticks(fontsize=20)
		plt.legend(fontsize=30,loc='upper left',bbox_to_anchor=(0.0,1.15),ncol=3,edgecolor='0.0')
		plt.xlim(-0.5,4.5)
	isolated=[0,0,0]
	pairs=[0,0,0]
	small=[0,0,0]  # 3-8
	medium=[0,0,0]  # 9-16
	large=[0,0,0]  # >16
	supercluster=[0,0,0]
	if big:
		for i in np.arange(len(C)):
			for j in np.arange(len(C[i])):
				if C[i][j]>0:
					q=Check_Charge(i,j)
					if (i+j)==1:
						isolated[q]=isolated[q]+C[i][j]
					if (i+j)==2:
						pairs[q]=pairs[q]+C[i][j]
					if (i+j)>3 and (i+j)<=8:
						small[q]=small[q]+C[i][j]
					if (i+j)>8 and (i+j)<=16:
						medium[q]=medium[q]+C[i][j]
					if (i+j)>16 and (i+j)<100:
						large[q]=large[q]+C[i][j]
					if (i+j)>100:
						supercluster[q]=supercluster[q]+C[i][j]
		total=np.sum(isolated)+np.sum(pairs)+np.sum(small)+np.sum(medium)+np.sum(large)+np.sum(supercluster)
		fig,(ax1,ax2)=plt.subplots(1,2,gridspec_kw={'width_ratios':[3, 1]})
		positive=[100*isolated[0]/total,100*pairs[0]/total,100*small[0]/total,100*medium[0]/total,100*large[0]/total,100*supercluster[0]/total]
		ax1.bar(np.arange(5)-0.2,positive[:5],color='r',width=0.15,label='positive')
		ax2.bar(-0.2,positive[5],color='r',width=0.15)
		neutral=[100*isolated[1]/total,100*pairs[1]/total,100*small[1]/total,100*medium[1]/total,100*large[1]/total,100*supercluster[1]/total]
		ax1.bar(np.arange(5),neutral[:5],color='k',width=0.15,label='neutral')
		ax2.bar(0,neutral[5],color='k',width=0.15)
		negative=[100*isolated[2]/total,100*pairs[2]/total,100*small[2]/total,100*medium[2]/total,100*large[2]/total,100*supercluster[2]/total]
		ax1.bar(np.arange(5)+0.2,negative[:5],color='b',width=0.15,label='negative')
		ax2.bar(0.2,negative[5],color='b',width=0.15)
		ax1.set_ylabel('Probability %',fontsize=20)
		ax2.yaxis.tick_right()
		ax2.yaxis.set_label_position("right")
		ax2.set_ylabel('Probability %',fontsize=20)
		ax1.set_xticks(np.arange(5))
		ax1.set_xticklabels(['Iso','Pairs','S-Clust \n (3-8 ions)','M-Clust \n (8-16 ions)','L-Clust \n (16-100 ions)'],fontsize=20)
		ax2.set_xticks([0])
		ax2.set_xticklabels(['Super-Clust'],fontsize=20)
		ax1.tick_params(axis='y', labelsize=20)
		ax2.tick_params(axis='y', labelsize=20)
		ax1.set_xlim(-0.5,4.5)
		ax2.set_xlim(-0.5,0.5)
		ax2.set_ylim(0,105)
		ax1.legend(fontsize=30,loc='upper left',bbox_to_anchor=(0.0,1.15),ncol=3,edgecolor='0.0')
	isolated=[0,0,0]
	pairs=[0,0,0]
	small=[0,0,0]  # 3-8
	medium=[0,0,0]  # 9-16
	large=[0,0,0]  # >16
	supercluster=[0,0,0]	
	if big:
		for i in np.arange(len(C)):
			for j in np.arange(len(C[i])):
				if C[i][j]>0:
					q=Check_Charge(i,j)
					if (i+j)==1:
						isolated[q]=isolated[q]+C[i][j]*(i+j)
					if (i+j)==2:
						pairs[q]=pairs[q]+C[i][j]*(i+j)
					if (i+j)>3 and (i+j)<=8:
						small[q]=small[q]+C[i][j]*(i+j)
					if (i+j)>8 and (i+j)<=16:
						medium[q]=medium[q]+C[i][j]*(i+j)
					if (i+j)>16 and (i+j)<100:
						large[q]=large[q]+C[i][j]*(i+j)
					if (i+j)>100:
						supercluster[q]=supercluster[q]+C[i][j]*(i+j)
		total=np.sum(isolated)+np.sum(pairs)+np.sum(small)+np.sum(medium)+np.sum(large)+np.sum(supercluster)
		fig,(ax1,ax2)=plt.subplots(1,2,gridspec_kw={'width_ratios':[3, 1]})
		positive=[100*isolated[0]/total,100*pairs[0]/total,100*small[0]/total,100*medium[0]/total,100*large[0]/total,100*supercluster[0]/total]
		ax1.bar(np.arange(5)-0.2,positive[:5],color='r',width=0.15,label='positive')
		ax2.bar(-0.2,positive[5],color='r',width=0.15)
		neutral=[100*isolated[1]/total,100*pairs[1]/total,100*small[1]/total,100*medium[1]/total,100*large[1]/total,100*supercluster[1]/total]
		ax1.bar(np.arange(5),neutral[:5],color='k',width=0.15,label='neutral')
		ax2.bar(0,neutral[5],color='k',width=0.15)
		negative=[100*isolated[2]/total,100*pairs[2]/total,100*small[2]/total,100*medium[2]/total,100*large[2]/total,100*supercluster[2]/total]
		ax1.bar(np.arange(5)+0.2,negative[:5],color='b',width=0.15,label='negative')
		ax2.bar(0.2,negative[5],color='b',width=0.15)
		ax1.set_ylabel('Probability %',fontsize=20)
		ax2.yaxis.tick_right()
		ax2.yaxis.set_label_position("right")
		ax2.set_ylabel('Probability %',fontsize=20)
		ax1.set_xticks(np.arange(5))
		ax1.set_xticklabels(['Iso','Pairs','S-Clust \n (3-8 ions)','M-Clust \n (8-16 ions)','L-Clust \n (16-100 ions)'],fontsize=20)
		ax2.set_xticks([0])
		ax2.set_xticklabels(['Super-Clust'],fontsize=20)
		ax1.tick_params(axis='y', labelsize=20)
		ax2.tick_params(axis='y', labelsize=20)
		ax1.set_xlim(-0.5,4.5)
		ax2.set_xlim(-0.5,0.5)
		ax2.set_ylim(0,105)
		ax1.legend(fontsize=30,loc='upper left',bbox_to_anchor=(0.0,1.15),ncol=3,edgecolor='0.0')
	plt.show(block=False)

# Function that checks the net charge of a cluster.
def Check_Charge(n_c,n_a):
	if n_c>n_a:
		q=0
	if n_c==n_a:
		q=1
	if n_c<n_a:
		q=2
	return(q)



''' CODE '''

# Name of the input file.
entrada=input('Cluster File Name: ')+'.txt'
dat=open(entrada)
datos=[]
for lin in dat:
	l=lin.split()
	datos.append(l)

# Ask for checking of the process by looking at the number of ions
# and clusters at each step.
print('\n')
check=input('Make checking (y/n): ')
clust=[]
t=[]
for i in np.arange(len(datos)):
	if len(datos[i])==0 or datos[i][0]=='Anion:' or datos[i][0]=='Cation:':
		continue
	elif datos[i][0]=='t':
		clust.append([])
		t.append(datos[i][2])
	else:
		clust[len(t)-1].append(datos[i])

# Index of the ions.
index_anion=np.arange(int(datos[0][1]),int(datos[0][3])+1)
index_cation=np.arange(int(datos[1][1]),int(datos[1][3])+1)	

# First checking by reading the file.
if check=='y':
	a=0
	b=0
	for i in np.arange(len(clust)):
		a=a+check_isolated(clust[i])[0]
		b=b+check_isolated(clust[i])[0]
		b=b+len(clust[i])
		for j in np.arange(len(clust[i])):
			a=a+len(clust[i][j])
	print('Fisrt Checking: \n')
	print('Total # Ions = %i' %(a))
	print('Total # Clusters = %i' %(b))
	print('\n')

# Structure of the clusters in matrix style.
'''
                  Columns = # Anions 
Files = # Cations                              
'''

# Transform the clusters in matrix style.
clust_matrix=[]
max_n=0
print('MEASURING CLUSTERS:')
if check=='y':
	c=0
	d=0
for i in np.arange(len(clust)):
	big=False
	for j in np.arange(len(clust[i])):
		if len(clust[i][j])>100:
			clust_matrix.append(np.zeros((len(index_cation)+1,len(index_cation)+1)))
			big=True
	if not big:
		clust_matrix.append(np.zeros((100,100)))
	iso_anion,iso_cation=check_isolated(clust[i])
	clust_matrix[i][0][1]=iso_anion
	clust_matrix[i][1][0]=iso_cation
	if check=='y':
		c=c+iso_cation+iso_anion
		d=d+iso_cation+iso_anion
	for j in np.arange(len(clust[i])):
		if len(clust[i][j])==1:
			continue
		if len(clust[i][j])>1:
			n_l=0
			n_a=0
			for k in np.arange(len(clust[i][j])):
				if int(clust[i][j][k]) in index_cation:
					n_l+=1
				if int(clust[i][j][k]) in index_anion:
					n_a+=1
		while n_l>len(clust_matrix[i])-1 or n_a>len(clust_matrix[i][0])-1:
			clust_matrix[i]=np.insert(clust_matrix[i],-1,np.zeros((len(clust_matrix[i][0]))),0)
			clust_matrix[i]=np.insert(clust_matrix[i],-1,np.zeros((len(clust_matrix[i]))),1)
		clust_matrix[i][n_l][n_a]+=1
		if check=='y':
			c=c+n_l+n_a
			d=d+1
		if n_l>max_n:
			max_n=n_l
		if n_a>max_n:
			max_n=n_a
	Progress_Bar(i,len(clust))
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %')
if check=='y':
	print('Second Checking: \n')
	print('Total # Ions = %i' %(c))
	print('Total # Clusters = %i' %(d))
	print('\n')

# Sum all the timesteps.
print('SUMMING ALL TIMESTEPS:')
total_clust_matrix=np.zeros((max_n+1,max_n+1))
if check=='y':
	e=0
	f=0
if not big:	
	for i in np.arange(len(clust_matrix)):
		for j in np.arange(max_n+1):
			for k in np.arange(max_n+1):
				total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
				if check=='y':
					e=e+clust_matrix[i][j][k]*(j+k)
					f=f+clust_matrix[i][j][k]
		Progress_Bar(i,len(clust_matrix))
if big:
	for i in np.arange(len(clust_matrix)):
		if np.sum(clust_matrix[i])==1:
			for j in np.arange(len(index_cation)-5,len(index_cation)+1):
				for k in np.arange(len(index_cation)-5,len(index_cation)+1):
					total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
					if check=='y':
						e=e+clust_matrix[i][j][k]*(j+k)
						if clust_matrix[i][j][k]>0:
							f=f+clust_matrix[i][j][k]
		if np.sum(clust_matrix[i])>1:
			for j in np.arange(int(np.sum(clust_matrix[i]))+100):
				for k in np.arange(int(np.sum(clust_matrix[i]))+100):
					total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
					if check=='y':
						e=e+clust_matrix[i][j][k]*(j+k)
						if clust_matrix[i][j][k]>0:
							f=f+clust_matrix[i][j][k]
			for j in np.arange(len(index_cation)-int(np.sum(clust_matrix[i]))-100,len(index_cation)+1):
				for k in np.arange(len(index_cation)-int(np.sum(clust_matrix[i]))-100,len(index_cation)+1):
					total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
					if check=='y':
						e=e+clust_matrix[i][j][k]*(j+k)
						if clust_matrix[i][j][k]>0:
							f=f+clust_matrix[i][j][k]
		Progress_Bar(i,len(clust_matrix))
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %')
if check=='y':
	print('Third Checking:')
	print('Total # Ions = %i' %(e))
	print('Total # Clusters = %i' %(f))
	print('\n')

# Plot the results.
Plot_Cluster_Map(total_clust_matrix)
Plot_Cluster_Probability(total_clust_matrix)
