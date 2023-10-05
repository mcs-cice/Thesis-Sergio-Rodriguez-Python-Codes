'''
    Author: Sergio Rodriguez Peña
'''



######################################################
#                                                    #
#   Compilation of functions to analyze the output   #
#   of TRAVIS for RDF and MSD analysis.              #
#                                                    #
######################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt



''' DEFINED FUNCTIONS '''

# Function that reads the MSD data. By default it does not read the time column.
# Time in ps to ns.
# MSD in nm^2 to A^2.
def lect_msd(input_file,time=False):
	file=open(input_file+'.csv')
	data=[]
	for line in file:
		l=line.split()
		data.append(l)
	if time==True:
		t=[]
		for i in np.arange(1,len(data)):
			t.append(float(data[i][0][:-1])/1000)
	msd=[]
	for i in range(1,len(data)):
		msd.append(float(data[i][1][:-1])/10000)
	if time==True:
		return(t,msd)
	else:
		return(msd)


# Functions that reads the RDF data. By default it does not read the distance column.
# Radius in pm to A.
def lect_rdf(input_file,distance=False):
	file=open(input_file+'.csv')
	data=[]
	for line in file:
		l=line.split()
		data.append(l)
	if distance==True:
		r=[]
		for i in np.arange(1,len(data)):
			r.append(float(data[i][0][:-1])/100)
	rdf=[]
	cn=[]
	for i in range(1,len(data)):
		rdf.append(float(data[i][1][:-1]))
		cn.append(float(data[i][2]))
	if distance==True:
		return(r,rdf,cn)
	else:
		return(rdf,cn)


# Function that computes the Diffusion Coefficient from the Fickian regime.
# To identify this region, it creates the straight line tangent with slope 1 to each point 
# in log-log scale and finds the section in which the data deviate less from the
# straight line.
# D in cm^2/s.
def compute_D(msd):
	delta_time=time[1]-time[0]
	time_corr=time[time.index(1):int(len(time))]
	msd_corr=msd[t.index(1):int(len(msd))]
	lower_limit=[]
	upper_limit=[]
	for i in range(len(msd_corr)):
	    lower_limit.append(msd_corr[i]-error_range)
	    upper_limit.append(msd_corr[i]+error_range)
	region_inside_limits=[]    
	for i in range(len(time_corr)):
	    ref_y_lineal=[]
	    start=t_corr[i]
	    for j in time_corr:
	        ref_y_lineal.append(j/(start/msd_corr[i]))
	    region_inside_limits.append([])
	    for j in range(len(t_corr)):
	        if ref_y_lineal[j]>=lower_limit[j] and ref_y_lineal[j]<=upper_limit[j]:
	            region_inside_limits[i].append(t_corr[j]) 
	region_corr=[]  
	for i in range(len(region_inside_limits)):
	    center=t_corr[i]
	    region_corr.append([])
	    for j in range(len(region_inside_limits[i])):
	        if region_inside_limits[i][j]==center:
	            start_j=j
	            break
	    for j in range(start_j):
	        if region_inside_limits[i][j]==(region_inside_limits[i][start_j]-(start_j-j)*delta_t):
	            region_corr[i].append(region_inside_limits[i][j])
	    for j in range(len(region_inside_limits[i])-start_j):
	        if region_inside_limits[i][start_j+j]==(region_inside_limits[i][start_j]+j*delta_t):
	            region_corr[i].append(region_inside_limits[i][j+start_j])
	l=[]
	for i in range(len(region_corr)):
	    l.append(len(region_corr[i]))
	pos_max=l.index(max(l))
	ref_y_lineal=[]
	for i in range(len(t_corr)):
	    ref_y_lineal.append(t_corr[i]/(t_corr[pos_max]/msd_corr[pos_max]))
	x_fit=region_corr[pos_max]
	y_fit=[]
	for i in x_fit:
		y_fit.append(msd_corr[t_corr.index(i)])   
	coef,V=np.polyfit(x_fit,y_fit,1,cov=True)
	print('Computed D')
	return(coef[0]*1e-7/6,np.sqrt(V[0][0])*1e-7/6)


# Function that computes the mean value and error of a set of data.
def make_mean(data):
	mean=[]
	error=[]
	for i in range(len(data[0])):
		a=[]
		for j in range(len(data)):
			a.append(data[j][i])
		mean.append(np.mean(a))
		b=0
		for j in range(len(data)):
			b=b+(data[j][i]-mean[i])**2
		error.append(np.sqrt(b/len(data)))
	return(mean,error)


# Function that finds the Coordination Number from the first minimum of the RDF.
# We start to look from 2 A.
def compute_CN(rdf,cn):
	b=0
	for i in range(len(rdf)):
		if rdf[i]>2:
			a=rdf[i-20:i+20]
			if rdf[i]==min(a):
				b=cn[i]
				break
		if rdf[i]>8:
			for j in range(len(r)):
				if r[j]>3:
					b=cn[j]
					break
			break
	return(b)
