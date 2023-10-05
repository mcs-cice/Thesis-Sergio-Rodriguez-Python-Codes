'''
    Author: Sergio Rodriguez Peña
'''



#####################################################################################
#                                                                                   #
#   Code that analyzes the movement of the reference atoms through the simulation   #
#   by showing the points of coordination with the observable molecules/atoms.      #
#                                                                                   #
#   This code takes as input a process trajectory by the code                       #
#   Process_Trajectory_PDB_to_TXT.py.                                               #   
#                                                                                   #
#####################################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt



''' DEFINED FUNCTIONS '''

# Function that reads the first snapshot of the simulation and analyzes its composition.
# Then asks the user fot the reference atom/molecule (only 1) and the observable atoms/molecules (up to 5 different).
def Type_of_Interaction():
    molecule_id=[]
    atom_id=[]
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='t' and float(l[2])!=0:
                break
            if len(l)==7:
                if l[3] not in molecule_id:
                    molecule_id.append(l[3])
                    atom_id.append([])
                if l[1] not in atom_id[molecule_id.index(l[3])]:
                    atom_id[molecule_id.index(l[3])].append(l[1])
    print('\n')
    print('Molecule ID: Atom ID\n')
    for i in range(len(molecule_id)):
        print('%s: '%(molecule_id[i]), *atom_id[i])
    print('\n')
    reference=input('Choose Reference Molecule ID: ')
    obs=input('Choose Observable Atom/Molecule ID (up to 5 elements separated by comas): ')
    print('\n')
    observable=[]
    o=[]
    for i in range(len(obs)):
        if obs[i]==',':
            observable.append(''.join(o))
            o=[]
            continue
        o.append(obs[i])
    observable.append(''.join(o))
    return(reference,observable,molecule_id,atom_id)


# Function that reads the trajectory and returns all the esential information of the system:
# time of the snapshots, size of the box, number of atoms, indexes of the molecules and atoms.
def System_Analysis():
    time=[]
    snap_offset=[]
    offset=0
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='t':
                time.append(float(l[2]))  
                snap_offset.append(offset)
            offset+=(len(line)+1)
    molecule_id=[]
    molecule_index=[]
    atom_id=[]
    atom_index=[]
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='Box':
                a_box=float(l[3])
            if len(l)==7:
                if l[3] not in molecule_id:
                    molecule_id.append(l[3])
                    molecule_index.append([])
                    atom_id.append([])
                    atom_index.append([])
                if l[1] not in atom_id[molecule_id.index(l[3])]:
                    atom_id[molecule_id.index(l[3])].append(l[1])
                    atom_index[molecule_id.index(l[3])].append([])
                if int(l[2]) not in molecule_index[molecule_id.index(l[3])]:
                    molecule_index[molecule_id.index(l[3])].append(int(l[2]))
                atom_index[molecule_id.index(l[3])][atom_id[molecule_id.index(l[3])].index(l[1])].append(int(l[0]))
            if l[0]=='t' and float(l[2])!=0:
                break 
    obs_poly=np.zeros((len(poly)))
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if len(l)==7:
                if l[3] in poly and l[1] in observable:
                 for i in range(len(poly)):
                     if int(l[2])==molecule_index[molecule_id.index(poly[i])][0]:
                         obs_poly[i]+=1
            if l[0]=='t' and float(l[2])!=0:
                break 
    n_atoms=0
    for i in range(len(atom_index)):
        for j in range(len(atom_index[i])):
            n_atoms+=len(atom_index[i][j])
    return(time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms,obs_poly)


# Function that creates the hypercube. It gets the position of the reference atom/molecule without applying PBC
# and the position of the observable atoms/molecules applying PBC.
def Create_Hyper_Cube(t_index):
    hyper_cube=[]
    pos_ref=[]
    with open(input_file,'r') as r:
        r.seek(snap_offset[t_index])
        correct_snap=False
        for line in r:
            l=line.split()
            if len(l)<3:
                    continue
            if l[0]=='t' and float(l[2])==time[t_index]:
                correct_snap=True
            if l[0]=='t' and float(l[2])!=time[t_index]:
                correct_snap=False
            if not correct_snap:
                break
            if correct_snap:
                if len(l)<7:
                    continue
                if l[3] in observable or l[1] in observable:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    hyper_cube.append([x,y,z,l[0],l[1],l[2],l[3]])
                    pbc=Check_PBC([x,y,z])
                    if pbc!=[0,0,0]:
                        if pbc[0]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y,z,l[0],l[1],l[2],l[3]])
                        if pbc[1]!=0:
                            hyper_cube.append([x,y+a_box*pbc[1],z,l[0],l[1],l[2],l[3]])
                        if pbc[2]!=0:
                            hyper_cube.append([x,y,z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[1]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z,l[0],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[2]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y,z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                        if pbc[1]!=0 and pbc[2]!=0:
                            hyper_cube.append([x,y+a_box*pbc[1],z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[1]!=0 and pbc[2]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                if l[1]==reference or l[3]==reference:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    pos_ref.append([x,y,z,l[0],l[1],l[2],l[3]])
    return(hyper_cube,pos_ref)


# Function that checks if a position is close enought to the edges of the box and needs to be applied PBC.
def Check_PBC(pos):
    pbc=[]
    for i in pos:
        if i>(a_box-cutoff+0.01):
            pbc.append(-1)
        elif i<(cutoff+0.01):
            pbc.append(1)
        else:
            pbc.append(0)
    return(pbc)


# Function that checks if the observable atom/molecules are inside the cutoff radius of the reference atom.
def Check_Observable(i,j,pos_ref,coordination,coordination_info,hyper_cube):
    for k in range(len(hyper_cube)):
        pos_obs=hyper_cube[k][:3]
        if abs(pos_obs[0]-pos_ref[0])>cutoff:
            continue
        if abs(pos_obs[1]-pos_ref[1])>cutoff:
            continue
        if abs(pos_obs[2]-pos_ref[2])>cutoff:
            continue
        dis=np.sqrt((pos_ref[0]-pos_obs[0])**2+(pos_ref[1]-pos_obs[1])**2+(pos_ref[2]-pos_obs[2])**2)
        if dis<=cutoff:
            coordination[j][i].append(int(hyper_cube[k][3]))
            coordination_info[j][i].append([hyper_cube[k][4],hyper_cube[k][5],hyper_cube[k][6]])
    return(coordination,coordination_info)


# Function that creates the progress bar shown in the terminal.
def Progress_Bar(progress,total):
    percent=100*(progress/float(total))
    bar='|'+'#'*int(percent/2)+'-'*(50-(int(percent/2)))+'|'
    print(f'%s %.1f' %(bar,percent) +' %' , end='\r')


# Function that creates a correlation between the index of the atoms in the input file
# and the index it will be used to plot the results.   
def Correlate_Index():
    correlation_coordination=[[],[]]
    actual_index=np.zeros((len(observable)))
    ref_mol=np.zeros((len(observable)))
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='t' and float(l[2])!=0:
                break
            if l[0]=='t' or l[0]=='Box':
                continue
            if l[1] in observable:
                correlation_coordination[0].append(int(l[0]))
                if l[3] not in poly:
                    if int(l[2])!=ref_mol[observable.index(l[1])]:
                        ref_mol[observable.index(l[1])]=int(l[2])
                        actual_index[observable.index(l[1])]+=1
                    correlation_coordination[1].append(int(actual_index[observable.index(l[1])]))
                if l[3] in poly:
                    actual_index[observable.index(l[1])]+=1
                    correlation_coordination[1].append(int(actual_index[observable.index(l[1])]))
            if l[3] in observable:
                correlation_coordination[0].append(int(l[0]))
                if int(l[2])!=ref_mol[observable.index(l[3])]:
                    ref_mol[observable.index(l[3])]=int(l[2])
                    actual_index[observable.index(l[3])]+=1
                correlation_coordination[1].append(int(actual_index[observable.index(l[3])]))
    return(correlation_coordination)


# Function that changes the index of the coordination points to the correlated ones,
# and creates the variables to plot the results.
def Change_Index():
    x=[]
    y=[]
    for i in range(len(in_ref)):
        x.append([])
        y.append([])
        nums_used=[]
        for j in range(len(observable)):
            x[i].append([])
            y[i].append([])
            nums_used.append([])
        for j in range(len(coordination[i])):
            for k in range(len(coordination[i][j])):
                if coordination_info[i][j][k][2] in observable:
                    if int(coordination_info[i][j][k][1]) not in nums_used[observable.index(coordination_info[i][j][k][2])]:
                        nums_used[observable.index(coordination_info[i][j][k][2])].append(int(coordination_info[i][j][k][1]))
                    if coordination_info[i][j][k][2] not in poly:
                        x[i][observable.index(coordination_info[i][j][k][2])].append(time[j]/1000)
                        new=correlation_coordination[1][correlation_coordination[0].index(coordination[i][j][k])]
                        y[i][observable.index(coordination_info[i][j][k][2])].append(nums_used[observable.index(coordination_info[i][j][k][2])].index(int(coordination_info[i][j][k][1]))+1)
                    if coordination_info[i][j][k][2] in poly:
                        x[i][observable.index(coordination_info[i][j][k][2])].append(time[j]/1000)
                        new=correlation_coordination[1][correlation_coordination[0].index(coordination[i][j][k])]
                        if new%obs_poly[poly.index(coordination_info[i][j][k][2])]==0:
                            y[i][observable.index(coordination_info[i][j][k][2])].append(obs_poly[poly.index(coordination_info[i][j][k][2])]+obs_poly[poly.index(coordination_info[i][j][k][2])]*nums_used[observable.index(coordination_info[i][j][k][2])].index(int(coordination_info[i][j][k][1])))
                        else:
                            y[i][observable.index(coordination_info[i][j][k][2])].append(new%obs_poly[poly.index(coordination_info[i][j][k][2])]+obs_poly[poly.index(coordination_info[i][j][k][2])]*nums_used[observable.index(coordination_info[i][j][k][2])].index(int(coordination_info[i][j][k][1])))           
                if coordination_info[i][j][k][0] in observable:
                    if int(coordination_info[i][j][k][1]) not in nums_used[observable.index(coordination_info[i][j][k][0])]:
                        nums_used[observable.index(coordination_info[i][j][k][0])].append(int(coordination_info[i][j][k][1]))
                    if coordination_info[i][j][k][2] not in poly:
                        x[i][observable.index(coordination_info[i][j][k][0])].append(time[j]/1000)
                        new=correlation_coordination[1][correlation_coordination[0].index(coordination[i][j][k])]
                        y[i][observable.index(coordination_info[i][j][k][0])].append(nums_used[observable.index(coordination_info[i][j][k][0])].index(int(coordination_info[i][j][k][1]))+1)
                    if coordination_info[i][j][k][2] in poly:
                        x[i][observable.index(coordination_info[i][j][k][0])].append(time[j]/1000)
                        new=correlation_coordination[1][correlation_coordination[0].index(coordination[i][j][k])]
                        if new%obs_poly[poly.index(coordination_info[i][j][k][2])]==0:
                            y[i][observable.index(coordination_info[i][j][k][0])].append(obs_poly[poly.index(coordination_info[i][j][k][2])]+obs_poly[poly.index(coordination_info[i][j][k][2])]*nums_used[observable.index(coordination_info[i][j][k][0])].index(int(coordination_info[i][j][k][1])))
                        else:
                            y[i][observable.index(coordination_info[i][j][k][0])].append(new%obs_poly[poly.index(coordination_info[i][j][k][2])]+obs_poly[poly.index(coordination_info[i][j][k][2])]*nums_used[observable.index(coordination_info[i][j][k][0])].index(int(coordination_info[i][j][k][1])))           
    return(x,y,nums_used)


# Function that plots the results.
def Plot_Results(x,y):
    num_subplots=0
    for i in range(len(y)):
        if len(y[i])>0:
            num_subplots+=1
    if num_subplots>1:
        fig,axs=plt.subplots(num_subplots,1,sharex=True)
        color=['k','r','b','g','m']
        n=0
        for i in range(len(x)):
            if len(y[i])==0:
                continue
            for j in range(len(atom_id)):
                if observable[i] in atom_id[j]:
                    ch=molecule_id[j]
            if observable[i] in molecule_id:
                ch=molecule_id[molecule_id.index(observable[i])]
            if ch in poly:
                if max(y[i])*1.2<obs_poly[poly.index(molecule_id[i])]:
                    del_mol=np.arange(0.5,obs_poly[poly.index(molecule_id[i])]+5,obs_poly[poly.index(molecule_id[i])])
                else:
                    del_mol=np.arange(0.5,max(y[i])*1.2,obs_poly[poly.index(molecule_id[i])])
                axs[n].scatter(x[i],y[i],s=2,color=color[i])
            else:
                del_mol=np.arange(0.5,max(y[i])+1)
                axs[n].scatter(x[i],y[i],s=5,color=color[i])
            if i==len(observable)-1:
                axs[n].set_xlabel('t / ns',fontsize=30)
            if ch in poly:
                axs[n].set_ylabel('%s index' %(observable[i]),fontsize=30)
            else:
                axs[n].set_ylabel('%s index' %(ch),fontsize=30)
            axs[n].set_xlim(0,time[-1]/1000)
            axs[n].tick_params(axis='both',labelsize=30)
            for j in range(len(del_mol)):
                axs[n].hlines(del_mol[j],time[0],time[-1]/1000,color='black',linestyle='dashed')
            axs[n].set_ylim(0,max(del_mol)+1)
            n+=1
    if num_subplots==1:
        plt.figure()
        color=['k','r','b','g','m']
        n=0
        for i in range(len(x)):
            if len(y[i])==0:
                continue
            for j in range(len(atom_id)):
                if observable[i] in atom_id[j]:
                    ch=molecule_id[j]
            if observable[i] in molecule_id:
                ch=molecule_id[molecule_id.index(observable[i])]
            if ch in poly:
                del_mol=range(0.5,max(y[i])*1.2,obs_poly[poly.index(molecule_id[i])])
            else:
                del_mol=range(0.5,max(y[i])+1)
            plt.scatter(x[i],y[i],s=2,color=color[i])
            plt.xlabel('t / ns',fontsize=30)
            if ch in poly:
                plt.ylabel('%s index' %(observable[i]),fontsize=30)
            else:
                plt.ylabel('%s index' %(ch),fontsize=30)
            plt.xlim(0,time[-1]/1000)
            plt.tick_params(axis='both',labelsize=30)
            for j in range(len(del_mol)):
                plt.hlines(del_mol[j],time[0],time[-1]/1000,color='black',linestyle='dashed')
            plt.ylim(0,max(del_mol)+1)



''' CODE '''

# Name of the input file.
print('\n')
entr=input('Trayectory File Name: ')
input_file=entr+'.txt'

# Ask for the reference and observables.
reference,observable,molecule_id,atom_id=Type_of_Interaction()

# Ask if any of the molecules requires a polymer analysis.
print('\n')
is_poly=input('Is any of the molecules a Polymer? (y/n): ')
print('\n')
poly=['']
if is_poly=='y':
    for i in range(len(molecule_id)):
        print('%s = %i' %(molecule_id[i],i))
    print('\n')
    which_poly=input('Which one?: ')
    pol=[]
    p=[]
    for i in range(len(which_poly)):
        if which_poly[i]==',':
            pol.append(''.join(p))
            p=[]
            continue
        p.append(which_poly[i])
    pol.append(''.join(p))
    if is_poly=='y':
        poly=[]
        for i in range(len(pol)):
            poly.append(molecule_id[int(pol[i])])

# Get the information of the system.
time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms,obs_poly=System_Analysis()

# Randomly choose which reference atoms analyze. 
print('\n')
num_ref=int(input('How many %s analyze?: ' %(reference)))
for i in range(1000):
    if reference in molecule_id:
        in_ref=list(np.random.randint(1,len(molecule_index[molecule_id.index(reference)]),num_ref))
    else:
        in_ref=list(np.random.randint(1,len(atom_index[molecule_id.index(reference)][atom_id[molecule_id.index(reference)].index(reference)]),num_ref))
    if len(set(in_ref))==num_ref:
        break

# Ask for the cutoff radius.
print('\n')
cutoff=float(input('Cutoff Distance (A): '))

# Analysis of the trajectory.
# For each time step, create the hypercube and check the distances between reference and observables.
print('\n')
print('ANALYZING COORDINATION ENVIRONMENT:')
coordination=[]
coordination_info=[]
for i in range(len(in_ref)):
    coordination.append([])
    coordination_info.append([])
for i in range(len(time)):
    for j in range(len(in_ref)):
        coordination[j].append([])
        coordination_info[j].append([])
    hyper_cube,position_ref=Create_Hyper_Cube(i)
    for j in range(len(in_ref)):
        pos_ref=position_ref[in_ref[j]]
        coordination,coordination_info=Check_Observable(i,j,pos_ref,coordination,coordination_info,hyper_cube)
    Progress_Bar(i,len(time))  
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %' , end='\r')
print('\n')

# Correlate the index to order them.
correlation_coordination=Correlate_Index()

# Change the indexes.
x,y,nums_used=Change_Index()

# Plot the results
for i in range(len(x)):
    Plot_Results(x[i],y[i])
plt.show(block=False)
