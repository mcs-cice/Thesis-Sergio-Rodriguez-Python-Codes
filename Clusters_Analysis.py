'''
    Author: Sergio Rodriguez Peña
'''



##############################################################################################
#                                                                                            #
#   Code that creates a file with the ionic clusters of the simulation. To do so, starting   #
#   from each Lithium that has not previously been placed in a cluster, in looks for the     #
#   near anions, and then it takes these anaions as starting points for new Lithiums.        #
#                                                                                            #
#   This code takes as input a process trajectory by the code                                #
#   Process_Trajectory_PDB_to_TXT.py.                                                        #
#                                                                                            #
#   This code writes the results in a file Clusters.txt                                      #
#                                                                                            #
##############################################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



''' DEFINED FUNCTIONS '''

# Function that reads the first snapshot of the simulation and analyzes its composition.
# Then asks the user fot the reference atom/molecule (only 1) and the observable atoms/molecules (up to 5 different).
def Type_of_Interaction():
    molecule_id=[]
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='t' and float(l[2])!=0:
                break
            if len(l)==7:
                if l[3] not in molecule_id:
                    molecule_id.append(l[3])
    print('\n')
    print('Molecule ID\n')
    for i in range(len(molecule_id)):
        print('%s'%(molecule_id[i]))
    print('\n')
    cation=input('Choose Molecule ID of the cation: ')
    anion=input('Choose Molecule ID of the anion: ')
    return(cation,anion)


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
    n_atoms=0
    for i in range(len(atom_index)):
        for j in range(len(atom_index[i])):
            n_atoms+=len(atom_index[i][j])
    return(time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms)


# Function that creates the hypercube. It gets the position of the reference atom/molecule without applying PBC
# and the position of the observable atoms/molecules applying PBC.
def Create_Hyper_Cube(t_index):
    anion_hyper_cube=[]
    cation_position=[]
    with open(input_file,'r') as r:
        r.seek(snap_offset[t_index])
        correct_snap=False
        for line in r:
            l=line.split()
            if l[0]=='t' and float(l[2])==time[t_index]:
                correct_snap=True
            if l[0]=='t' and float(l[2])!=time[t_index]:
                correct_snap=False
            if not correct_snap:
                break
            if correct_snap:
                if l[3]==anion:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    anion_hyper_cube.append([x,y,z,l[2]])
                    pbc=Check_PBC([x,y,z])
                    if pbc!=[0,0,0]:
                        if pbc[0]!=0:
                            anion_hyper_cube.append([x+a_box*pbc[0],y,z,l[2]])
                        if pbc[1]!=0:
                            anion_hyper_cube.append([x,y+a_box*pbc[1],z,l[2]])
                        if pbc[2]!=0:
                            anion_hyper_cube.append([x,y,z+a_box*pbc[2],l[2]])
                        if pbc[0]!=0 and pbc[1]!=0:
                            anion_hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z,l[2]])
                        if pbc[0]!=0 and pbc[2]!=0:
                            anion_hyper_cube.append([x+a_box*pbc[0],y,z+a_box*pbc[2],l[2]])
                        if pbc[1]!=0 and pbc[2]!=0:
                            anion_hyper_cube.append([x,y+a_box*pbc[1],z+a_box*pbc[2],l[2]])
                        if pbc[0]!=0 and pbc[1]!=0 and pbc[2]!=0:
                            anion_hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z+a_box*pbc[2],l[2]])
                if l[3]==cation:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    cation_position.append([x,y,z,l[2]])
    return(anion_hyper_cube,cation_position)


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


# Function that checks for new elements in the cluster.
def Check_New_Elements(Clusters,num_cl,s):
    if int(Clusters[num_cl][s]) in molecule_index[molecule_id.index(cation)]:
        for j in range(len(cation_position)):
            if cation_position[j][3]==Clusters[num_cl][s]:
                pos_ref=cation_position[j][:3]
        for k in range(len(anion_hyper_cube)):
            if anion_hyper_cube[k][3] in Clusters[num_cl]:
                continue
            elif abs(pos_ref[0]-anion_hyper_cube[k][0])>cutoff or abs(pos_ref[1]-anion_hyper_cube[k][1])>cutoff or abs(pos_ref[2]-anion_hyper_cube[k][2])>cutoff:
                continue
            else:
                d=np.sqrt((pos_ref[0]-anion_hyper_cube[k][0])**2+(pos_ref[1]-anion_hyper_cube[k][1])**2+(pos_ref[2]-anion_hyper_cube[k][2])**2)
                if d<cutoff:
                    Clusters[num_cl].append(anion_hyper_cube[k][3])
        j=0
        while j<len(cation_position):
            if cation_position[j][3]==Clusters[num_cl][s]:
                del cation_position[j]
            else:
                j+=1
    if int(Clusters[num_cl][s]) in molecule_index[molecule_id.index(anion)]:
        for k in range(len(anion_hyper_cube)):
            if anion_hyper_cube[k][3]==Clusters[num_cl][s]:
                pos_ref=anion_hyper_cube[k][:3]
                for m in range(len(cation_position)):
                    if cation_position[m][3] in Clusters[num_cl]:
                        continue
                    elif abs(pos_ref[0]-cation_position[m][0])>cutoff or abs(pos_ref[1]-cation_position[m][1])>cutoff or abs(pos_ref[2]-cation_position[m][2])>cutoff:
                        continue
                    else:
                        d=np.sqrt((pos_ref[0]-cation_position[m][0])**2+(pos_ref[1]-cation_position[m][1])**2+(pos_ref[2]-cation_position[m][2])**2)
                        if d<cutoff:
                            Clusters[num_cl].append(cation_position[m][3])  
        j=0
        while j<len(anion_hyper_cube):
            if anion_hyper_cube[j][3]==Clusters[num_cl][s]:
                del anion_hyper_cube[j]
            else:
                j+=1
    return(Clusters)


# Function that creates the progress bar shown in the terminal.
def Progress_Bar(progress,total):
    percent=100*(progress/float(total))
    bar='|'+'#'*int(percent/2)+'-'*(50-(int(percent/2)))+'|'
    print(f'%s %.1f' %(bar,percent) +' %' , end='\r')



''' CODE '''

# Name of the input file.
print('\n')
entr=input('Trayectory File Name: ')
input_file=entr+'.txt'

# Ask for the reference and observables.
cation,anion=Type_of_Interaction()

# Get the information of the system.
time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms=System_Analysis()

# Ask for the initial and final time of analysis.
print('\n')
print('Simulation of %i ns composed by %i snapshots every %i ps\n' %(time[-1]/1000,len(time),time[1]-time[0]))
start_time=float(input('Initial Time (ns): '))*1000
end_time=float(input('Final Time (ns): '))*1000 
print('\n')

# Ask for the cutoff radius.
cutoff=float(input('Cutoff Distance (A): '))

# Analysis of the trajectory.
# For each time step, create the hypercube and check the distances between the molecules that form the clusters.
# Write the Cluster.txt file.
print('\n')
print('COMPUTING COORDINATION ANALYSIS:')
file=open('Clusters.txt','w')
for i in range(len(time)):
    if time[i]<start_time:
        continue
    elif time[i]>end_time:
        break
    else:
        file.write('Anion: %i - %i \n' %(molecule_index[molecule_id.index(anion)][0],molecule_index[molecule_id.index(anion)][-1]))
        file.write('Cation: %i - %i \n' %(molecule_index[molecule_id.index(cation)][0],molecule_index[molecule_id.index(cation)][-1]))
        Clusters=[]
        num_cl=-1
        anion_hyper_cube,cation_position=Create_Hyper_Cube(i)
        while len(cation_position)>0:
            Clusters.append([cation_position[0][3]])
            num_cl+=1
            s=0
            while s<len(Clusters[num_cl]):
                Clusters=Check_New_Elements(Clusters,num_cl,s)
                s+=1
    file.write('t = %i ps\n' %(time[i]))
    for j in range(len(Clusters)):
        for k in range(len(Clusters[j])):
            file.write('%s ' %(Clusters[j][k]))
        file.write('\n')
        file.write('\n')
    Progress_Bar(i-time.index(start_time),len(range(time.index(start_time),time.index(end_time)+1)))
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %' , end='\r')
print('\n')

file.close()

