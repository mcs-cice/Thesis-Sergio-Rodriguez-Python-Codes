'''
    Author: Sergio Rodriguez Peña
'''



######################################################################
#                                                                    #
#   Code that analyzes the saturation level of of the different      # 
#   molecules of the system with respect to a observable molecule.   #
#                                                                    #
#   This code takes as input a process trajectory by the code        #
#   Process_Trajectory_PDB_to_TXT.py.                                #
#                                                                    #
######################################################################



''' IMPORTED PACKAGES '''

import numpy as np
 


''' DEFINED FUNCTIONS '''

# Function that reads the first snapshot of the simulation and analyzes its composition.
# Then asks the user fot the observable molecules.
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
    observable=input('Choose Observable Molecule ID: ')
    obs=input('Choose Points of Coordination (1 atom ID per molecule): ')
    print('\n')
    poc=[]
    p=[]
    for i in range(len(obs)):
        if obs[i]==',':
            poc.append(''.join(p))
            p=[]
            continue
        p.append(obs[i])
    poc.append(''.join(p))
    return(observable,poc)

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
    return(time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index)


# Function that creates the hypercube. It gets the position of the reference atom/molecule without applying PBC
# and the position of the observable atoms/molecules applying PBC.
def Create_Hyper_Cube(t_index):
    hyper_cube=[]
    ref_position=[]
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
                if l[3]==observable and len(l)>5:
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
                if len(l)>5:
                    if l[3] in poc or l[1] in poc:
                        x=float(l[4])
                        y=float(l[5])
                        z=float(l[6])
                        ref_position.append([x,y,z,l[0],l[1],l[2],l[3]])
    return(hyper_cube,ref_position)


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


# Function that checks if the observable atoms/molecules are inside the cutoff radius of the reference atom/molecule.
def Check_Observable(i,ref_in,coordination,ref_pos,hyper_cube):
    for j in range(len(hyper_cube)):
        pos_obs=hyper_cube[j][:3]
        if abs(pos_obs[0]-ref_pos[0])>cutoff:
            continue
        if abs(pos_obs[1]-ref_pos[1])>cutoff:
            continue
        if abs(pos_obs[2]-ref_pos[2])>cutoff:
            continue
        dis=np.sqrt((pos_obs[0]-ref_pos[0])**2+(pos_obs[1]-ref_pos[1])**2+(pos_obs[2]-ref_pos[2])**2)
        if dis<=cutoff:
            coordination[ref_in][i]+=1
            break
    return(coordination)


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

# Ask for the observable molecule.
observable,poc=Type_of_Interaction()

# Get the information of the system.
time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index=System_Analysis()

# Ask for the initial and final time of analysis.
print('Simulation of %i ns composed by %i snapshots every %i ps\n' %(time[-1]/1000,len(time),time[1]-time[0]))
start_time=float(input('Initial Time (ns): '))*1000
end_time=float(input('Final Time (ns): '))*1000 
print('\n')

# Ask for the cutoff radius.
cutoff=float(input('Cutoff Distance (A): '))

# Analysis of the trajectory.
# For each time step, create the hypercube and check the distances between reference and observables.
print('\n')
print('COMPUTING SATURATION LEVEL:')
coordination=[]
for i in range(len(poc)):
    coordination.append([])
for i in range(len(time)):
    if time[i]<start_time:
        continue
    elif time[i]>end_time:
        break
    else:
        for j in range(len(coordination)):
            coordination[j].append(0)
        hyper_cube,ref_position=Create_Hyper_Cube(i)
        mol_alr_coord=[]
        for j in range(len(ref_position)):
            ref_pos=ref_position[j][:3]
            ref_in=molecule_id.index(ref_position[j][6])
            if ref_position[j][5] in mol_alr_coord:
                continue
            coordination=Check_Observable(i,ref_in,coordination,ref_pos,hyper_cube)
            if ref_position[j][6] in poc:
                mol_alr_coord.append(ref_position[j][5])
        Progress_Bar(i,len(range(time.index(start_time),time.index(end_time)+1))) 
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %' , end='\r')

# Saturation of each molecule at each snapshot.
saturation=[]
for i in range(len(coordination)):
    saturation.append([])
    if molecule_id[i]==observable:
        continue
    for j in range(len(coordination[i])):
        if poc[i] not in molecule_id:
            saturation[i].append(coordination[i][j]/len(atom_index[i][atom_id[i].index(poc[i])]))
        else:
            saturation[i].append(coordination[i][j]/len(molecule_index[i]))

# Average saturation.
mean_saturation=[]
for i in range(len(saturation)):
    if molecule_id[i]!=observable:
        mean_saturation.append(np.mean(saturation[i]))

# Print the results.
print('\n')
for i in range(len(mean_saturation)):
    if molecule_id[i]!=observable:
        print('Saturation of %s by %s = %.3f %%' %(molecule_id[i],observable,mean_saturation[i]*100))
