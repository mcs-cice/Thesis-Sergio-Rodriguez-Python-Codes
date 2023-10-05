'''
    Author: Sergio Rodriguez Peña
'''



###############################################################################################
#                                                                                             #
#   Code that analyzes the composition of the first coordination shell of a reference         #
#   molecule of the simulation by searching the observable atoms/molecules                    #
#   inside a cutoff radius. Then it plots the results in a Venn diagram if possible.          #
#   The code takes a maximum of 1 reference atom/molecule and 5 observable atoms/molecules.   #
#   The Venn diagram is plotted if 2 or 3 observables are analyzed.                           #
#                                                                                             #
#   This code takes as input a process trajectory by the code                                 #
#   Process_Trajectory_PDB_to_TXT.py.                                                         #
#                                                                                             #
#   This code writes a file Coord.txt with the exact results.                                 #
#                                                                                             #
###############################################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn as vplt
import pandas as pd
import itertools



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
    return(reference,observable)


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


# Function that creates the ''hypercube''. It gets the position of the reference atom/molecule 
# without applying PBC and the position of the observable atoms/molecules applying PBC.
def Create_Hyper_Cube(t_index):
    hyper_cube=[]
    pos_ref=[]
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
                if l[1] in observable or l[3] in observable:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    hyper_cube.append([x,y,z,l[1],l[2],l[3]])
                    pbc=Check_PBC([x,y,z])
                    if pbc!=[0,0,0]:
                        if pbc[0]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y,z,l[1],l[2],l[3]])
                        if pbc[1]!=0:
                            hyper_cube.append([x,y+a_box*pbc[1],z,l[1],l[2],l[3]])
                        if pbc[2]!=0:
                            hyper_cube.append([x,y,z+a_box*pbc[2],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[1]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z,l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[2]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y,z+a_box*pbc[2],l[1],l[2],l[3]])
                        if pbc[1]!=0 and pbc[2]!=0:
                            hyper_cube.append([x,y+a_box*pbc[1],z+a_box*pbc[2],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[1]!=0 and pbc[2]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z+a_box*pbc[2],l[1],l[2],l[3]])
                if l[1]==reference or l[3]==reference:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    pos_ref.append([x,y,z,l[2]])
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


# Function that checks if the observable atoms/molecules are inside the cutoff radius of the reference atom/molecule.
def Check_Observable(pos_ref,actual_ref,hyper_cube):
    m_n=[]
    for i in range(len(observable)):
        m_n.append([])
    for i in range(len(hyper_cube)):
        if hyper_cube[i][4]==actual_ref:
            continue
        pos_obs=hyper_cube[i][0:3]
        if abs(pos_obs[0]-pos_ref[0])>cutoff:
            continue
        if abs(pos_obs[1]-pos_ref[1])>cutoff:
            continue
        if abs(pos_obs[2]-pos_ref[2])>cutoff:
            continue
        dis=np.sqrt((pos_obs[0]-pos_ref[0])**2+(pos_obs[1]-pos_ref[1])**2+(pos_obs[2]-pos_ref[2])**2)
        if dis<=cutoff:
            if hyper_cube[i][3] in observable:
                obs_id=observable.index(hyper_cube[i][3])
            if hyper_cube[i][5] in observable:
                obs_id=observable.index(hyper_cube[i][5])
            m_n[obs_id].append(hyper_cube[i][4])
    return(m_n)


# Function that updates the variable coordination where the number of the different compositions of the 
# first coordination shell are saved.
def coordination_save(coordination,molecule_number):
    lon=[]
    for i in range(len(molecule_number)):
        lon.append(len(set(molecule_number[i])))
    if np.sum(lon)==0:
        coordination[0]+=1
    else:
        coord_id,num_coord=indexes(lon)
        if num_coord==1:
            while len(coordination[1][coord_id[0]])<len(set(molecule_number[coord_id[0]])):
                coordination[1][coord_id[0]].append(0)
            coordination[1][coord_id[0]][len(set(molecule_number[coord_id[0]]))-1]+=1
        if num_coord==2:
            x=int(find_double(coord_id,molecule_number))
            while len(coordination[2][x])<len(set(molecule_number[coord_id[0]])):
                coordination[2][x].append([0])
            while len(coordination[2][x][len(set(molecule_number[coord_id[0]]))-1])<len(set(molecule_number[coord_id[1]])):
                coordination[2][x][len(set(molecule_number[coord_id[0]]))-1].append(0)
            coordination[2][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1]+=1
        if num_coord==3:
            x=int(find_triple(coord_id,molecule_number))
            while len(coordination[3][x])<len(set(molecule_number[coord_id[0]])):
                coordination[3][x].append([[0]])
            while len(coordination[3][x][len(set(molecule_number[coord_id[0]]))-1])<len(set(molecule_number[coord_id[1]])):
                coordination[3][x][len(set(molecule_number[coord_id[0]]))-1].append([0])
            while len(coordination[3][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1])<len(set(molecule_number[coord_id[2]])):
                coordination[3][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1].append(0)
            coordination[3][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1]+=1
        if num_coord==4:
            x=int(find_quadruple(coord_id,molecule_number))
            while len(coordination[4][x])<len(set(molecule_number[coord_id[0]])):
                coordination[4][x].append([[[0]]])
            while len(coordination[4][x][len(set(molecule_number[coord_id[0]]))-1])<len(set(molecule_number[coord_id[1]])):
                coordination[4][x][len(set(molecule_number[coord_id[0]]))-1].append([[0]])
            while len(coordination[4][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1])<len(set(molecule_number[coord_id[2]])):
                coordination[4][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1].append([0])
            while len(coordination[4][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1])<len(set(molecule_number[coord_id[3]])):
                coordination[4][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1].append(0)
            coordination[4][x][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1][len(set(molecule_number[coord_id[3]]))-1]+=1
        if num_coord==5:
            while len(coordination[5])<len(set(molecule_number[coord_id[0]])):
                coordination[5].append([[[[0]]]])
            while len(coordination[5][len(set(molecule_number[coord_id[0]]))-1])<len(set(molecule_number[coord_id[1]])):
                coordination[5][len(set(molecule_number[coord_id[0]]))-1].append([[[0]]])
            while len(coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1])<len(set(molecule_number[coord_id[2]])):
                coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1].append([[0]])
            while len(coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1])<len(set(molecule_number[coord_id[3]])):
                coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1].append([0])
            while len(coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1][len(set(molecule_number[coord_id[3]]))-1])<len(set(molecule_number[coord_id[4]])):
                coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1][len(set(molecule_number[coord_id[3]]))-1].append(0)
            coordination[5][len(set(molecule_number[coord_id[0]]))-1][len(set(molecule_number[coord_id[1]]))-1][len(set(molecule_number[coord_id[2]]))-1][len(set(molecule_number[coord_id[3]]))-1][len(set(molecule_number[coord_id[4]]))-1]+=1
    return(coordination)  


# Function that gets the observable indexes of the specific composition analyzed.
def indexes(lon):
    coord_id=[]
    num_coord=0
    for i in range(len(lon)):
        if lon[i]>0:
            coord_id.append(i)
            num_coord+=1
    return(coord_id,num_coord)


# Function that gets the coord_id of a double coordination.
def find_double(coord_id,molecule_number):
    a=0
    for i in range(len(molecule_number)-coord_id[0],len(molecule_number)):
        a+=i
    x=a+(coord_id[1]-coord_id[0]-1)
    return(x)


# Function that gets the coord_id of a triple coordination.
def find_triple(coord_id,molecule_number):
    a=0
    for i in range(len(molecule_number)-coord_id[0],len(molecule_number)):
        a=a+(np.math.factorial(i)/(np.math.factorial(2)*np.math.factorial(i-2)))
    for i in range(len(molecule_number)-(coord_id[1]-coord_id[0]),len(molecule_number)-1):
        a=a+i
    x=a+(coord_id[2]-coord_id[1]-1)
    return(x)


# Function that gets the coord_id of a quadrupole coordination.
def find_quadruple(coord_id,molecule_number):
    a=0
    for i in range(len(molecule_number)-coord_id[0],len(molecule_number)):
        a=a+(np.math.factorial(i)/(np.math.factorial(3)*np.math.factorial(i-3)))
    for i in range(len(molecule_number)-(coord_id[1]-coord_id[0]),len(molecule_number)-1):
        a=a+(np.math.factorial(i)/(np.math.factorial(2)*np.math.factorial(i-2)))
    x=a+(coord_id[3]-coord_id[2]-1)
    return(x)


# Function that creates the progress bar shown in the terminal.
def Progress_Bar(progress,total):
    percent=100*(progress/float(total))
    bar='|'+'#'*int(percent/2)+'-'*(50-(int(percent/2)))+'|'
    print(f'%s %.1f' %(bar,percent) +' %' , end='\r')


# Function that writes the results of the analysis in the Coord.txt file.
def write_coordination(coordination,combination):
    print('\n')
    f=open('Coord.txt','w')
    f.write('COORDINATION RESULTS: \n')
    f.write('\n')
    f.write('\n')
    # Isolated
    f.write('Isolated %s Atoms: %i \n' %(reference,coordination[0]))
    f.write('\n')
    f.write('\n')
    # Single Coordination
    f.write('Coordination with 1 type of molecule: \n')
    f.write('\n')
    for i in range(len(coordination[1])):
        f.write('%s coordinating only with %s: \n' %(reference,combination[1][i]))
        for j in range(len(coordination[1][i])):
            f.write('%i different molecule: %i \n' %(j+1,coordination[1][i][j]))
        f.write('\n')
    # Double Coordination
    f.write('\n')
    f.write('Coordination with 2 types of molecules: \n')
    for i in range(len(coordination[2])):
        if len(coordination[2][i])==0:
            f.write('\n')
            f.write('%s coordinating with %s and %s: \n' %(reference,combination[2][i][0],combination[2][i][1]))
            f.write('\n')
            continue
        else:
            f.write('\n')
            f.write('%s coordinating with %s and %s: \n' %(reference,combination[2][i][0],combination[2][i][1]))
            a=[]
            b=[]
            x=[]
            for j in range(len(coordination[2][i])):
                x.append(len(coordination[2][i][j]))
            x=max(x)
            for j in range(x):
                a.append('Diff. %s' %(combination[2][i][1]))
                b.append('%i' %(j+1))
            header=[np.array(a),np.array(b)]
            c=[]
            d=[]
            for j in range(len(coordination[2][i])):
                c.append('Diff. %s' %(combination[2][i][0]))
                d.append('%i' %(j+1))
            coord_id=[np.array(c),np.array(d)]
            data=np.zeros((len(coordination[2][i]),x))
            for j in range(len(coordination[2][i])):
                for k in range(len(coordination[2][i][j])):
                    data[j][k]=coordination[2][i][j][k]
            df=pd.DataFrame(data,index=coord_id,columns=header)
            df_string = df.to_string(header=True, index=True)
            f.write(df_string)
            f.write('\n')
    # Triple Coordination
    if len(coordination)>3:
        f.write('\n')
        f.write('\n')
        f.write('Coordination with 3 types of molecules: \n')
        for i in range(len(coordination[3])):
            if len(coordination[3][i])==0:
                f.write('\n')
                f.write('%s coordinating with %s, %s and %s: \n' %(reference,combination[3][i][0],combination[3][i][1],combination[3][i][2]))
                f.write('\n')
                continue
            else:
                for j in range(len(coordination[3][i])):
                    f.write('\n')
                    f.write('%s coordinating with %i %s, %s and %s: \n' %(reference,j+1,combination[3][i][0],combination[3][i][1],combination[3][i][2]))
                    a=[]
                    b=[]
                    x=[]
                    for k in range(len(coordination[3][i][j])):
                        x.append(len(coordination[3][i][j][k]))
                    if len(x)!=0:
                        x=max(x)
                    for k in range(x):
                        a.append('Diff. %s' %(combination[3][i][1]))
                        b.append('%i' %(k+1))
                    header=[np.array(a),np.array(b)]
                    c=[]
                    d=[]
                    for k in range(len(coordination[3][i][j])):
                        c.append('Diff. %s' %(combination[3][i][2]))
                        d.append('%i' %(k+1))
                    coord_id=[np.array(c),np.array(d)]
                    data=np.zeros((len(coordination[3][i][j]),x))
                    for k in range(len(coordination[3][i][j])):
                        for m in range(len(coordination[3][i][j][k])):
                            data[k][m]=coordination[3][i][j][k][m]
                    df=pd.DataFrame(data,index=coord_id,columns=header)
                    df_string = df.to_string(header=True, index=True)
                    f.write(df_string)
                    f.write('\n')
    # Quadrupole Coordination
    if len(coordination)>4:
        f.write('\n')
        f.write('\n')
        f.write('Coordination with 4 types of molecules: \n')
        for i in range(len(coordination[4])):
            if len(coordination[4][i])==0:
                f.write('%s coordinating with %s, %s, %s and %s: \n' %(reference,combination[4][i][0],combination[4][i][1],combination[4][i][2],combination[4][i][3]))
                f.write('\n')
                f.write('\n')
                continue
            else:
                for j in range(len(coordination[4][i])):
                    for k in range(len(coordination[4][i][j])):
                        f.write('\n')
                        f.write('%s coordinating with %i %s, %i %s, %s and %s: \n' %(reference,i+1,combination[4][i][0],j+1,combination[4][i][1],combination[4][i][2],combination[4][i][3]))
                        a=[]
                        b=[]
                        x=[]
                        for m in range(len(coordination[4][i][j][k])):
                            x.append(len(coordination[4][i][j][k][m]))
                        if len(x)!=0:
                            x=max(x)
                        for m in range(x):
                            a.append('Diff. %s' %(combination[4][i][2]))
                            b.append('%i' %(k+1))
                        header=[np.array(a),np.array(b)]
                        c=[]
                        d=[]
                        for m in range(len(coordination[4][i][j][k])):
                            c.append('Diff. %s' %(combination[4][i][3]))
                            d.append('%i' %(k+1))
                        coord_id=[np.array(c),np.array(d)]
                        data=np.zeros((len(coordination[4][i][j][k]),x))
                        for m in range(len(coordination[4][i][j][k])):
                            for n in range(len(coordination[4][i][j][k][m])):
                                data[m][n]=coordination[4][i][j][k][m][n]
                        df=pd.DataFrame(data,index=coord_id,columns=header)
                        df_string = df.to_string(header=True, index=True)
                        f.write(df_string)
                        f.write('\n')
    # Quintuple Coordination
    if len(coordination)>5:
        f.write('\n')
        f.write('\n')
        f.write('Coordination with 5 types of molecules: \n')
        for i in range(len(coordination[5])):
            if len(coordination[5][i])==0:
                f.write('%s coordinating with %s, %s, %s, %s and %s: \n' %(reference,combination[5][i][0],combination[5][i][1],combination[5][i][2],combination[5][i][3],combination[5][i][4]))
                f.write('\n')
                f.write('\n')
                continue
            else:
                for j in range(len(coordination[5][i])):
                    for k in range(len(coordination[5][i][j])):
                        f.write('\n')
                        f.write('%s coordinating with %i %s, %i %s, %i %s, %s and %s: \n' %(reference,i+1,combination[5][0][0],j+1,combination[5][0][1],k+1,combination[5][0][2],combination[5][0][3],combination[5][0][4]))
                        a=[]
                        b=[]
                        x=[]
                        for m in range(len(coordination[5][i][j][k])):
                            x.append(len(coordination[5][i][j][k][m]))
                        if len(x)!=0:
                            x=max(x)
                        for m in range(x):
                            a.append('Diff. %s' %(combination[5][0][3]))
                            b.append('%i' %(k+1))
                        header=[np.array(a),np.array(b)]
                        c=[]
                        d=[]
                        for m in range(len(coordination[5][i][j][k])):
                            c.append('Diff. %s' %(combination[5][0][4]))
                            d.append('%i' %(k+1))
                        coord_id=[np.array(c),np.array(d)]
                        data=np.zeros((len(coordination[5][i][j][k]),x))
                        for m in range(len(coordination[5][i][j][k])):
                            for n in range(len(coordination[5][i][j][k][m])):
                                data[m][n]=coordination[5][i][j][k][m][n]
                        df=pd.DataFrame(data,index=coord_id,columns=header)
                        df_string = df.to_string(header=True, index=True)
                        f.write(df_string)
                        f.write('\n')
    f.close()


# Function that sums each tipe of coordination and writes the results in the Coord.txt file.
def write_summary(coordination,combination):
    total=suma_coordination(coordination)
    f=open('Coord.txt','a')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('Isolated %s: %i (%.2f' %(reference,coordination[0],100*coordination[0]/total)+' %) \n')
    f.write('\n')
    for i in range(len(coordination[1])):
        f.write('%s coordinated with %s: %i (%.2f' %(reference,combination[1][i],sum(coordination[1][i]),100*sum(coordination[1][i])/total) +' %) \n')
    f.write('\n')
    for i in range(len(coordination[2])):
        s=0
        for j in range(len(coordination[2][i])):
            s=s+sum(coordination[2][i][j])
        f.write('%s coordinated with %s and %s: %i (%.2f' %(reference,combination[2][i][0],combination[2][i][1],s,100*s/total) +' %) \n')
    if len(coordination)>3:
        f.write('\n')
        for i in range(len(coordination[3])):
            s=0
            for j in range(len(coordination[3][i])):
                for k in range(len(coordination[3][i][j])):
                    s=s+sum(coordination[3][i][j][k])
            f.write('%s coordinated with %s, %s and %s: %i (%.2f' %(reference,combination[3][i][0],combination[3][i][1],combination[3][i][2],s,100*s/total) +' %) \n')
    if len(coordination)>4:
        f.write('\n')
        for i in range(len(coordination[4])):
            s=0
            for j in range(len(coordination[4][i])):
                for k in range(len(coordination[4][i][j])):
                    for m in range(len(coordination[4][i][j][k])):
                        s=s+sum(coordination[4][i][j][k][m])
            f.write('%s coordinated with %s, %s, %s and %s: %i (%.2f' %(reference,combination[4][i][0],combination[4][i][1],combination[4][i][2],combination[4][i][3],s,100*s/total) +' %) \n')
    if len(coordination)>5:
        f.write('\n')
        s=0
        for i in range(len(coordination[5])):
            for j in range(len(coordination[5][i])):
                for k in range(len(coordination[5][i][j])):
                    for m in range(len(coordination[5][i][j][k])):
                        s=s+sum(coordination[5][i][j][k][m])
        f.write('%s coordinated with %s, %s, %s, %s and %s: %i (%.2f' %(reference,combination[5][i][0],combination[5][i][1],combination[5][i][2],combination[5][i][3],combination[5][i][4],s,100*s/total) +' %) \n')
    f.close()

# Function that sums all the coordinations to find the total.
def suma_coordination(coordination):
    total=0
    total=total+coordination[0]
    if len(coordination)>2:
        for i in range(len(coordination[1])):
            total=total+sum(coordination[1][i])
        for i in range(len(coordination[2])):
            for j in range(len(coordination[2][i])):
                total=total+sum(coordination[2][i][j])
    if len(coordination)>3:
        for i in range(len(coordination[3])):
            for j in range(len(coordination[3][i])):
                for k in range(len(coordination[3][i][j])):
                    total=total+sum(coordination[3][i][j][k])
    if len(coordination)>4: 
        for i in range(len(coordination[4])):
            for j in range(len(coordination[4][i])):
                for k in range(len(coordination[4][i][j])):
                    for m in range(len(coordination[4][i][j][k])):
                        total=total+sum(coordination[4][i][j][k][m])  
    if len(coordination)>5:     
        for i in range(len(coordination[5])):
            for j in range(len(coordination[5][i])):
                for k in range(len(coordination[5][i][j])):
                    for m in range(len(coordination[5][i][j][k])):
                        total=total+sum(coordination[5][i][j][k][m])       
    return(total)


# Function that plots the results in a Venn diagram (only available for 2 or 3 observables).
def plot_venn(coordination,observable):
    total=suma_coordination(coordination)
    l=len(coordination)
    if l==3:
        indi=np.zeros((len(coordination[1])))
        for j in range(len(coordination[1])):
            indi[j]=sum(coordination[1][j])
        doble=0
        for j in range(len(coordination[2])):
            for k in range(len(coordination[2][j])):
                doble=doble+sum(coordination[2][j][k])
        labels=observable
        plt.figure()        
        subset=float("{:.1f}".format(indi[0]/total*100)),float("{:.1f}".format(indi[1]/total*100)),float("{:.1f}".format(doble/total*100))      
        v2=vplt.venn2_unweighted(subsets=subset,set_labels=labels)
        for text in v2.set_labels:
            text.set_fontsize(30)
            text.set_fontweight('bold')
        for text in v2.subset_labels:
            text.set_fontsize(26)
            text.set_fontweight('bold')
        v2.get_patch_by_id('10').set_color('r')
        v2.get_patch_by_id('01').set_color('b')
        v2.get_patch_by_id('11').set_color('m')
        v2.get_patch_by_id('10').set_alpha(0.8)
        v2.get_patch_by_id('01').set_alpha(0.8)
        v2.get_patch_by_id('11').set_alpha(0.8)
        plt.show()
    if l==4:
        indi=np.zeros((len(coordination[1])))
        for j in range(len(coordination[1])):
            indi[j]=sum(coordination[1][j])
        doble=np.zeros((len(coordination[2])))
        for j in range(len(coordination[2])):
            d=0
            for k in range(len(coordination[2][j])):
                d=d+sum(coordination[2][j][k])
            doble[j]=d
        triple=0
        for j in range(len(coordination[3][0])):
            for k in range(len(coordination[3][0][j])):
                for m in range(len(coordination[3][0][j][k])):
                    triple=triple+coordination[3][0][j][k][m]
        labels=observable
        plt.figure()    
        subset=float("{:.1f}".format(indi[0]/total*100)),float("{:.1f}".format(indi[1]/total*100)),float("{:.1f}".format(doble[0]/total*100)),float("{:.1f}".format(indi[2]/total*100)),float("{:.1f}".format(doble[1]/total*100)),float("{:.1f}".format(doble[2]/total*100)),float("{:.1f}".format(triple/total*100))              
        v3=vplt.venn3_unweighted(subsets=subset,set_labels=labels)
        for text in v3.set_labels:
            text.set_fontsize(30)
            text.set_fontweight('bold')
        for text in v3.subset_labels:
            text.set_fontsize(26)
            text.set_fontweight('bold')
        v3.get_patch_by_id('100').set_color('r')
        v3.get_patch_by_id('010').set_color('b')
        v3.get_patch_by_id('001').set_color('g')
        v3.get_patch_by_id('110').set_color('m')
        v3.get_patch_by_id('101').set_color('gold')
        v3.get_patch_by_id('011').set_color('c')
        v3.get_patch_by_id('111').set_color('w')
        v3.get_patch_by_id('100').set_alpha(0.8)
        v3.get_patch_by_id('010').set_alpha(0.8)
        v3.get_patch_by_id('001').set_alpha(0.8)
        v3.get_patch_by_id('110').set_alpha(0.8)
        v3.get_patch_by_id('101').set_alpha(0.8)
        v3.get_patch_by_id('011').set_alpha(0.8)
        v3.get_patch_by_id('111').set_alpha(0.8)
    plt.show()


# Function that plots the probability map if we only have two observables.
def plot_prob_map(coordination,observable):
    total=suma_coordination(coordination)
    size=0
    for i in range(len(coordination[1])):
        if len(coordination[1][i])>size:
            size=len(coordination[1][i])
    if len(coordination[2])>size:
        size=len(coordination[2])
    for i in range(len(coordination[2])):
        if len(coordination[2][i])>size:
            size=len(coordination[2][i])
    size=size+1
    prob_map=np.zeros((size,size))
    for i in range(len(coordination[1][0])):
        prob_map[0][i+1]=100*coordination[1][0][i]/total
    for i in range(len(coordination[1][1])):
        prob_map[i+1][0]=100*coordination[1][1][i]/total
    for i in range(len(coordination[2][0])):
        for j in range(len(coordination[2][0][i])):
            prob_map[i+1][j+1]=100*coordination[2][0][i][j]/total
    limit=0
    for i in range(len(prob_map)):
        if max(prob_map[i])>limit:
            limit=max(prob_map[i])
    for i in range(len(atom_id)):
        if observable[0] in atom_id[i]:
            x_label=molecule_id[i]
    for i in range(len(atom_id)):
        if observable[1] in atom_id[i]:
            y_label=molecule_id[i]
    plt.figure()
    plt.imshow(prob_map,cmap='hot_r',vmin=0,vmax=limit,origin='lower')
    cbar=plt.colorbar()
    cbar.ax.set_title('Prob / %',fontsize=22)
    cbar.ax.tick_params(labelsize=22)
    plt.xlabel('# %s molecules' %(x_label),fontsize=22)
    plt.ylabel('# %s molecules' %(y_label),fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.show()



''' CODE '''

# Name of the input file.
print('\n')
entr=input('Trayectory File Name: ')
input_file=entr+'.txt'

# Ask for the reference and observables.
reference,observable=Type_of_Interaction()

# Get the information of the system.
time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms=System_Analysis()

# Ask for the initial and final time of analysis.
print('Simulation of %i ns composed by %i snapshots every %i ps\n' %(time[-1]/1000,len(time),time[1]-time[0]))
start_time=float(input('Initial Time (ns): '))*1000
end_time=float(input('Final Time (ns): '))*1000 
print('\n')

# Ask for the cutoff radius.
cutoff=float(input('Cutoff Distance (A): '))

# Create the variable coordination with all the posible combinations.
coordination=[0]
combination=['']
for i in range(1,len(observable)+1):
    pos_comb=np.math.factorial(len(observable))/(np.math.factorial(i)*np.math.factorial(len(observable)-i))
    combination.append([])
    coordination.append([])
    for j in range(int(pos_comb)):
        coordination[i].append([])
        if i>1:
            combination[i].append([])
    m=0
    for j in itertools.combinations(observable,i):
        for k in j:
            if i>1:
                combination[i][m].append(k)
            else:
                combination[i].append(k)
        m+=1

# Analysis of the trajectory.
# For each time step, create the hypercube and check the distances between reference and observables.
print('\n')
print('COMPUTING COORDINATION ANALYSIS:')
for i in range(len(time)):
    if time[i]<start_time:
        continue
    elif time[i]>end_time:
        break
    else:
        molecule_number=[]
        for j in range(len(observable)):
            molecule_number.append([])
        hyper_cube,ref_position=Create_Hyper_Cube(i)
        actual_ref=ref_position[0][3]
        for j in range(len(ref_position)):
            if ref_position[j][3]!=actual_ref:
                coordination=coordination_save(coordination,molecule_number)
                actual_ref=ref_position[j][3]
                molecule_number=[]
                for k in range(len(observable)):
                    molecule_number.append([])
            pos_ref=ref_position[j][0:3]
            m_n=Check_Observable(pos_ref,actual_ref,hyper_cube)
            for k in range(len(m_n)):
                for m in m_n[k]:
                    molecule_number[k].append(m)
        Progress_Bar(i,len(range(time.index(start_time),time.index(end_time)+1)))
        coordination=coordination_save(coordination,molecule_number)
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %' , end='\r')

# Write the results in the Coord.txt file.
write_coordination(coordination,combination)
write_summary(coordination,combination)

# Plot the results in a Venn diagram.
if len(coordination)<5:
    plot_venn(coordination,observable)

# Plot the probability map.
if len(coordination)==3:
    plot_prob_map(coordination,observable)