'''
    Author: Sergio Rodriguez Peña
'''



###################################################################################################
#                                                                                                 #
#   Code that takes as input the trajectory of a MD simulation in .pdb format and extracts the    #
#   important information for the future analyses in a file with the same name and .txt format.   #
#                                                                                                 #
#   Structure of the output file:                                                                 #
#   #_atom   ID_atom   #_molecule   ID_molecule  X   Y   Z                                        #
#                                                                                                 #
#   The .txt output of this code serves as input for the rest of analyses codes.                  #
#                                                                                                 #
###################################################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import os



''' DEFINED FUNCTIONS '''

# Function that identifies the different types of molecules and atoms of the system and 
# change the name of the atoms if they are repeated by adding letters of the molecule.
def System_Identification():
    molecule_id=[]
    atom_id=[]
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if 't=' in l and float(l[l.index('t=')+1])!=0:
                break
            if l[0]=='ATOM':
                mol=l[3][:4]
                if mol not in molecule_id:
                    molecule_id.append(mol)
                    atom_id.append([])
                at=l[2]
                if at not in atom_id[molecule_id.index(mol)]:
                    atom_id[molecule_id.index(mol)].append(at)
    atom_id_mod=[]
    for i in range(len(molecule_id)):
        atom_id_mod.append([])
        for j in range(len(atom_id[i])):
            a=atom_id[i][j]
            b=0
            for k in range(i):
                while a in atom_id_mod[k]:
                    a=a+molecule_id[i][:b]
                    b+=1
            atom_id_mod[i].append(a)
    return(molecule_id,atom_id,atom_id_mod)

    
# Function that reads the input .pdb file, line by line, saves the important information and writes the 
# output .txt file with the same name.
def Process_PDB_to_TXT():
    with open(input_file,'r') as r, open(output_file,'w') as w:
        s=0
        m=0
        text=[]
        print('PROCESSING INPUT FILE:')
        for line in r:
            s+=1
            l=line.split()
            if l[0]=='TITLE':
                time=float(l[l.index('t=')+1])
                text.append('t = %.3f ps \n' %(time))
            if l[0]=='CRYST1':
                box_size=float(l[1])
                text.append('Box lenght = %.3f A \n' %(box_size))
            if l[0]=='ATOM':
                atom_number,atom_name,molecule_name=l[1:4]
                molecule_name=molecule_name[:4]
                molecule_number=l[l.index('1.00')-4]
                x,y,z=l[l.index('1.00')-3:l.index('1.00')]
                atom_name=atom_id_mod[molecule_id.index(molecule_name)][atom_id[molecule_id.index(molecule_name)].index(atom_name)]
                text.append('%s %s %s %s %s %s %s\n' %(atom_number,atom_name,molecule_number,molecule_name,x[:-1],y[:-1],z[:-1]))
            if s==int(num_lines/1000):
                Progress_Bar(m*s,num_lines)
                s=0
                m+=1
                for i in text:
                    w.write(i)
                text=[]
        if text!=[]:
            for i in text:
                w.write(i)
        w.close()
        bar='|'+'#'*50+'|'
        print(f'%s 100.0' %(bar) +' %')


# Function that creates the progress bar shown in the terminal.
def Progress_Bar(progress,total):
    percent=100*(progress/float(total))
    bar='|'+'#'*int(percent/2)+'-'*(50-(int(percent/2)))+'|'
    print(f'%s %.1f' %(bar,percent) +' %' , end='\r')



''' CODE '''

# Name of the input .pdb file.
print('\n')
file_name=input('Input File Name (without .pdb extention): ')
input_file=file_name+'.pdb'
output_file=file_name+'.txt'

# Read the size and number of lines of the input file.
num_lines=sum(1 for _ in open(input_file))
size_input_file=os.stat(input_file).st_size
print('\n')
print('Size of Initial .pdb trajectory: %.2f MB \n' %(size_input_file/(1024*1024)))

# Analyze the composition of the system.
molecule_id,atom_id,atom_id_mod=System_Identification()

# Process the .pdb file to .txt.
Process_PDB_to_TXT()

# End the process printing the results.
size_output_file=os.stat(output_file).st_size
print('\n')
print('Size of Processed .txt trajectory: %.2f MB \n' %(size_output_file/(1024*1024)))
print('Space Reduction: %.2f' %(100*(size_input_file-size_output_file)/size_input_file) +' %')
