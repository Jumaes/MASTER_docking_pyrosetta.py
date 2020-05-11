# This script is supposed to make use of the notion, that to some extend, protein-protein interfaces are just not-yet connected domain interfaces within  a protein.
# Therefore, e.g. for two alpha helices or a beta sheet and a helix or ... in the PDB, there will be numerous examples of those secondary structures to have good packing interactions in the interior of a protein.
# Those conformations of the two motifs towards each other will likely be very well packable, while a slighltly different conformation might just be a bit too far away, too close, too tilted.
# This usually is taking care of in design algorithms by a docking step, but the question is, if empirical examples are not the better idea.

#Step1 : Take input conformation of two structural motifs, which are a simplified depiction of two full proteins. Take a number of conformations, which should be tested (thinkinkg about how much can be handled in reasonable time). Then split this number over the 6 degrees of freedom to calculate many steps in a predefined range are going to be cheked.
# These inputs are then used by prosetta, to construct poses going through all the 6 degrees of freedom of the chain-chain jump and changing them stepwise, exhausting the given design space. Write out pdb files with names indicating the values of all 6 DF.

#Step2: Use MASTER to search for hits within a very narrow RMSD for every one of the pdb files from Step1. The Library against which is searched will be a non-redundant subset of the PDB with roughly 40 000 structures. Write out a number for the hits for each structure (potentially several different RMSD cutoffs for this?).

#Step3: Combined all the numbers of hits with the different parameters of DF into a massive pandas DataFrame and produce some maybe graphical output.

import math
from math import *
from pyrosetta import *
from rosetta import *
from rosetta.numeric import xyzVector_double_t
from rosetta.numeric import z_rotation_matrix_degrees
from pyrosetta.rosetta.protocols.geometry import centroids_by_jump

import numpy as np
import pandas as pd

import glob

import subprocess

init()

#Step 1:


#Input Parameters
input_pdb = 'NCFus_model2_parts.pdb'  #File with input conformation of the minimized structural motifs.
input_number_of_conformations = 5000 #Number of different conformations that should be generated and tested. This will be crucial to determine how long it will take.
subfolder = 'pdb_models/'
master_program_folder = '/media/tezcanlab/Data/Julian/programs/master/master-v1.3.1/'
query_pds_subfolder = 'query_pds_files/'
query_library_path = '/media/tezcanlab/Data/PDBdata/70ident_20200423/master_pds_database/'

#Note: For the moment the axes and angles will just be applied to the jump and they are most likely not global x,y,z coordinates or angles, but might be local. If control over axes/angles is desired, might need to be taken care off later.
# If a DF should not be designable, just put 0.
dev_x_translate = 6 #Total translation in x axis in A. Together with number_of_conformations this will also determine the stepsize for this DF. 4 means 2A in each direction from starting conformation.
dev_y_translate = 6
dev_z_translate = 6

dev_alpha_rotate = 10 # Similar to translate, just this is the total angle in degree to be screened for this DF.
dev_beta_rotate = 10
dev_gamma_rotate = 10

jump_to_move = 1

# Now calculating a few numbers, which derive from that.

designable_DF = 0
for a in [dev_x_translate,dev_y_translate,dev_z_translate,dev_alpha_rotate,dev_beta_rotate,dev_gamma_rotate]:
    if a != 0: designable_DF += 1 # Just counts how many DF are supposed to be designable.

steps_per_DF = int(math.floor(input_number_of_conformations ** (1. /designable_DF))) # This basically takes the root to the power of designable_DF and then returns the integer always rounding down to see with the given number_of_conformations, how many steps can be done per DF.

number_of_conformations = steps_per_DF ** designable_DF

print ('Requested %s conformations to be checked using %s degrees of freedom. %s steps per degree of freedom can be checked giving a total of %s conformations. With %s steps per DF, %s conformations would need to be checked.' %(input_number_of_conformations, designable_DF,steps_per_DF, number_of_conformations,steps_per_DF+1,(steps_per_DF+1) ** designable_DF ))

#Now we need to make the different conformations as pdb files using pyrosetta:
start = Pose()
pose_from_file(start, input_pdb)
jumppoint1 = start.pdb_info().pdb2pose('A',105)  # Residue one first chain to jump from to the other chain.
jumppoint2 = start.pdb_info().pdb2pose('E',56) # Reside on second Chain to jump to.
for i in range(2,start.__len__()):
    if start.residue(i).is_upper_terminus():
        end1 = i
        start2 = i+1
        break
#print ('end1:' +str(end1))
#print ('start2:' +str(start2))
#print ('jumppoint1:' +str(jumppoint1))
#print ('jumppoint2:' +str(jumppoint2))

ft=FoldTree()
ft.add_edge(jumppoint1,jumppoint2, 1)
ft.add_edge(jumppoint1,1,-1)
if jumppoint1 != end1:
    ft.add_edge(jumppoint1,end1,-1)

ft.add_edge(jumppoint2,start.__len__(), -1)
if jumppoint2 != start2:
    ft.add_edge(jumppoint2,start2, -1)
#print (ft)
if ft.check_fold_tree(): #Just making sure this is old done when foldtree is ok, because otherwise it immediately segfaults.
    start.fold_tree(ft)

workpose = Pose()
workpose.assign(start)
# Now need to apply the DF changes to the pose.

def one_or_zero(numberin):
    if numberin.__nonzero__():
        output = 1
    else:
        output = 0
    return output

# Sounds a bit insane, but in the latest rosetta version the functions to make a rotation matrix for rotation around x,y or z axis from a degree value are not there anymore. So make that myself after looking it up.
def rotation_matrix_from_alpha_degree(angle):
    radangle= radians(angle)
    matrix = rosetta.numeric.xyzMatrix_double_t.cols(1,0,0,0,cos(radangle),-sin(radangle),0,sin(radangle),cos(radangle))
    return matrix

def rotation_matrix_from_beta_degree(angle):
    radangle= radians(angle)
    matrix = rosetta.numeric.xyzMatrix_double_t.cols(cos(radangle),0,sin(radangle),0,1,0,-sin(radangle),0,cos(radangle))
    return matrix

def rotation_matrix_from_gamma_degree(angle):
    radangle= radians(angle)
    matrix = rosetta.numeric.xyzMatrix_double_t.cols(cos(radangle),-sin(radangle),0,sin(radangle),cos(radangle),0,0,0,1)
    return matrix

def matrixmultiplication(m1,m2):
	#Yes, seriously. The rosetta.numeric.xyzMatrix_double_t object has no matrix multiplication function ?!?!?!?!?!?! It should have, but it doesn't.
	#So just converting that object into a list of lists, then using numpy to do the matrix multiplication and then convert that back into that stupid object.
	#I mean, really...!?!?!?
	#Ah and also that object doesn't support indexing, so that elements have to be called with matrix.row() and that vector object can then be indexed like any normal python object.
	list1=[['a']*3,['a']*3,['a']*3]
	for i in range(0,3):
		for k in range(0,3):
			list1[i][k]=m1.row(i+1)[k]
	list2=[['a']*3,['a']*3,['a']*3]
	for i in range(0,3):
		for k in range(0,3):
			list2[i][k]=m2.row(i+1)[k]
	#Actual multiplcation is done by using .dot of numpy.
	p=np.dot(list1,list2)
    #And transform it back into that ridiculous rosetta.numeric.xyzMatrix_double_t object
	product=rosetta.numeric.xyzMatrix_double_t.cols(p[0][0],p[0][1],p[0][2],p[1][0],p[1][1],p[1][2],p[2][0],p[2][1],p[2][2])
	return product


def rotate_pose(inputpose,inputjump,rota,rotb,rotg):
    functionpose = Pose()
    functionpose.assign(inputpose)
    # Relatively complicated movers.
    # We need the jumpobject corresponding to the given jumpnumber.
    flexible_jump=functionpose.jump(inputjump)
    #We need the stubs of the protein on either side if the jump, if we wanna move them.
    upstream_stub=functionpose.conformation().upstream_jump_stub(inputjump)
    downstream_stub=functionpose.conformation().downstream_jump_stub(inputjump)
    #Need to find the rotation centres for rotating the upstream and downstream stubs respectively; they also depend on the jump and the workpose.
    upstream_dummy=xyzVector_double_t() #Just need to define that these are vectors of that weird type.
    downstream_dummy=xyzVector_double_t()
    centroids_by_jump(functionpose,inputjump,upstream_dummy,downstream_dummy) #This function does the job for us.
    # After testing a rotation by axis function, which didn't really work well, now back to Maxtrix rotation from the Negative Docking script.
    # Except rosetta.numeric seems to not have the functions z_rotation_matrix_degrees_double_t anymore. z_rotation_matrix_degrees seems to be the thing now.
    # But there are no equivalents for x and y. So just made my own.
    matrixA=rotation_matrix_from_alpha_degree(rota) #This is  the pyrosetta function. The original mover used a c++ function called numeric.z_rotation_matrix_degrees, which doesn't seem to exist in pyrosetta. Hope this here is the replacement.
    matrixB=rotation_matrix_from_beta_degree(rotb)
    matrixG=rotation_matrix_from_gamma_degree(rotg)
    # And now using my home-made matrix-multiplication function to get the total resulting rotation matrix.
    rotation_matrix=matrixmultiplication(matrixmultiplication(matrixA,matrixB),matrixG)
    flexible_jump.rotation_by_matrix(upstream_stub, upstream_dummy,  rotation_matrix)
    functionpose.set_jump(inputjump, flexible_jump)
    centroids_by_jump(functionpose,inputjump,upstream_dummy,downstream_dummy) #seems necessary to do that after the rotation. Maybe updates something.
    return functionpose


def translate_pose(inputpose,inputjump,movex,movey,movez):
    functionpose = Pose()
    functionpose.assign(inputpose)
    # Relatively complicated movers.
    # We need the jumpobject corresponding to the given jumpnumber.
    flexible_jump=functionpose.jump(inputjump)
    upstream_dummy=xyzVector_double_t() #Just need to define that these are vectors of that weird type.
    downstream_dummy=xyzVector_double_t()
    centroids_by_jump(functionpose,inputjump,upstream_dummy,downstream_dummy)
    # Had lots of problems with the translation by axis function. Instead use get_translation as vector for this jump, then just add the intended translation, and then set the sum as new translation vector for this jump.
    # This also allows for one operation at the same time for all three directions. Much better.
    before_translate = flexible_jump.get_translation()
    translate_vec = xyzVector_double_t(movex, movey, movez)
    after_translate = before_translate + translate_vec
    #print ('New translation vector for the jump is: ')
    #print (after_translate)
    flexible_jump.set_translation(after_translate)
    functionpose.set_jump(inputjump, flexible_jump ) # This is actually applying the jump to the pose.
    #print (functionpose.residue(20))
    centroids_by_jump(functionpose,inputjump,upstream_dummy,downstream_dummy) #Maybe necessary to do that here. Maybe updates something.
    return functionpose





def make_steps(translate,steps):
    if translate == 0:
        return [0] # Otherwise even for a DF with 0 rotation/translation, it still is going to make a list with several members all being 0.
    else:
        return [(-float(translate)/2)+float(translate)/(float(steps)-1)*number for number in range(0,steps)]

x_steps = make_steps(dev_x_translate,steps_per_DF)
y_steps = make_steps(dev_y_translate,steps_per_DF)
z_steps = make_steps(dev_z_translate,steps_per_DF)
a_steps = make_steps(dev_alpha_rotate,steps_per_DF)
b_steps = make_steps(dev_beta_rotate, steps_per_DF)
g_steps = make_steps(dev_gamma_rotate, steps_per_DF)

try:
    os.mkdir(subfolder)
except:
    pass

structure_dict = {}

print ('Now preparing structures with modified geometry.')
counter = 0
for x in x_steps:
    for y in y_steps:
        for z in z_steps:
            trans_pose = translate_pose(workpose,jump_to_move,x,y,z)
            #trans_pose.dump_pdb('%s%s_%s.pdb' %(subfolder, input_pdb[:-4],counter))
            for a in a_steps:
                for b in b_steps:
                    for g in g_steps:
                        rota_trans_pose = rotate_pose(trans_pose,jump_to_move,a,b,g)
                        rota_trans_pose.dump_pdb('%s%s_%s.pdb' %(subfolder, input_pdb[:-4],counter))
                        structure_dict[counter] = [x,y,z,a,b,g]
                        with open('%s%s_%s.pdb' %(subfolder, input_pdb[:-4],counter),'a') as fout:
                            fout.write('REMARK %f %f %f %f %f %f \n' %(x,y,z,a,b,g))
                        counter += 1


#Step2: Use MASTER to search for hits within a very narrow RMSD for every one of the pdb files from Step1. The Library against which is searched will be a non-redundant subset of the PDB with roughly 40 000 structures. Write out a number for the hits for each structure (potentially several different RMSD cutoffs for this?).
query_pdbs_list = glob.glob(subfolder+'*.pdb') # Get a list of all actually created pdb files.
try:
    os.mkdir(query_pds_subfolder) # As usual make subfolder for the pds query files, that master needs to actually then search.
except:
    pass

for file in query_pdbs_list: # Next line makes a query pds file from every pdb file and puts the pds files into the new subfolder. That's pretty fast with several files per seconds.
    os.system('%screatePDS --type query --pdb %s --pds %s%s.pds' %(master_program_folder,file, query_pds_subfolder,file[:-4].split('/')[-1]))

query_pds_list = glob.glob(query_pds_subfolder+'*.pds')
library_file_list = glob.glob('%s*.pds' %(query_library_path))
with open('pds_library.txt','w') as lib_out:
    for x in range(0,len(library_file_list),10):
        lib_out.write(library_file_list[x] + '\n')

def extract_number_of_hits(stream_in):
    outlist = [item.split()[0] for item in stream_in.split('\n')[:-1]]
    Hits = []
    for item in outlist:
        try:
            rmsd = float(item)
            Hits.append(rmsd)
        except:
            pass
    return (Hits)

def get_search_time(stream_in):
    outlist = stream_in.split('\n')[:-1]
    for item in outlist:
        if 'Search ' in item:
            break
    return item


hitlist = []
print('Starting acutal search through master library.')
searchcounter = 0
for file in query_pds_list:
    print ('Starting search for query sequence %s' %(file))
    struct_number = int(file[:-4].split('_')[-1])
    search = subprocess.Popen(['%smaster' %(master_program_folder), '--query', '%s' %(file), '--targetList', 'pds_library.txt', '--rmsdCut', '1.2'], stdout = subprocess.PIPE) #optional --bbRMSD to search full backbone RMSD vs just C alpha
    searchcounter += 1
    out,errors = search.communicate()
    print get_search_time(out) #This just prints out the search time. More interesting for debugging. Sacrifice for speed later.
    try:
        outhits = extract_number_of_hits(out)
        hitlist.append([len(outhits),np.array(outhits).mean(),struct_number])
        hitlist[-1].extend(structure_dict[struct_number])
    except:
        hitlist.append([0,0,struct_number])
        hitlist[-1].extend(structure_dict[struct_number])
    if divmod(searchcounter,100)[1] == 0: #This spits out some partial solutions every 100 searches. Just to get a bit of an output while it is running.
        df = pd.DataFrame(hitlist, columns=['hits','RMSD','number','x','y','z','alpha','beta','gamma']).set_index('number')
        df.to_csv('hitlist.csv')

df = pd.DataFrame(hitlist, columns=['hits','RMSD','number','x','y','z','alpha','beta','gamma']).set_index('number')
df.to_csv('hitlist.csv')



#TODO: Need to take the system output from master, then throw away all until 'Search completed' then count the lines until 'Output completed'
# Or keep those lines, split at whitespace, take 0 element and make mean and store number of lines plus mean for each structure, that was analysed. Maybe spit out a results table for every 100 strctures, just to not just lose the results.
