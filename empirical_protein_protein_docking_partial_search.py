#This is only the second part of the script, which is meant to be run after running the "setup" part.
#Therefore all the sampling structures should already be created as well as their search pds files.
#And the search library.txt file, which determines against which structures to search.
#This script can be run with two integers part and parts, allowing to split the full set of query structures into a number of subsets and running the search only on one of those subsets.
# If the script is just started for all subsets in parallel, then this is a very easy way at parallelizing the search.
#Of course this way, the output will be seperate csv files, but they can just super easily be combined when analyzing them by using Pandas.concat


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
import numpy as np
import pandas as pd
import sys
import glob

import subprocess

if len(sys.argv)=1:
    print ('This is the subscript, which only makes the actual search. So the setup script should have been run already.')
    print ('This script needs to be run with SCRIPTNAME X Y with Y being the number of parts the whole querylist shold be split in and X being the Xth part.')
    print ('So SCRIPTNAME 1 10 would split the querylist into 10 parts and run only on the 1st of them. Clearly, you would want to start parallel runs with SCRIPTNAME 2 10 .... etc until SCRIPTNAME 10 10')

part = sys.argv[1]
parts = sys.argv[2]


#Input Parameters
input_pdb = 'NCFus_model2_parts.pdb'  #File with input conformation of the minimized structural motifs.
input_number_of_conformations = 5000 #Number of different conformations that should be generated and tested. This will be crucial to determine how long it will take.
subfolder = 'pdb_models/'
master_program_folder = '/media/tezcanlab/Data/Julian/programs/master/master-v1.3.1/'
query_pds_subfolder = 'query_pds_files/'
query_library_path = '/media/tezcanlab/Data/PDBdata/70ident_20200423/master_pds_database/'

query_pds_list = glob.glob(query_pds_subfolder+'*.pds')

def get_number(infile):
    number = int(infile.split('_')[-1][:-4])
    return number

query_pds_list.sort(key = get_number)




def extract_rmsd_each_hit(stream_in):
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

def mastersearch(pdsin):
    print ('Starting search for query sequence %s' %(pdsin))
    struct_number = int(pdsin[:-4].split('_')[-1])
    # The following needed to make a subprocess, for which I can actually read out the STDOUT. Very strange way to parse the command as list of options, but that's how it has to be.
    search = subprocess.Popen(['%smaster' %(master_program_folder), '--query', '%s' %(pdsin), '--targetList', 'pds_library.txt', '--rmsdCut', '1.2'], stdout = subprocess.PIPE) #optional --bbRMSD to search full backbone RMSD vs just C alpha
    out,errors = search.communicate()
    print get_search_time(out) #This just prints out the search time. More interesting for debugging. Sacrifice for speed later.
    try:
        outhits = extract_rmsd_each_hit(out)
        outlist = [len(outhits),np.array(outhits).mean(),struct_number,structure_dict[struct_number]]
    except:
        outlist = [0,0,struct_number,structure_dict[struct_number]]
    return outlist

per_part = len(query_pds_list)/parts #Keep in mind that this is actually a floor division if two integers are passed; means it returnts the smallest integer, which still divides.
if not part == parts:
    partial_query_pds_list = query_pds_list[(part-1)*per_part:part*per_part]


hitlist = []
print('Starting actual search through master library.')
searchcounter = 0
for file in partial_query_pds_list:
    results_query_one_structure_vs_database = mastersearch(file)
    searchcounter += 1
    hitlist.append(results_query_one_structure_vs_database)
    if divmod(searchcounter,100)[1] == 0: #This spits out some partial solutions every 100 searches. Just to get a bit of an output while it is running.
        df = pd.DataFrame(hitlist, columns=['hits','RMSD','number','x','y','z','alpha','beta','gamma']).set_index('number')
        df.to_csv('hitlist.csv')




df = pd.DataFrame(hitlist, columns=['hits','RMSD','number','x','y','z','alpha','beta','gamma']).set_index('number')
df.to_csv('hitlist.csv')
