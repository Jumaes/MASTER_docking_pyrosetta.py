# This script is supposed to make use of the notion, that to some extend, protein-protein interfaces are just not-yet connected domain interfaces within  a protein.
# Therefore, e.g. for two alpha helices or a beta sheet and a helix or ... in the PDB, there will be numerous examples of those secondary structures to have good packing interactions in the interior of a protein.
# Those conformations of the two motifs towards each other will likely be very well packable, while a slighltly different conformation might just be a bit too far away, too close, too tilted.
# This usually is taking care of in design algorithms by a docking step, but the question is, if empirical examples are not the better idea.

#Step1 : Take input conformation of two structural motifs, which are a simplified depiction of two full proteins. Take a number of conformations, which should be tested (thinkinkg about how much can be handled in reasonable time). Then split this number over the 6 degrees of freedom to calculate many steps in a predefined range are going to be cheked.
# These inputs are then used by prosetta, to construct poses going through all the 6 degrees of freedom of the chain-chain jump and changing them stepwise, exhausting the given design space. Write out pdb files with names indicating the values of all 6 DF.

#Step2: Use MASTER to search for hits within a very narrow RMSD for every one of the pdb files from Step1. The Library against which is searched will be a non-redundant subset of the PDB with roughly 40 000 structures. Write out a number for the hits for each structure (potentially several different RMSD cutoffs for this?).

#Step3: Combined all the numbers of hits with the different parameters of DF into a massive pandas DataFrame and produce some maybe graphical output.

import math

#Step 1:


#Input Parameters
input_pdb = ''  #File with input conformation of the minimized structural motifs.
input_number_of_conformations = 10000 #Number of different conformations that should be generated and tested. This will be crucial to determine how long it will take.

#Note: For the moment the axes and angles will just be applied to the jump and they are most likely not global x,y,z coordinates or angles, but might be local. If control over axes/angles is desired, might need to be taken care off later.
# If a DF should not be designable, just put 0.
dev_x_translate = 4 #Total translation in x axis in A. Together with number_of_conformations this will also determine the stepsize for this DF. 4 means 2A in each direction from starting conformation.
dev_y_translate = 4
dev_z_translate = 4

dev_alpha_rotate = 10 # Similar to translate, just this is the total angle in degree to be screened for this DF.
dev_beta_rotate = 10
dev_gamma_rotate = 10

designable_DF = 0
for a in [dev_x_translate,dev_y_translate,dev_z_translate,dev_alpha_rotate,dev_beta_rotate,dev_gamma_rotate]:
    if a != 0: designable_DF += 1 # Just counts how many DF are supposed to be designable.

steps_per_DF = int(math.floor(input_number_of_conformations ** (1. /designable_DF))) # This basically takes the root to the power of designable_DF and then returns the integer always rounding down to see with the given number_of_conformations, how many steps can be done per DF.

number_of_conformations = steps_per_DF ** designable_DF

print ('Requested %s conformations to be checked using %s degrees of freedom. %s steps per degree of freedom can be checked giving a total of %s conformations. With %s steps per DF, %s conformations would need to be checked.' %(input_number_of_conformations, designable_DF,steps_per_DF, number_of_conformations,steps_per_DF+1,(steps_per_DF+1) ** designable_DF ))
