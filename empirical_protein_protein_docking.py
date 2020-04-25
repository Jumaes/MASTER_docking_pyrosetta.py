# This script is supposed to make use of the notion, that to some extend, protein-protein interfaces are just not-yet connected domain interfaces within  a protein.
# Therefore, e.g. for two alpha helices or a beta sheet and a helix or ... in the PDB, there will be numerous examples of those secondary structures to have good packing interactions in the interior of a protein.
# Those conformations of the two motifs towards each other will likely be very well packable, while a slighltly different conformation might just be a bit too far away, too close, too tilted.
# This usually is taking care of in design algorithms by a docking step, but the question is, if empirical examples are not the better idea.

#Step1 : Take input conformation of two structural motifs, which are a simplified depiction of two full proteins. Take a number of conformations, which should be tested (thinkinkg about how much can be handled in reasonable time). Then split this number over the 6 degrees of freedom to calculate many steps in a predefined range are going to be cheked.
# These inputs are then used by prosetta, to construct poses going through all the 6 degrees of freedom of the chain-chain jump and changing them stepwise, exhausting the given design space. Write out pdb files with names indicating the values of all 6 DF.

#Step2: Use MASTER to search for hits within a very narrow RMSD for every one of the pdb files from Step1. The Library against which is searched will be a non-redundant subset of the PDB with roughly 40 000 structures. Write out a number for the hits for each structure (potentially several different RMSD cutoffs for this?).

#Step3: Combined all the numbers of hits with the different parameters of DF into a massive pandas DataFrame and produce some maybe graphical output.
