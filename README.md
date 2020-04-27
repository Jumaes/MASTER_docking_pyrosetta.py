# MASTER_docking_pyrosetta.py
Using Pyrosetta to create a large number of similar docking geometries as input protein dimer, then use master to find the interfaces, which are most likely designable according to examples from proteins in the pdb.

## Problems:
1. MASTER is only fast after the whole target database is loaded into the RAM once. With 1/10th of current library, search takes maybe 10 seconds after first loading of library. But with full library, probably too slow. 
=> Need to make a smaller library... maybe include only monomers or so.

2. At the moment didn't find any hits. Need to rerun the whole script with a part of the actualy cb562, which should give some hits, to test, if it works reasonably and with what RMSD.

3. Not sure what the output is going to be. Writing out too much sequences/structures, will just take way to much IO time. Ideally just capture them in a list or so.
