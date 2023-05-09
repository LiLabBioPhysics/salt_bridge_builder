import structure as st
import distances as ds
import pandas as pd

"""This is a file which shows how to use my code.
The basic compoenents of my code so far are:
1. A structure file (structure.py) defining proteins and/or trajectories
2. A distance file (distances.py) which defines the distances that can
    be applied to different residues/atoms/proteins."""

#You can create a trajectory object, traj which reads a multi-frame pdb file.
traj = st.TRAJECTORY("bimc_traj.pdb")

#You can index by a frame to create a complex.
frame0 = traj[0]
print(type(frame0))

#You can also create a complex from a single-frame pdb file if you don't have a trajectory.
cmplx = st.COMPLEX("6ta4.pdb")

#Once you read in a file, you can index to obtain chains, complexes, and residues

#Indexing chains A and B, respectively
chainA = cmplx["A"]
chainB = cmplx["B"]

#Indexing residues 1 and 233 from the coresponding chains
res1 = chainA[1]
res233 = chainB[233]

#Indexing the backbone nitrogen atoms from the coresponding residues
atomX = res1[0]
atomY = res233[0]

#Having saved various structure objects in memory, you can find different distances

#Atom distance
ds.atom_dist(atomX, atomY)

#Residue distances
ds.sidechain_centroid_dist(res1, res233)

#You can also find all of the pairwise distances from one chain to all the other chains
#for a cutoff and distance function of your choice.
info, dist = ds.all_dist(cmplx, "A", 8, ds.sidechain_centroid_dist)

#Also, you can find all existing salt bridges within a complex.
info, dist = ds.all_salt_bridges(cmplx, 3.2)

#Lastly, you can find all possible mutations for the uncharged residues of a chain which
#are within the vicinity of charged residue from another chain (measured by an specified residue distance).
#Optionally, you request that the charged residue not be invovled in any existing salt bridges.
info, dist = ds.find_putative_sb(cmplx, "K", 7, 3.2, ds.sidechain_centroid_dist, True)
colnames = ["chain1", "resname1", "resnum1", "chain2", "resname2", "resnum2"]
results_df = pd.DataFrame(data=info, columns=colnames)
results_df["sidechain_centroid_dist"] = dist
results_df.to_csv("results.csv", index=False)