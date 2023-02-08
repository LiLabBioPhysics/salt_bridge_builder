from structure import *
from distances import *
from filters import *
import pandas as pd




"""This is a file which shows how to use my code.
The basic compoenents of my code so far are:
1. A structure file (structure.py) defining proteins and/or trajectories
2. A distance file (distances.py) which defines the distances that can
    be applied to different residues/atoms.
3. A filter file (filters.py) which can be used to filter results by whether
    or not the residue pair might be a good candiate for mutation. Optionally,
    you can specify a dont_consider file to ignore certain results from the filter."""


#First you create a trajectory object, traj1.
#Note, you can read entire trajectories (like with traj2) but this takes a lot more time (X1000) so I commented it out.
traj1 = TRAJECTORY("frame1.pdb")
#traj2 = TRAJECTORY("150mM_2000_3000.pdb")

#Once you have a trajectory, you need to specify what frame to read. I am reading frame 0 here. Doing so outputs a complex.
cmplx1 = traj1[0]

#Once you have a complex, you can specify what chains to apply a distance to.
#For example, the following line uses the centroid distance (centroid_dist)
# defined in distances to find all residues which are 7 angstroms apart 
centroid1  = within_dist(cmplx1, "A", "K", 7, centroid_dist)

#From these results, you can filter those pairs which are likely to be good
#candidates for mutations. Specifically, we care about those pairs for which one
#residue is charged and the other is neutral. You can optionally
#input a dont_consider file to ignore certain results (like in filtered_centroid2 below).
#If you do so, make sure the dont_consider file has exactly the same columnames as the 
#example.

filtered_centroid1 = filter_charged_neutral(centroid1)
filtered_centroid2 = filter_charged_neutral(centroid1, dont_consider="dont_consider.csv")

print(filtered_centroid1)
print(filtered_centroid2)

#Because the output from filter_charged_neutral is a pandas dataframe,
#you can easily save the dataframe to a csv file.
filtered_centroid1.to_csv("total_results.csv", index=False)
filtered_centroid2.to_csv("reduced_results.csv", index=False)