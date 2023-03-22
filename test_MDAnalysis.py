import MDAnalysis as mda
import numpy as np

negative = {"ASP", "GLU"}
positive = {"HIS", "ARG", "LYS"}

charged = negative | positive

cutoff = 7

mutationsA = []
mutationsB = []

u = mda.Universe("6ta4.pdb")

chainA = u.select_atoms("protein and segid A")
chainB = u.select_atoms("protein and segid B")
chainK = u.select_atoms("protein and segid K")

centroidA = chainA.center_of_geometry(compound="residues")
centroidB = chainB.center_of_geometry(compound="residues")
centroidK = chainK.center_of_geometry(compound="residues")


for res, centroid in zip(chainA.residues, centroidA):
    resname = res.resname
    resid = res.resid

    for res2, centroid2 in zip(chainK.residues, centroidK):
        resname2 = res2.resname
        resid2 = res2.resid

        dist = np.linalg.norm(centroid - centroid2)

        if (dist < cutoff) and (resname in charged) and (resname2 not in charged):
            mutationsA.append([dist, "chainA", resname, resid, "chainK", resname2, resid2])


for res, centroid in zip(chainB.residues, centroidB):
    resname = res.resname
    resid = res.resid

    for res2, centroid2 in zip(chainK.residues, centroidK):
        resname2 = res2.resname
        resid2 = res2.resid

        dist = np.linalg.norm(centroid - centroid2)

        if (dist < cutoff) and (resname in charged) and (resname2 not in charged):
            mutationsB.append([dist, "chainB", resname, resid, "chainK", resname2, resid2])

print(mutationsA, "\n", mutationsB)