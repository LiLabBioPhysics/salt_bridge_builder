<h1>
<p align="center">
    <img src="sbb.png" alt="Salt Bridge Builder Logo" width="500"/>
</p>
</h1>

[![Requires Python 3.12+](https://img.shields.io/badge/Python-3.12+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)

# Salt Bridge Builder
Salt Bridge Builder (SBB) is a python-based tool that identifies uncharged-to-charged mutations that result in interprotein salt bridges.

## Table of Contents
- [Dependencies](#dependencies)
- [Installation](#installation)
- [File Input and Output](#file-input-and-output)
- [Script Minimal Working Example](#script-minimal-working-example)
- [Terminal Minimal Working Example](#terminal-minimal-working-example)

## Dependencies
This project uses the following main dependencies:
*   [Python](https://www.python.org) (Version >= 3.12)
*   [Pip](https://pypi.org/project/pip/) (Version >= 26)
*   [NumPy](https://numpy.org) (Version >= 2.4.2)
*   [Pandas](https://pandas.pydata.org) (Version >= 3.0.1)
  
## Installation

The preferred method installing python and pip is via the [conda website](https://conda.io).

```bash
conda create -n sbb-env python=3.12 pip=26
conda activate sbb-env
pip install sbb
```
## File Input and Output
SBB can take single-frame or multi-frame PDB files as inputs. For PDB file specifications visit the [RCSB website](https://www.rcsb.org).

SBB outputs results as numpy arrays. These can be convereted to .csv files via pandas. See the examle in the [Script Minimal Working Example](#script-minimal-working-example) section.


## Script Minimal Working Example
SBB can be run within python.
```python
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
```

## Terminal Minimal Working Example
The same example can be run by invoking the test.py file from the terminal.
```bash
python test.py
```
