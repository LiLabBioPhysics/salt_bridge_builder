from structure import *
from constants import *
import numpy as np

"""Distance building blocks beyond this point"""
def get_centroid(res: RESIDUE) -> np.array:
    atomis = np.array(list(res.heavyatoms.keys()))
    if np.size(atomis) == 0:
        raise ValueError("Input a residue with at least one atom")
    
    x_tot = y_tot = z_tot = 0

    ac = res.atomcount

    for i in atomis:
        cur_atom = res[int(i)]
        x_tot += cur_atom.x
        y_tot += cur_atom.y
        z_tot += cur_atom.z
        
    cx = x_tot/ac
    cy = y_tot/ac
    cz = z_tot/ac

    return np.array([cx, cy, cz])

def get_sidechain_centroid(res: RESIDUE) -> np.array:
    atomis = np.array(list(res.heavyatoms.keys()))
    if np.size(atomis) == 0:
        raise ValueError("Input a residue with at least one atom")
    
    x_tot = y_tot = z_tot = 0

    ac = res.atomcount

    for i in atomis[atomis >= SIDE_CHAIN_START]:
        cur_atom = res[int(i)]
        x_tot += cur_atom.x
        y_tot += cur_atom.y
        z_tot += cur_atom.z
        
    cx = x_tot/(ac-SIDE_CHAIN_START)
    cy = y_tot/(ac-SIDE_CHAIN_START)
    cz = z_tot/(ac-SIDE_CHAIN_START)

    return np.array([cx, cy, cz])

"""Atom distance beyond this point"""

#Function returns the Euclidean distance between two atoms
def atom_dist(atom1:ATOM, atom2:ATOM) -> float:
    x1= atom1.x
    y1= atom1.y
    z1= atom1.z

    x2= atom2.x
    y2= atom2.y
    z2= atom2.z
    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(0.5)

"""Residue distances beyond this point"""

#Function returns all pairwise atom distancess for two residues. 
def pairwise_dist(res1: RESIDUE, res2: RESIDUE) -> np.ndarray:
    atom1is = np.array(list(res1.heavyatoms.keys()))
    atom2is = np.array(list(res2.heavyatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one atom")

    a1c = list(range(res1.atomcount))
    a2c = list(range(res2.atomcount))

    arr = np.zeros((len(a1c), len(a2c)))

    for i in atom1is:
        for j in atom2is:
            a_dist = atom_dist(res1[int(i)], res2[int(j)])
            arr[i,j] = a_dist

    return arr

#Function returns the minimum pairwise distance for two residues.
def minimum_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    pairwise_arr = pairwise_dist(res1, res2)
    return np.min(pairwise_arr)

#Function returns the alpha carbon distance between two residues.
def alpha_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    dist = atom_dist(res1[ALPHA_CARBON_POS], res2[ALPHA_CARBON_POS])
    return dist

#Function returns the minimum distance for all side chain atom pairs between two residues.
def sidechain_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    n1 = len(HEAVY_ATOM_POS[res1.name][SIDE_CHAIN_START:])
    n2 = len(HEAVY_ATOM_POS[res2.name][SIDE_CHAIN_START:])
    arr = np.empty((n1,n2))

    for i in range(n1):
        for j in range(n2):
            arr[i,j] = atom_dist(res1[i+SIDE_CHAIN_START], res2[j+SIDE_CHAIN_START])
    
    return np.min(arr)

#Function returns the Eulidean distance between the centroids of two residues.
def centroid_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    atom1is = np.array(list(res1.heavyatoms.keys()))
    atom2is = np.array(list(res2.heavyatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one atom")

    ac1 = res1.atomcount
    ac2 = res2.atomcount
 
    x1_tot = y1_tot = z1_tot = x2_tot = y2_tot = z2_tot = 0
        
    for i in atom1is:
        cur_atom = res1[int(i)]
        x1_tot += cur_atom.x
        y1_tot += cur_atom.y
        z1_tot += cur_atom.z
        
    c1x = x1_tot/ac1
    c1y = y1_tot/ac1
    c1z = z1_tot/ac1

    for i in atom2is:
        cur_atom = res2[int(i)]
        x2_tot += cur_atom.x
        y2_tot += cur_atom.y
        z2_tot += cur_atom.z
        
    c2x = x2_tot/ac2
    c2y = y2_tot/ac2
    c2z = z2_tot/ac2
    
    return ((c2x-c1x)**2+(c2y-c1y)**2+(c2z-c1z)**2)**(0.5)

#Function returns the Euclidean distance between the centroids formed
#by the side chain atoms of two residues.
def sidechain_centroid_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:    
    atom1is = np.array(list(res1.heavyatoms.keys()))
    atom2is = np.array(list(res2.heavyatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one atom")
    
    ac1 = res1.atomcount
    ac2 = res2.atomcount

    x1_tot = y1_tot = z1_tot = x2_tot = y2_tot = z2_tot = 0
        
    for i in atom1is[atom1is >= SIDE_CHAIN_START]:
        cur_atom = res1[int(i)]
        x1_tot += cur_atom.x
        y1_tot += cur_atom.y
        z1_tot += cur_atom.z
        
    c1x = x1_tot/(ac1-SIDE_CHAIN_START)
    c1y = y1_tot/(ac1-SIDE_CHAIN_START)
    c1z = z1_tot/(ac1-SIDE_CHAIN_START)


    for i in atom2is[atom2is >= SIDE_CHAIN_START]:
        cur_atom = res2[int(i)]
        x2_tot += cur_atom.x
        y2_tot += cur_atom.y
        z2_tot += cur_atom.z
        
    c2x = x2_tot/(ac2-SIDE_CHAIN_START)
    c2y = y2_tot/(ac2-SIDE_CHAIN_START)
    c2z = z2_tot/(ac2-SIDE_CHAIN_START)
    
    return ((c2x-c1x)**2+(c2y-c1y)**2+(c2z-c1z)**2)**(0.5)

#Function takes in two residues and outputs their minimum N-O distance
#if those residues are charged.
def NO_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    condition1 = (res1.name not in POSITIVE and res2.name not in NEGATIVE)
    condition2 = (res1.name not in NEGATIVE and res2.name not in POSITIVE)
    if  condition1 and condition2:
        raise ValueError("You must input a postive and negative residue")
    
    else:
        atom1is = np.array(list(res1.heavyatoms.keys()))
        atom2is = np.array(list(res2.heavyatoms.keys()))
        if np.size(atom1is) == 0 or np.size(atom2is) == 0:
            raise ValueError("Input residues with at least one atom")
        
        distances = []
        for atom1 in atom1is:
            cur_atom1 = res1[int(atom1)]
            if cur_atom1.sb_able == "+":
                for atom2 in atom2is:
                    cur_atom2 = res2[int(atom2)]
                    if cur_atom2.sb_able == "-":
                        distances.append(atom_dist(cur_atom1, cur_atom2))
            elif cur_atom1.sb_able =="-":
                for atom2 in atom2is:
                    cur_atom2 = res2[int(atom2)]
                    if cur_atom2.sb_able == "+":
                        distances.append(atom_dist(cur_atom1, cur_atom2))

        if distances == []:
            raise ValueError("{} {} and/or {} {} have missing charged atoms.".format(res1.name, res1.index, res2.name, res2.index))
        return min(distances)
    
#Function takes in one complex, a chainname, a resid, a residue distance, and a cutoff and outputs
#all of those residues that are within the cutoff distance.
def within_distance(complex: COMPLEX, chainname: str, resid: int, func, cutoff: float, uncharged = True) -> np.array:
    distances = []
    info = []
    
    avail_funcs = [minimum_dist, sidechain_dist, centroid_dist, sidechain_centroid_dist, alpha_dist]
    if func not in avail_funcs:
        raise Exception("Use one of minimum_dist, alpha_dist, sidechain_dist, centroid_dist, or sidechain_centroid_dist functions")
    
    query_res = complex[chainname][resid]
    
    for c in complex.chainnames:
        cur_chain = complex[c]
        for r in cur_chain.residues.keys():
            if c == chainname and r == resid:
                continue
            cur_res = cur_chain[r]

            dist = func(cur_res, query_res)

            if dist <= cutoff:
                if uncharged:
                    info.append([c,cur_res.name, cur_res.index])
                    distances.append(dist)
                elif cur_res.name in CHARGED:
                    info.append([c,cur_res.name, cur_res.index])
                    distances.append(dist)

    info = np.array(info)
    distances = np.array(distances)

    return info, distances
    
"""Complex distances beyond this point"""

#Function takes in a complex, a chainname, a cutoff, and a residue distance function as inputs
#and outputs all of the residue distances from the input chain to all of the other chains
#in the complex.
def all_dist(complex:COMPLEX, chainname1: str, cutoff: float, func) -> tuple:
    chain1 = complex[chainname1]
    heuristic = []
    distances = []
    info = []
    avail_funcs = [minimum_dist, sidechain_dist, centroid_dist, sidechain_centroid_dist, alpha_dist]
    if func not in avail_funcs:
        raise Exception("Use one of minimum_dist, alpha_dist, sidechain_dist, centroid_dist, or sidechain_centroid_dist functions")
    for chainame2 in complex.chainnames:
        if chainame2 != chainname1:
            chain2 = complex[chainame2]
            for res1i in chain1.residues.keys():
                cur_res1 = chain1[res1i]
                if cur_res1.name not in AA:
                    continue
                for res2i in chain2.residues.keys():
                    cur_res2 = chain2[res2i]
                    if cur_res2.name not in AA:
                        continue
                    cur_dist = centroid_dist(cur_res1, cur_res2)
                    if cur_dist < cutoff + 10:
                        heuristic.append([cur_res1.name, cur_res1.index, cur_res2.name, cur_res2.index, cur_dist, chainame2])

    for pairs in heuristic:
        chain2 = complex[pairs[5]]
        cur_res1 = chain1[pairs[1]]
        cur_res2 = chain2[pairs[3]]
        try:
            cur_dist = func(cur_res1, cur_res2)
        except:
            continue
        if cur_dist < cutoff:
            info.append([chainname1,cur_res1.name, cur_res1.index, pairs[5], cur_res2.name, cur_res2.index])
            
            distances.append(cur_dist)

    info = np.array(info)
    distances = np.array(distances)

    return info, distances

#Function returns all salt bridges for a complex, within a cutoff
def all_salt_bridges(complex: COMPLEX, sb_cutoff:float = 3.2) -> tuple:
    
    distances = []
    info = []

    for chain1 in complex.chainnames:
        cur_chain1 = complex[chain1]
        for res1 in cur_chain1.residues.keys():
            cur_res1 = cur_chain1[res1]
            if cur_res1.name not in AA:
                continue
            if cur_res1.charge == "+":
                    for chain2 in complex.chainnames:
                        cur_chain2 = complex[chain2]
                        for res2 in cur_chain2.residues.keys():
                            cur_res2 = cur_chain2[res2]
                            if cur_res2.name not in AA:
                                continue
                            if cur_res2.charge == "-":
                                try:
                                    dist = NO_dist(cur_res1, cur_res2)
                                    if dist <= sb_cutoff:
                                        info.append([chain1, cur_res1.name, cur_res1.index, chain2, cur_res2.name, cur_res2.index])
                                        distances.append(dist)
                                except Exception as e:
                                    continue


    distances = np.array(distances)
    info = np.array(info)

    return info, distances

#Function takes in a complex, a chain name, a distance cutoff, a salt-bridge cutoff
#and outputs all those uncharged residues which are in the vicinity of a charged residue.
#Optionally, the user can specify that the charged residues must not be involved in potential
#salt bridges. 
def find_putative_sb(complex: COMPLEX, chainname: str, cutoff: float, sb_cutoff: float, func, mind_sb=True) -> tuple:
    info = []
    distances = []

    all_info, all_distances = all_dist(complex, chainname, cutoff, func)

    if mind_sb:
        sb_info, _ = all_salt_bridges(complex, sb_cutoff)

        for cur_i, cur_d in zip(all_info, all_distances):
            if cur_i[1] in NEUTRAL:
                stop_iterate = False
                query_chain = cur_i[3]
                query_resname = cur_i[4]
                query_resid = cur_i[5]

                if query_resname in POSITIVE:
                    for sb_i in sb_info:
                        sb_chain = sb_i[0]
                        sb_resid = sb_i[2]

                        if (sb_chain == query_chain and sb_resid == query_resid):
                            stop_iterate = True
                            break

                elif query_resname in NEGATIVE:
                    for sb_i in sb_info:
                        sb_chain = sb_i[3]
                        sb_resid = sb_i[5]

                        if (sb_chain == query_chain and sb_resid == query_resid):
                            stop_iterate = True
                            break

                if stop_iterate:
                    continue
                if query_resname in NEUTRAL:
                    continue
                else:
                    info.append(cur_i)
                    distances.append(cur_d)

    else:
        for cur_i, cur_d in zip(all_info, all_distances):
            if cur_i[1] in NEUTRAL and cur_i[4] in CHARGED:
                info.append(cur_i)
                distances.append(cur_d)

    info = np.array(info)
    distances = np.array(distances)

    return info, distances

"""Data structures beyond this point"""

res_func_dict = {"minimum_distance": minimum_dist, "sidechain_distance": sidechain_dist,
                "alpha_carbon_distance": alpha_dist, "centroid_distance": centroid_dist,
                "sidechain_centroid_distance": sidechain_centroid_dist}