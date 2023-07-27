"""
The distances module contains definitions for all residue distances.
It also contains various useful utility functions.
"""

from structure import *
from constants import *
import numpy as np


def atom_dist(atom1:ATOM, atom2:ATOM) -> float:
    """
    Parameters
    ----------

    atmo1 : ATOM
        The first atom.

    
    atmo2 : ATOM
        The second atom.
    
    Returns
    ------
    float
        Euclidean distance between two atoms (points).
    """
    diff = atom2.coordinates - atom1.coordinates
    return (np.linalg.norm(diff))

def get_centroid(res: RESIDUE) -> np.ndarray:
    """
    Parameters
    ----------

    res : RESIDUE
        The residue for which the centroid is calculated.
    
    Raises
    ------
    ValueError
        If the residue does not have at least one atom.

    Returns
    ------
    np.ndarray
        The residue centroid.
    """
    atomis = np.array(list(res.heavyatoms.keys()))
    if np.size(atomis) == 0:
        raise ValueError("Input a residue with at least one atom")
    
    centroid = np.array([0,0,0])
    
    for a in res.heavyatoms.values():
        centroid = centroid + a.coordinates

    centroid = centroid / res.atomcount

    return(centroid)


def get_sidechain_centroid(res: RESIDUE) -> np.ndarray:
    """
    Parameters
    ----------

    res : RESIDUE
        The residue for which the side chain centroid is calculated.
    
    Raises
    ------
    ValueError
        If the residue does not have at least one side chain atom.

    Returns
    ------
    np.ndarray
        The residue side chain centroid.
    """
    atomis = np.array(list(res.sideatoms.keys()))
    if np.size(atomis) == 0:
        raise ValueError("Input a residue with at least one side chain atom")
    
    sc_centroid = np.array([0,0,0])
    
    for a in res.sideatoms.values():
        sc_centroid = sc_centroid + a.coordinates

    sc_centroid = sc_centroid / res.sideatomcount

    return(sc_centroid)


def pairwise_dist(res1: RESIDUE, res2: RESIDUE) -> set:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.
    
    Raises
    ------
    ValueError
        If the residues do not have at least one atom.

    Returns
    ------
    set
        Pairwise atom distances are calculated between all atoms. A set of these floats is returned.
    """
    atom1is = np.array(list(res1.heavyatoms.keys()))
    atom2is = np.array(list(res2.heavyatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one atom")
    
    pairwise_set = set()

    for i in atom1is:
        for j in atom2is:
            a_dist = float(atom_dist(res1[int(i)], res2[int(j)]))
            pairwise_set.add(a_dist)

    return pairwise_set

def pairwise_sidechain_dist(res1: RESIDUE, res2: RESIDUE) -> set:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.
    
    Raises
    ------
    ValueError
        If the residues do not have at least one atom.

    Returns
    ------
    set
        Pairwise atom distances for side chain atoms are calculated. A set of these floats is returned.
    """
    atom1is = np.array(list(res1.sideatoms.keys()))
    atom2is = np.array(list(res2.sideatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one side chain atom")

    side_set = set()

    for i in atom1is:
        for j in atom2is:
            a_dist = float(atom_dist(res1[int(i)], res2[int(j)]))
            side_set.add(a_dist)
    return side_set

def minimum_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.

    Returns
    ------
    float
        Pairwise atom distances for all atoms are calculated. The minimum value of this set is returned.
    """
    pairwise_set = pairwise_dist(res1, res2)
    return min(pairwise_set)

def sidechain_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.

    Returns
    ------
    float
        Pairwise atom distances for side chain atoms are calculated. The minimum value of this set is returned.
    """
    side_set = pairwise_sidechain_dist(res1, res2) 
    return min(side_set)

def alpha_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.

    Returns
    ------
    float
        The alpha carbon distance between two residues.
    """
    dist = atom_dist(res1[ALPHA_CARBON_POS], res2[ALPHA_CARBON_POS])
    return float(dist)

def centroid_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.

    Returns
    ------
    float
        Function returns the Eulidean distance between the centroids of two residues.
    """
    atom1is = np.array(list(res1.heavyatoms.keys()))
    atom2is = np.array(list(res2.heavyatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one atom")
    
    diff = get_centroid(res2) - get_centroid(res1)

    return float(np.linalg.norm(diff))

def sidechain_centroid_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.

    Returns
    ------
    float
        Function returns the Eulidean distance between the side chain centroids of two residues.
    """    
    atom1is = np.array(list(res1.heavyatoms.keys()))
    atom2is = np.array(list(res2.heavyatoms.keys()))

    if np.size(atom1is) == 0 or np.size(atom2is) == 0:
        raise ValueError("Input residues with at least one atom")
    
    diff = get_sidechain_centroid(res1) - get_sidechain_centroid(res2)

    return np.linalg.norm(diff)

def NO_dist(res1: RESIDUE, res2: RESIDUE) -> float:
    """
    Parameters
    ----------
    res1 : RESIDUE
        The first residue.

    res1 : RESIDUE
        The second residue.
    
    Raises
    ------
    ValueError
        If the residues are not charged or if there is not at least one charged atom in both residues.

    Returns
    ------
    float
        The minimum value of a set of pairwise differences between charged side chain atoms.
    """
    condition1 = (res1.name not in POSITIVE and res2.name not in NEGATIVE)
    condition2 = (res1.name not in NEGATIVE and res2.name not in POSITIVE)
    if  condition1 and condition2:
        raise ValueError("You must input a postive and negative residue")
    
    else:
        atom1is = np.array(list(res1.heavyatoms.keys()))
        atom2is = np.array(list(res2.heavyatoms.keys()))
        if np.size(atom1is) == 0 or np.size(atom2is) == 0:
            raise ValueError("Input residues with at least one atom")
        
        distances = {}
        for atom1 in atom1is:
            cur_atom1 = res1[int(atom1)]
            if cur_atom1.sb_able == "+":
                for atom2 in atom2is:
                    cur_atom2 = res2[int(atom2)]
                    if cur_atom2.sb_able == "-":
                        distances.add(atom_dist(cur_atom1, cur_atom2))
            elif cur_atom1.sb_able =="-":
                for atom2 in atom2is:
                    cur_atom2 = res2[int(atom2)]
                    if cur_atom2.sb_able == "+":
                        distances.add(atom_dist(cur_atom1, cur_atom2))

        if distances == {}:
            raise ValueError("{} {} and/or {} {} have missing charged atoms.".format(res1.name, res1.index, res2.name, res2.index))
        return min(distances)
    
def within_distance(complex: COMPLEX, chainname: str, resid: int, func: function, cutoff: float, uncharged = True) -> tuple:
    """
    Parameters
    ----------
    complex : COMPLEX
        The complex to iterate over.

    chainname : str
        The query chain.

    resid : res
        The query residue index.

    func : function
        The function used to calculate residue distances.

    cutoff : float
        The residue distance cutoff.

    uncharged : bool
        Weather to include uncharged residues in the output.
    
    Raises
    ------
    Exception
        If the chosen distance is not defined.

    Returns
    ------
    tuple
        A tuple of values where where the two objects in the tuple are information and distance np.ndarrays.
        The objects represent all residues which are within the cutoff distance to the query residue.
        The information array contains information about the chain name, residue name, and reside index.
        The distance array contains the residue distances/s.
    """

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

def all_dist(complex:COMPLEX, chainname1: str, cutoff: float, func) -> tuple:
    """
    Parameters
    ----------
    complex : COMPLEX
        The complex to iterate over.

    chainname1 : str
        The query chain.

    cutoff : float
        The residue distance cutoff.

    func : function
        The function used to calculate residue distances.
    
    Raises
    ------
    Exception
        If the chosen distance is not defined.

    Returns
    ------
    tuple
        A tuple of values where where the two objects in the tuple are information and distance np.ndarrays.
        The objects represent all pairiwise residue distances from all residues in the query chain to all residues in all other chains.
        All residues are within the residue distance cutoff. An alpha carbon distance of 16 Å is used to shortlist residue distances.
        The information array contains information about the chain names, residue names, and reside indexes. 
        The first three columns represent the query residues. The last three columns represent all other chains.
        The distance array contains the residue distances/s.
    """
    chain1 = complex[chainname1]
    interface = []
    info = []
    distances = []
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
                    cur_dist = alpha_dist(cur_res1, cur_res2)
                    if cur_dist <= INTERFACE_CUTOFF*2:
                        if func == alpha_dist:
                            if cur_dist <= cutoff:
                                info.append([chainname1, cur_res1.name, cur_res1.index, chainame2, cur_res2.name, cur_res2.index])
                                distances.append(cur_dist)
                        
                        else:
                            interface.append([cur_res1.name, cur_res1.index, cur_res2.name, cur_res2.index, cur_dist, chainame2])
    
    if func == alpha_dist:
        info = np.array(info)
        distances = np.array(distances)

        return info, distances

    for pairs in interface:
        chain2 = complex[pairs[5]]
        cur_res1 = chain1[pairs[1]]
        cur_res2 = chain2[pairs[3]]
        try:
            cur_dist = func(cur_res1, cur_res2)
        except:
            continue
        if cur_dist <= cutoff:
            info.append([chainname1,cur_res1.name, cur_res1.index, pairs[5], cur_res2.name, cur_res2.index])
            
            distances.append(cur_dist)

    info = np.array(info)
    distances = np.array(distances)

    return info, distances

def all_salt_bridges(complex: COMPLEX, sb_cutoff:float = 3.2) -> tuple:
    """
    Parameters
    ----------
    complex : COMPLEX
        The complex to iterate over.

    sb_cutoff : float
        The salt bridge distance cutoff.

    Returns
    ------
    tuple
        A tuple of values where where the two objects in the tuple are information and distance np.ndarrays.
        The objects represent all salt bridges in the complex (including intraprotein salt bridges).
        The information array contains information about the chain names, residue names, and reside indexes. 
        The first three columns represent the query residues. The last three columns represent all other chains.
        The distance array contains the residue distances/s.
    """
    
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
    """
    Parameters
    ----------
    complex : COMPLEX
        The complex to iterate over.

    chainname1 : str
        The query chain.

    cutoff : float
        The residue distance cutoff.

    sb_cutoff : float
        The salt bridge distance cutoff.

    func : function
        The function used to calculate residue distances.

    mind_sb : bool
        Whether to disregard those putative salt bridges where the wild-type charged residue is already involved in a salt bridge. 
    
    Raises
    ------
    Exception
        If the chosen distance is not defined.

    Returns
    ------
    tuple
        A tuple of values where where the two objects in the tuple are information and distance np.ndarrays.
        The objects represent putative salt bridges based on the residue distance cutoff.
        The information array contains information about the chain names, residue names, and reside indexes. 
        The first three columns represent the query residues. The last three columns represent all other chains.
        The distance array contains the residue distances/s.
    """
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

res_func_dict = {"minimum_distance": minimum_dist, "sidechain_distance": sidechain_dist,
                "alpha_carbon_distance": alpha_dist, "centroid_distance": centroid_dist,
                "sidechain_centroid_distance": sidechain_centroid_dist}
