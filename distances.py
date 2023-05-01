from structure import *
import numpy as np

def atom_dist(atom1:ATOM, atom2:ATOM) -> float:
    x1= atom1.x
    y1= atom1.y
    z1= atom1.z

    x2= atom2.x
    y2= atom2.y
    z2= atom2.z
    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(0.5)


def pairwise_dist(res1: RESIDUE, res2: RESIDUE) -> np.ndarray:
    a1index = list(range(res1.atomcount))
    a2index = list(range(res2.atomcount))


    arr = np.zeros((a1index[-1]+1, a2index[-1]+1))
 
    for i in a1index:
        for j in a2index:
            a_dist = atom_dist(res1[i+1], res2[j+1])
            arr[i,j] = a_dist

    return arr

def minimum_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    pairwise_arr = pairwise_dist(res1, res2)
    return np.min(pairwise_arr)

def alpha_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    res1.get_alpha()
    res2.get_alpha()
    
    dist = atom_dist(res1[res1.alpha], res2[res2.alpha])

    return dist

def sidechain_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    res1.get_side()
    res2.get_side()

    s1index = sorted(list(res1.sideatoms.keys()))
    s2index = sorted(list(res2.sideatoms.keys()))

    arr = np.zeros((len(s1index), len(s2index)))
 
    for i in range(len(s1index)):
        for j in range(len(s2index)):
            a_dist = atom_dist(res1[s1index[i]], res2[s2index[j]])
            arr[i,j] = a_dist

    return np.min(arr)


def centroid_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    res1.get_centroid()
    res2.get_centroid()

    x1= res1.centroidx
    y1= res1.centroidy
    z1= res1.centroidz

    x2= res2.centroidx
    y2= res2.centroidy
    z2= res2.centroidz
    
    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(0.5)

def sidechain_centroid_dist(res1: RESIDUE, res2: RESIDUE) -> np.float64:
    res1.get_side_centroid()
    res2.get_side_centroid()

    x1= res1.sidecentroidx
    y1= res1.sidecentroidy
    z1= res1.sidecentroidz

    x2= res2.sidecentroidx
    y2= res2.sidecentroidy
    z2= res2.sidecentroidz

    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(0.5)


def within_dist(complex:COMPLEX, chainname1: str, cutoff: float, func) -> list:
    chain1 = complex[chainname1]
    heuristic = []
    distances = []
    avail_funcs = [minimum_dist, sidechain_dist, centroid_dist, sidechain_centroid_dist, alpha_dist]
    if func not in avail_funcs:
        raise Exception("Use one of minimum_dist, alpha_dist, sidechain_dist, centroid_dist, or sidechain_centroid_dist functions")

    for chainame2 in complex.chainnames:
        if chainame2 != chainname1:
            chain2 = complex[chainame2]
            for res1i in chain1.residues.keys():
                cur_res1 = chain1[res1i]
                for res2i in chain2.residues.keys():
                    cur_res2 = chain2[res2i]
                    cur_dist = centroid_dist(cur_res1, cur_res2)
                    if cur_dist < cutoff + 10:
                        heuristic.append([cur_res1.name, cur_res1.index, cur_res2.name, cur_res2.index, cur_dist, chainame2])


    for pairs in heuristic:
        chain2 = complex[pairs[5]]
        cur_res1 = chain1[pairs[1]]
        cur_res2 = chain2[pairs[3]]
        cur_dist = func(cur_res1, cur_res2)
        if cur_dist < cutoff:
            distances.append([chainname1,cur_res1.name, cur_res1.index, pairs[5], cur_res2.name, cur_res2.index, cur_dist])
                
    return distances


