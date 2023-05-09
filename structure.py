from constants import *

"""
An atom is the basic building block for all Bridge Builder
objects. An atom has x/y/z coordinates, a type (based on its parent residue),
a charge, and a pdb serial number. If a side chain atom has the possibility to be ionized, then it will be
considered either positve (nitrogen) or negative (oxygen), regardless of the ambient
conditions. 
"""

class ATOM:    
    def __init__(self, x:float, y:float, z:float, pdbpos: int, type: str):
        self.x = x
        self.y = y
        self.z = z
        self.type = type
        self.pdbpos = pdbpos
        
        if self.type in NEGATOM:
            self.sb_able = "-"
        elif self.type in POSATOM:
            self.sb_able = "+"
        else:
            self.sb_able = "0"

"""A residue is a collection of atoms and represents one of the cannonical
amino acids. A residue can be postive, negative, or charged
according to the definitions in constants.py."""

class RESIDUE:
    def __init__(self, rname : str, rposition: int):
        self.name = rname
        self.index = rposition
        self.heavyatoms = dict()
        self.atomcount = 0
        self.otheratoms = dict()
        if self.name in POSITIVE:
            self.charge = "+"
        elif (self.name in NEGATIVE):
            self.charge = "-"
        elif (self.name in NEUTRAL):
            self.charge = "0"
        else:
            self.charge = None
    
    #The internal numbering of the PDB file is stored in the atom object but not used 
    #to index atoms. Instead, heavy atom numbering is done according to the atom positions outlined
    #in constants.py. These tend to flow from the atom backbone away to the 
    #side chain.
    def add_atom(self, atom: ATOM):
        if atom.type in HEAVY_ATOM_NAMES:
            try:
                aposition = HEAVY_ATOM_POS[self.name].index(atom.type)
               
            except:
                aposition = HEAVY_ATOM_POS_ALT[self.name][atom.type]
            finally:
                self.heavyatoms[aposition] = atom
                self.atomcount += 1
        #Functionality will be added later for hydrogen numbering.
        else:
            self.otheratoms[atom.pdbpos] = atom

    def __getitem__(self, key: int) -> ATOM:
        if not isinstance(key, int):
            raise TypeError("You must index atoms by an integer")
        return self.heavyatoms[key]
    
"""Chains are peptides or protein monomers. 
Residues can be indexed by their pdb file residue indexes."""

class CHAIN:
    def __init__(self, chainame: str):
        self.residues = dict()
        self.name = chainame
    
    def add_residue(self, residue: RESIDUE):
        self.residues[residue.index] = residue

    def __getitem__(self, key: int) -> RESIDUE:
        if not isinstance(key, int):
            raise TypeError("You must index residues by an integer")
        return self.residues[key]

"""
A complex is defined by a single pdb file containing one or more protein chains.
A complex can be indexed by the chain name.
"""  
class COMPLEX:
    def add_chain(self, chain: CHAIN):
        self.chains[chain.name] = chain
        self.chainnames.append(chain.name)

    def __getitem__(self, key: str) -> CHAIN:
        if not isinstance(key, str):
            raise TypeError("You must index chains by the chain letter")
        return self.chains[key]

    #Here we create a COMPLEX directly from a single pdb file or 
    #from the data passed by the TRAJECTORY class. 
    def __init__(self, pdbfile=None, frame_data=None):
        self.chains = dict()
        self.chainnames = []
        
        #Creating from the single pdb file.
        if pdbfile is not None:
            frame_data = []
            with open(pdbfile, "r") as f:
                data = f.readlines()
            for line in data:
                if line[0:4] == "ATOM":
                    frame_data.append(line)
        
        #Creating from the multi-frame pdb file.
        else:
            tmp_data = []
            for line in frame_data:
                 if line[0:4] == "ATOM":
                    tmp_data.append(line)
            frame_data = tmp_data

        #Column defintions from PDB files are used to obtain atom information
        #and create the COMPLEX.
        for line in frame_data:
            atomnum = int(line[6:11].strip())
            atomtype = line[12:16].strip()
            resname = line[17:20]
            chainname = line[21]
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            cur_atom = ATOM(x, y, z, atomnum, atomtype)

            no_chain = (not bool(self.chains)) or (chainname not in self.chainnames)

            if no_chain:
                self.add_chain(CHAIN(chainname))
            
            cur_chain = self.chains[chainname]

            no_res = (not bool (cur_chain.residues)) or (resnum not in cur_chain.residues.keys())

            if no_res:
                cur_res = RESIDUE(resname, resnum)
                cur_chain.add_residue(cur_res)

            self.chains[chainname].residues[resnum].add_atom(cur_atom)

"""
A trajectory is defined by a multi-frame pdb file where each frame represents a protein  
complex at a specific timestep. Note: trajectories are zero-indexed.
"""
#Remember to add auxiliary pdb file cleaning functions.           
class TRAJECTORY:
    #read_start_points() returns the first line number in a multi-frame pdb file for each frame.
    #Doing so ensures that not the whole multipdb file is read into memory when you create a complex.
    def read_start_points(self, multipdbfile: str) -> list:
        count = 1
        start_points = [1]
        with open(multipdbfile, "r") as f:
            curr_line= f.readline()
            while (curr_line):
                curr_line = f.readline()
                if curr_line[0:3] == "END":
                    start_points.append(count+2)
                count += 1

        start_points.pop(-1)
        return start_points

    #read_frame() parses the multipdb file at the intervals found by start_points() and returns parsed data
    #in list form.
    def read_frame(self, start) -> list:
        data = []
        
        if len(self.start_points) == 1:
            with open(self.filename, "r") as f:
                data = f.readlines()
            data.pop(-1)

        else:
            data = []
            dt = self.start_points[1] - self.start_points[0] 
            specified_lines = range(start*dt, start*dt+1 +dt)
            with open(self.filename, "r") as f:
                for pos, l_num in enumerate(f):
                    if pos in specified_lines:
                        data.append(l_num)
            data.pop(-1)

        return data

    #Remember to add the ability to read more frames into memory.
    def __init__(self, multipdbfile:str):
        self.filename = multipdbfile
        self.start_points = self.read_start_points(multipdbfile)
        self.frame = 0
        self.add_complex(self.read_frame(self.frame))

    def add_complex(self, data: list):
        self.complex = COMPLEX(frame_data = data)
        
    def __len__(self):
        return len(self.start_points)

    def __getitem__(self, key: int) -> COMPLEX:
        if not isinstance(key, int):
            raise TypeError("You must index frames by an integer")
        if key == self.frame:
            return self.complex
        
        elif (key >= len(self)):
            raise IndexError
        
        else:
            self.frame = key
            self.add_complex(self.read_frame(key))
            return self.complex