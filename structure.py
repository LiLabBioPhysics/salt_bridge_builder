import subprocess
import os

class ATOM:
    def __init__(self, x:float, y:float, z:float, type: str = None):
        self.x = x
        self.y = y
        self.z = z
        self.type = type

class RESIDUE:
    CTERATOMS = {"C"}
    NTERATOMS = {"HN", "N", "HT3"}
    ALPHAATOM = "CA"

    def __init__(self):
        self.name = None
        self.index = None
        self.charge = None
        self.atoms = dict()
        self.atomcount = 0
        self.atommap = dict()
        self.sideatoms = dict()
        self.alpha = None

    def add_atom(self, atom: ATOM, aposition: int):
        self.atoms[aposition] = atom
        self.atomcount += 1
        self.atommap[self.atomcount] = aposition

    def __getitem__(self, key):
        return self.atoms[self.atommap[key]]
    
    def get_alpha(self):
        for i in range(self.atomcount, 0, -1):
            if self.atoms[self.atommap[i]].type == self.ALPHAATOM:
                self.alpha = i
                break

    def get_centroid(self):
        aindex = list(range(self.atomcount))
        n = len(aindex)
        x_tot = y_tot = z_tot = 0
        
        for i in aindex:
            cur_atom = self.atoms[self.atommap[i+1]]
            x_tot += cur_atom.x
            y_tot += cur_atom.y
            z_tot += cur_atom.z
        
        self.centroidx = x_tot/n
        self.centroidy = y_tot/n
        self.centroidz = z_tot/n

    def get_side(self):
        for i in range(self.atomcount, 0, -1):
            if self.atoms[self.atommap[i]].type in self.CTERATOMS:
                continue
            if self.atoms[self.atommap[i]].type in self.NTERATOMS:
                break
            else:
                self.sideatoms[i] = self.atoms[self.atommap[i]]

    def get_side_centroid(self):
        self.get_side()

        aindex = sorted(list(self.sideatoms.keys()))
        n = len(aindex)
        x_tot = y_tot = z_tot = 0
        
        for i in aindex:
            cur_atom = self.atoms[self.atommap[i]]
            x_tot += cur_atom.x
            y_tot += cur_atom.y
            z_tot += cur_atom.z
        
        self.sidecentroidx = x_tot/n
        self.sidecentroidy = y_tot/n
        self.sidecentroidz = z_tot/n

class CHAIN:

    positive = {"HIS", "LYS", "ARG"}
    negative = {"ASP", "GLU"}
    charged = positive | negative
    neutral = {"GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "MET", "PHE",
    "TRP", "PRO", "TYR", "GLN", "ASN", "CYS", "THR", "SER"}

    def charge(self, name:str):
        if name in self.positive:
            return "+"
        if name in self.negative:
            return "-"
        if name in self.neutral:
            return "0"
        else:
            return None

    def __init__(self):
        self.residues = dict()
    
    def add_residue(self, residue: RESIDUE, rname : str, rposition: int):
        self.residues[rposition] = residue
        self.residues[rposition].index = rposition
        self.residues[rposition].name = rname
        self.residues[rposition].charge = self.charge(rname)
    
    def __getitem__(self, key):
        return self.residues[key]

"""
A complex is defined by a single pdb file containing one or more protein chains.
"""       

class COMPLEX:

    def add_chain(self, chain: CHAIN, chainname: str):
        self.chains[chainname] = chain

    def __getitem__(self, key):
        return self.chains[key]

    #Creating a COMPLEX directly from a single pdb file or 
    # from the data passed by the TRAJECTORY class 
    def __init__(self, pdbfile=None, frame_data=None):

        self.chains = dict()
        self.chainnames = self.chains.keys()
        
        if pdbfile is not None:
            frame_data = []
            # self.pdb_clean = ""
            with open(pdbfile, "r") as f:
                data = f.readlines()
            for line in data:
                if line[0:4] == "ATOM":
                    frame_data.append(line)

        else:
            tmp_data = []
            for line in frame_data:
                 if line[0:4] == "ATOM":
                    tmp_data.append(line)
            frame_data = tmp_data

        for line in frame_data:
            atomnum = int(line[6:11].strip())
            atomtype = line[12:16].strip()
            resname = line[17:20]
            chainname = line[21]
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            cur_atom = ATOM(x, y, z, atomtype)

            no_chain = (not bool(self.chains)) or (chainname not in self.chains.keys())

            if no_chain:
                self.add_chain(CHAIN(), chainname)
            
            cur_chain = self.chains[chainname]

            no_res = (not bool (cur_chain.residues)) or (resnum not in cur_chain.residues.keys())

            if no_res:
                cur_res = RESIDUE()
                cur_chain.add_residue(cur_res, resname, resnum)

            self.chains[chainname].residues[resnum].add_atom(cur_atom, atomnum)

    # def split_chain(self, chain: str, write=False):
    #     CHAIN_COLUMN = 21
    #     self.remove_chain = ""
    #     self.keep_chains = "" 

    #     for line in self.pdb_clean.splitlines():
    #         if line[CHAIN_COLUMN] == chain:
    #             self.remove_chain = self.remove_chain+line+"\n"
    #         else:
    #             self.keep_chains=self.keep_chains+line+"\n"
        
    #     if self.remove_chain == "":
    #         raise ValueError("Your chain should exist in the pdb file.")
        
    #     if write:
    #         with open("seperated_chain.pdb", "w") as f:
    #             f.write(self.remove_chain)

    #         with open("unseperated_chains.pdb", "w") as f:
    #             f.write(self.keep_chains)

    #     if not write:
    #         os.remove("seperated_chain.pdb")
    #         os.remove("seperated_chain.pqr")
    #         os.remove("unseperated_chain.pdb")
    #         os.remove("unseperated_chain.pqr")

    # def run_DelPhiForce(self, chain, output):
    #     self.make_pqr_files(chain)
    #     subprocess.run(["bash", "DelPhiForce/bin/DelPhiForce.sh", "-1", "unseperated_chains.pqr", "-2", "seperated_chain.pqr", "-e", "./DelPhiForce/bin/delphicpp", "-o", output])

"""
A trajectory is defined by a multi-pdb file where each frame represents a protein  
complex at a specific timestep. Note: trajectories are zero-indexed.
"""           
class TRAJECTORY:
    #read_start_points() returns the first line number in a multipdb file for each frame.
    #Doing so ensures that not the whole multipdb file is read into memory.
    def read_start_points(self, multipdbfile: str) -> None:
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

    def read_frame(self, start) -> None:
        #read_frame() parses the multipdb file at the intervals found by start_points() and returns parsed data
        #in list form.
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

    def __init__(self, multipdbfile:str) -> None:
        self.filename = multipdbfile
        self.start_points = self.read_start_points(multipdbfile)
        self.frame = 0
        self.add_complex(self.read_frame(self.frame))

    def add_complex(self, data):
        self.complex = COMPLEX(frame_data = data)
        
    def __len__(self):
        return len(self.start_points)

    def __getitem__(self, key):
        if key == self.frame:
            return self.complex
        
        elif (key >= len(self)):
            raise IndexError
        
        else:
            self.frame = key
            self.add_complex(self.read_frame(key))
            return self.complex