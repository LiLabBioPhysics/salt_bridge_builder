class ATOM:
    def __init__(self, x:float, y:float, z:float, type: str = None):
        self.x = x
        self.y = y
        self.z = z
        self.type = type

class RESIDUE:
    CTERATOMS = {"C"}
    NTERATOMS = {"HN", "N", "HT3"}

    def __init__(self):
        self.name = None
        self.index = None
        self.charge = None
        self.atoms = dict()
        self.atomcount = 0
        self.atommap = dict()
        self.sideatoms = dict()

    def add_atom(self, atom: ATOM, aposition: int):
        self.atoms[aposition] = atom
        self.atomcount += 1
        self.atommap[self.atomcount] = aposition

    def __getitem__(self, key):
        return self.atoms[self.atommap[key]]

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


class COMPLEX:

    def add_chain(self, chain: CHAIN, chainname: str):
        self.chains[chainname] = chain

    def __getitem__(self, key):
        return self.chains[key]

    def __init__(self, frame_data):

        self.chains = dict()
        self.chainnames = self.chains.keys()

        for line in frame_data:
            atomnum = int(line[0])
            atomtype = line[1]
            resname = line[2]
            chainname = line[3]
            resnum = int(line[4])
            x = float(line[5])
            y = float(line[6])
            z = float(line[7])

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

                
    
class TRAJECTORY:
    def read_start_points(self, multipdbfile: str):
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


    def read_frame(self, start = 0):
        data = []
        if len(self.start_points) == 1:
            with open(self.filename, "r") as f:
                data = f.readlines()
            data.pop(-1)

        else:
            data = []
            dt = self.start_points[1] - self.start_points[0] 
            specified_lines = range(start, dt)
            with open(self.filename, "r") as f:
                for pos, l_num in enumerate(f):
                    if pos in specified_lines:
                        data.append(l_num)
            data.pop(-1)

        atoms = [x.split()[1:9] for x in data if x[0:4]=="ATOM"]
        
        return atoms


    def __init__(self, multipdbfile:str, complexname = None) -> None:
        self.filename = multipdbfile
        self.start_points = self.read_start_points(multipdbfile)
        self.frame = 0
        self.add_complex(self.read_frame(self.frame))


    def add_complex(self, data):
        self.complex = COMPLEX(data)
        
    def __len__(self):
        return len(self.start_points)

    def __getitem__(self, key):
        if key == self.frame:
            return self.complex

        else:
            self.frame = key
            self.add_complex(self.read_frame(key))
            return self.complex
    


    