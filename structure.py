"""
The structure module is used to parse and load PDB files into memory.

The stucture hierarchy is atom -> residue -> chain -> complex -> trajectory.

Trajectories are indexed by frame number to return complexes.
Complexes are indexed by chain ID to return chains.
Chains are indexed by residue ID to return residues.
Residues are indexed by atom number (positions in the constants module) to 
return atoms.
Atoms contain coordinates used for distance calculations (numpy arrays).

PDB files are read as complexes (single-frame PDB) or trajectories 
(multi-frame PDB).
"""

from constants import *
import numpy as np
import warnings

class ATOM:
    """
    An atom is the basic building block for all of Bridge Builder
    objects.
    
    ...

    Attributes
    ----------

    coordinates : numpy.ndarray of shape (3,)
        An atom has x/y/z coordinates saved as a numpy array.

    type : str
        The type of atom within the residue (e.g. alpha carbon)

    pdbpos : int
        The number of the atom in the PDB file

    alt_location : False | str
        When multiple coordinates may exist for an atom, there
        can be alternate atom locations. These are denotated by
        single letter.

    sb_able : False | str
        Wheather an atom is considered positive, negative, or
        neutral is denoted by +, - , or 0, respectively.

    """

    def __init__(self, x:float, y:float, z:float, pdbpos: int, type: str, alt_location: str = False):
        """
        Parameters
        ----------

        x : float
            The x coordinate for an atom

        y : float
            The y coordinate for an atom

        z : float
            The z coordinate for an atom

        type : str
            The type of atom within the residue (e.g. alpha carbon)

        pdbpos : int
            The number of the atom in the PDB file

        alt_location : False | str
            When multiple coordinates may exist for an atom, there
            can me alternate atom locations. These are denotated by
            single letter.

        """
        self.coordinates = np.array([x, y, z])
        self.type = type
        self.pdbpos = pdbpos
        self.alt_location = alt_location
        if self.type in NEGATOM:
            self.sb_able = "-"
        elif self.type in POSATOM:
            self.sb_able = "+"
        else:
            self.sb_able = "0"

class RESIDUE:
    """A residue is a collection of atoms and represents one of the cannonical
    amino acids. A residue can be postive, negative, or charged
    according to the definitions in the constants module. A warning is displayed
    if a non-cannonical RESIDE object is created.

    ...
    
    Attributes
    ----------
    name : str
        The name of the residue.
    
    index : int
        The position of the residue from the PDB file

    heavyatoms: dict
        A dictionary where the keys are integers whose value is
        determined by the residue name and atom type (aposition). The positions
        are given in the constants module. The values are ATOM objects.

    atomcount = int
        The number of heavy non-alternate atoms in the residue. Simply,
        just the number of heavy atoms in the residue. Used for calculating
        the centroid.

    otheratoms = dict
        A dictionary where the keys are the PDB numbers of non-cannonical atoms
        (atoms not present in the 20 cannonical amino acids). The values are
        ATOM objects.

    unkownatoms = dict
        A dictionary with information about atoms from unkown residues. Sometimes
        residue identities cannot be established. The unkownatoms dictionary has keys
        whose values are the PDB numbers of atoms. The values are the atom objects.

    altatoms = dict
        A dictonary with information about alternate atoms.
        When a PDB file is read, if there are alternate atoms, only the first
        is kept as a heavy atom and as an alternate atom. 
        The key is the atom type and the value is the alternate location ATOM
        attribute.

    charge = str
        Wheather a residue is considered positive, negative, or
        neutral is dependent on the residue name
        and is denoted by +, - , or 0, respectively.

    Methods
    -------
    add_atom(atom)
        Used to populate atoms into residue.

    __getitem__(key)
        Used to index the atoms in a residue.

    """


    def __init__(self, rname : str, rposition: int):
        """
        Parameters
        ----------

        rname : str
            PDB name for the reside

        rname : int
            PDB index for the reside

        """
        self.name = rname
        self.index = rposition
        self.heavyatoms = dict()
        self.atomcount = 0
        self.otheratoms = dict()
        self.unkownatoms = dict()
        self.altatoms = dict()
        if self.name in POSITIVE:
            self.charge = "+"
        elif (self.name in NEGATIVE):
            self.charge = "-"
        elif (self.name in NEUTRAL):
            self.charge = "0"
        else:
            warnings.warn("There is at least one non-cannonical amino acid in the pdb file: {}".format(self.name))
            self.charge = None    

    def add_atom(self, atom: ATOM):
        """
        Parameters
        ----------

        atom : ATOM
            Atom object to be added

        """
        if isinstance(atom.alt_location, str):
            #Making sure the atom has no alternate locations.
            if atom.type in self.altatoms.keys():
                return None
            else:
                #Keeping the first alternate location.
                self.altatoms[atom.type] = atom.alt_location
        if self.name not in AA:
            #Making sure the residue is cannonical.
            self.unkownatoms[atom.pdbpos] = atom
        else:
            #When there are no special cases, we add assign a heavy atom position to the atom 
            #based on residue name and atom type.
            if atom.type in HEAVY_ATOM_NAMES:
                try:
                    aposition = HEAVY_ATOM_POS[self.name].index(atom.type)
                
                except:
                    aposition = HEAVY_ATOM_POS_ALT[self.name][atom.type]

                finally:
                    #Most atoms should be added here.
                    self.heavyatoms[aposition] = atom
                    self.atomcount += 1
            #Functionality will be added later for hydrogen numbering.
            else:
                self.otheratoms[atom.pdbpos] = atom

    def __getitem__(self, key: int) -> ATOM:
        """
        Parameters
        ----------

        key : int
            Used to index the atom. Overloaded for the [] operator.
        
        Raises
        ------
        TypeError
            If the atoms are not indexed by an integer.

        Returns
        ------
        ATOM
            The indexed atom.
        """
        if not isinstance(key, int):
            raise TypeError("You must index atoms by an integer")
        return self.heavyatoms[key]
    
class CHAIN:
    """
    Chains are peptides or protein monomers. 
    Residues can be indexed by their pdb file residue indexes.

    ...

    Attributes
    ----------

    residues : dict
        A dictionary where the keys are the residue indexes and values are the residues.
    name : str
        The chain name

    Methods
    -------
    add_residue(residue)
        Used to populate residues into the chain.

    __getitem__(key)
        Used to index the residues in a chain.

    """

    def __init__(self, chainame: str):
        """
        Parameters
        ----------

        chainname : str
            PDB name for the chain

        """
        self.residues = dict()
        self.name = chainame
    
    def add_residue(self, residue: RESIDUE):
        """
        Parameters
        ----------

        residue : RESIDUE
            RESIDUE object to be added

        """
        self.residues[residue.index] = residue

    def __getitem__(self, key: int) -> RESIDUE:
        """
        Parameters
        ----------

        key : int
            Used to index the residue. Overloaded for the [] operator.
        
        Raises
        ------
        TypeError
            If the residues are not indexed by an integer.

        Returns
        ------
        RESIDUE
            The indexed residue.
        """
        if not isinstance(key, int):
            raise TypeError("You must index residues by an integer")
        return self.residues[key]

class COMPLEX:
    """
    COMPLEXES are collections of chains that typically bind to one another. 
    Complexes can be indexed by their chain names.
    ...

    Attributes
    ----------

    chains : dict
        A dictionary where the keys are the chain names and values are the chains.
    chainnames : list
        A list of the chain names in the complex.
    filepath : str
        The path pointing to the single-frame PDB file.

    Methods
    -------
    add_chain(chain)
        Used to populate chains into the complex.

    __getitem__(key)
        Used to index the chains in a complex.
    """
    def add_chain(self, chain: CHAIN):
        """
        Parameters
        ----------

        chain : CHAIN
            CHAIN object to be added

        """
        self.chains[chain.name] = chain
        self.chainnames.append(chain.name)

    def __getitem__(self, key: str) -> CHAIN:
        """
        Parameters
        ----------

        key : str
            Used to index the chain. Overloaded for the [] operator.
        
        Raises
        ------
        TypeError
            If the residues are not indexed by a string.

        Returns
        ------
        CHAIN
            The indexed CHAIN.
        """
        if not isinstance(key, str):
            raise TypeError("You must index chains by the chain letter")
        return self.chains[key]

    def __init__(self, pdbfile: str = None, frame_data: list = None):
        """
        We can create a COMPLEX directly from a single pdb file or 
        from the data passed by the TRAJECTORY class. 

        Parameters
        ----------

        pdbfile : None | str
            Used to create the complex. If the value is a string, then pdfile
            is assumed to be the path leading to the PDB file. The file
            is parsed and the appropriate hierarchy of objects is created.

        frame_data : None | list
            Not intended to be used except by the trajectory class. 
            Used to create the complex. When the value is a list,
            the list is parsed and the appropriate hierarchy of 
            objects is created.
        
        """
        self.chains = dict()
        self.chainnames = []
        self.filepath = pdbfile
        
        #Creating from the single pdb file.
        if pdbfile is not None:
            frame_data = []
            with open(pdbfile, "r") as f:
                data = f.readlines()
            for line in data:
                if line[0:4] == "ATOM":
                    #We only care about lines with ATOM entries
                    frame_data.append(line)
        
        #Creating from the multi-frame pdb file.
        elif frame_data is not None:
            tmp_data = []
            for line in frame_data:
                 if line[0:4] == "ATOM":
                    #We only care about lines with ATOM entries
                    tmp_data.append(line)
            frame_data = tmp_data

        #Column defintions from PDB files are used to obtain atom information
        #and create the COMPLEX.
        for line in frame_data:
            atomnum = int(line[6:11].strip())
            atomtype = line[12:16].strip()
            alt_location = line[16]
            resname = line[17:20]
            chainname = line[21]
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            #Creating the ATOM
            cur_atom = ATOM(x, y, z, atomnum, atomtype, alt_location)

            #Making sure a chain exists. If not, create one and add to the complex.
            no_chain = (not bool(self.chains)) or (chainname not in self.chainnames)

            if no_chain:
                self.add_chain(CHAIN(chainname))
            
            cur_chain = self.chains[chainname]

            #Making sure a residue exists. If not, create on and add to the chain.
            no_res = (not bool (cur_chain.residues)) or (resnum not in cur_chain.residues.keys())

            if no_res:
                cur_res = RESIDUE(resname, resnum)
                cur_chain.add_residue(cur_res)

            #Populate atom information into the residue.
            self.chains[chainname].residues[resnum].add_atom(cur_atom)

 
class TRAJECTORY:
    """
    A trajectory is defined by a multi-frame pdb file where each frame represents a protein  
    complex at a specific timestep. Note: trajectories are zero-indexed.

    ...

    Attributes
    ----------

    filename : str
        A dictionary where the keys are the chain names and values are the chains.
    startpoints : list
        A list which contains the line numbers which mark the start of a frame in a multi-pdb file.
    frame : int
        A number which denotes the frame currently in memory 
    complex : COMPLEX
        The COMPLEX object currently in memory within the trajectory.

    Methods
    -------
    read_start_points(multipdbfile)
        Used to obtain the start points in a pdbfile.

    read_frame(start)
        Parses the data for use by the COMPLEX constructor.

    add_complex(start)
        A wrapper for the COMPLEX constructor.

    __getitem__(key)
        Used to index the frames in a trajectory by frame number.

    __len__(key)
        Used get the number of frames in the TRAJECTORY object.
    """     
    #read_start_points() returns the first line number in a multi-frame pdb file for each frame.
    #Doing so ensures that not the whole multipdb file is read into memory when you create a complex.

    def read_start_points(self, multipdbfile: str) -> list:
        """
        Parameters
        ----------

        multipdbfile : str
            The  path leading to the PDB file.
        
        Returns
        -------
        list
            A list of line numbers which mark the starting points for each frame.
        """
        #We assume at least one frame.
        count = 1
        start_points = [1]
        with open(multipdbfile, "r") as f:
            curr_line= f.readline()
            while (curr_line):
                curr_line = f.readline()
                if curr_line[0:3] == "END":
                    #The END flag marks the and of a frame.
                    start_points.append(count+2)
                count += 1
        start_points.pop(-1)
        return start_points


    def read_frame(self, start: int) -> list:
        """
        Parameters
        ----------

        start : int
            The first line to read for a specific frame in a multi-pdb file.
        
        Returns
        -------
        list
            A list where each item is a line from one frame of a PDB file. 
        The line can be parsed by the COMPLEX class.
        """
        data = []
        #We assume at least one frame.
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

    def __init__(self, multipdbfile:str):
        """
        Parameters
        ----------

        multipdbfile : str
            The  path leading to the PDB file. The file
            is parsed and the appropriate hierarchy of objects is created
            automatically for frame 0.
        """
        self.filename = multipdbfile
        self.start_points = self.read_start_points(multipdbfile)
        self.frame = 0
        self.add_complex(self.read_frame(self.frame))

    def add_complex(self, data: list):
        """
        Parameters
        ----------

        data : list
            The list which is used to create the complex from the COMPLEX constructor.
            The complex is loaded into memory as an attribute.
        """
        self.complex = COMPLEX(frame_data = data)
        
    def __len__(self):
        """
        Returns
        -------
        COMPLEX
            The number of frames in the trajectory.
        """
        
        return len(self.start_points)

    def __getitem__(self, key: int) -> COMPLEX:
        """
        Parameters
        ----------

        key : str
            Used to index the trajectory. Overloaded for the [] operator.
        
        Raises
        ------
        TypeError
            If the trajectory is not indexed by an integer.

        Returns
        -------
        COMPLEX
            The indexed complex.
        """
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