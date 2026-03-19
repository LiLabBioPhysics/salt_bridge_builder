"""Structure Constants"""

POSATOM = {"ND1", "NE2", "NE", "NH1", "NH2", "NZ"}
NEGATOM = {"OD1", "OD2", "OE1", "OE2"}

CTERATOMS = {"C"}
NTERATOMS = {"HN", "N", "HT3"}
OTERATOMS = "OXT"
ALPHAATOM = "CA"

AA = ["ALA", "GLY", "ILE", "LEU", "PRO", "VAL", 
      "PHE", "TRP", "TYR", 
      "SER", "THR", 
      "CYS", "MET",
      "ASN", "GLN",
      "ASP", "GLU", #Negative
      "ARG", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "LYS" #Positive
]

AA1 = ["A", "G", "I", "L", "P", "V", 
       "F", "W", "Y",
       "S", "T",
       "C", "M",
       "N", "Q",
       "D", "E", #Negative
       "R", "H", "H", "H", "H", "H", "H", "H" "K"] #Positive

AA_3_1 = dict(zip(AA, AA1))
AA_1_3 = dict(zip(AA1, AA))

NEUTRAL = set(AA[0:15])
NEGATIVE = set(AA[15:17])
POSITIVE = set(AA[17:])
CHARGED = NEGATIVE | POSITIVE

HEAVY_ATOM_POS = {
    AA[0]:["N", "C", "O", "CA", "CB"], #ALA
    AA[1]:["N", "C", "O", "CA"], #GLY
    AA[2]:["N", "C", "O", "CA", "CB", "CG1", "CG2", "CD1"], #ILE
    AA[3]:["N", "C", "O", "CA", "CB", "CG", "CD1", "CD2"], #LEU
    AA[4]:["N", "C", "O", "CA", "CB", "CG", "CD"], #PRO
    AA[5]:["N", "C", "O", "CA", "CB", "CG1", "CG2"], #VAL
    AA[6]:["N", "C", "O", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"], #PHE
    AA[7]:["N", "C", "O", "CA", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"], #TRP
    AA[8]:["N", "C", "O", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"], #TYR
    AA[9]:["N", "C", "O", "CA", "CB", "OG"], #SER
    AA[10]:["N", "C", "O", "CA", "CB", "OG1", "CG2"], #THR
    AA[11]:["N", "C", "O", "CA", "CB", "SG"], #CYS
    AA[12]:["N", "C", "O", "CA", "CB", "CG", "SD", "CE"], #MET
    AA[13]:["N", "C", "O", "CA", "CB", "CG", "OD1", "ND2"], #ASN
    AA[14]:["N", "C", "O", "CA", "CB", "CG", "CD", "OE1", "NE2"], #GLN
    AA[15]:["N", "C", "O", "CA", "CB", "CG", "OD1", "OD2"], #ASP
    AA[16]:["N", "C", "O", "CA", "CB", "CG", "CD", "OE1", "OE2"], #GLU
    AA[17]:["N", "C", "O", "CA", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"], #ARG
    AA[18]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HIS
    AA[19]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HID
    AA[20]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HIE
    AA[21]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HIP
    AA[22]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HSD
    AA[23]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HSE
    AA[24]:["N", "C", "O", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], #HSP
    AA[25]:["N", "C", "O", "CA", "CB", "CG", "CD", "CE", "NZ"], #LYS
}

HEAVY_ATOM_POS_ALT = {
    AA[2]: {"CD": 7}
}

HEAVY_ATOM_NAMES = set()
for aa in HEAVY_ATOM_POS:
    atom_types = HEAVY_ATOM_POS[aa]
    for a in atom_types:
        HEAVY_ATOM_NAMES.add(a)

SIDE_CHAIN_START = 3
ALPHA_CARBON_POS = 3

"""Distance Constants"""
INTERFACE_CUTOFF = 8