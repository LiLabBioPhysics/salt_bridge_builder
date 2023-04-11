import subprocess
import pandas as pd
from structure import *


def write_script(result, pdbfile):

    chain1 = result[0]
    resname1 = result[1]
    resnum1 = result[2]

    chain2 = result[3]
    resname2 = result[4]
    resnum2 = result[5]

    if resname1 in CHAIN.positive:
        commands="""open {}
swapaa #1/{}:{} {}
save mutation.pdb
        """.format(pdbfile, chain2, resnum2, "GLU")

    elif resname1 in CHAIN.negative:
        commands="""open {}
swapaa #1/{}:{} {}
save mutation.pdb
        """.format(pdbfile, chain2, resnum2, "LYS")


    elif resname2 in CHAIN.positive:
        commands="""open {}
swapaa #1/{}:{} {}
save mutation.pdb
        """.format(pdbfile, chain1, resnum1, "GLU")

    elif resname2 in CHAIN.negative:
        commands="""open {}
swapaa #1/{}:{} {}
save mutation.pdb
        """.format(pdbfile, chain1, resnum1, "LYS")


    else:
        raise ValueError("There is no charged residue in the results.")

    with open("commands.cxc", "w") as f:
        f.write(commands)


def run_script():
    subprocess.run(["chimerax", "--nogui", "--exit", "--silent", "commands.cxc"])

    subprocess.run(["rm", "commands.cxc"])

def clean_pdb(pdbfile):
    output=""
    with open(pdbfile) as f:
        data = f.readlines()
    for line in data:
        if line[0:4]=="ATOM":
            output+=line
    with open("mutation.pdb", "w") as f:
        f.write(output[0:-1])



# df = pd.read_csv("total_results.csv")
# write_script(df.iloc[0,:], "6ta4.pdb")
# run_script()
clean_pdb("mutation.pdb")
