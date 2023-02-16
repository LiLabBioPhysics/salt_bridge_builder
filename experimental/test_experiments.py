import sys
sys.path.append("..")

from structure import *
from distances import *
from filters import *

c_6bbn = COMPLEX("kinesin_folders/2A/6bbn_modelled.pdb")

r_6bbn = within_dist(c_6bbn, "A", "E", 7, alpha_dist)
f_6bbn = filter_charged_neutral(r_6bbn)
print(f_6bbn)