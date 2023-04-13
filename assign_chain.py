with open("6ta4_k_mut.pqr", "r") as f:
    data = f.readlines()

output=""

for line in data:
    if line[0:4]=="ATOM":
        output=output+line[0:21] + "K" + line[22:]
    else:
        output=output+line

with open("6ta4_k_mut_2.pqr", "w") as f:
    f.write(output)
