import pandas as pd

df = pd.read_csv("hs_kinesins.csv")

df_choose = df[(df.PDB == 1) & (df.Complex == 1) & (df.Resolution == 1)]

df_choose.to_csv("chosen.csv", index=False)