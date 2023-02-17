from structure import *
import pandas as pd
import numpy as np

def filter_charged_neutral(results: list, dont_consider=None) -> list:
    keep = []

    columns = ["chain1", "resname1","resnum1", "chain2", "resname2", "resnum2", "distance"]

    for pair in results:
        hit = (pair[1] in CHAIN.charged and pair[4] in CHAIN.neutral) or (pair[1] in CHAIN.neutral and pair[4] in CHAIN.charged)
        if hit:
            keep.append(pair)
    
    df = pd.DataFrame(columns=columns, data=keep)

    if dont_consider is not None:
        keep_dont_consider = np.full(len(df), True)
        df_drop = pd.read_csv(dont_consider)
        for i in range(len(df_drop)):
            row = df_drop.iloc[i]
            c = str(row["chain"])
            rn  = str(row["resname"])
            ri  = int(row["resnum"])

            for j in range(len(df)):
                row2 = df.iloc[j]
                c1 = str(row2["chain1"])
                r1n= str(row2["resname1"])
                r1i = int(row2["resnum1"])
                c2 = str(row2["chain2"])
                r2n = str(row2["resname2"])
                r2i = int(row2["resnum2"])

                if (c1 == c and r1n == rn and r1i == ri):
                    keep_dont_consider[j] = False

                elif ( 2 == c and r2n == rn and r2i == ri):
                    keep_dont_consider[j]= False
        
        return df[keep_dont_consider].reset_index(drop=True)

    return df
