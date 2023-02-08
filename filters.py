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
            c1f = str(row["chain1"])
            r1nf  = str(row["resname1"])
            r1if  = int(row["resnum1"])
            c2f  = str(row["chain2"])
            r2nf  = str(row["resname2"])
            r2if  = int(row["resnum2"])

            for j in range(len(df)):
                row2 = df.iloc[j]
                c1 = str(row2["chain1"])
                r1n= str(row2["resname1"])
                r1i = int(row2["resnum1"])
                c2 = str(row2["chain2"])
                r2n = str(row2["resname2"])
                r2i = int(row2["resnum2"])

                if (c1 == c1f and r1n == r1nf and r1i == r1if and c2 == c2f and r2n == r2nf and r2i == r2if):
                    keep_dont_consider[j] = False

                elif (c1 == c2f and r1n == r2nf and r1i == r2if and c2 == c1f and r2n == r1nf and r2i == r1if):
                    keep_dont_consider[j]= False
        
        return df[keep_dont_consider].reset_index(drop=True)

    return df
