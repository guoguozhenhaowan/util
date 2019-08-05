#!/usr/bin/env python3

from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import pandas as pd
import sys
import os

if len(sys.argv) <= 2:
    print(sys.argv[0], "<anarpt1.xlsx>", "<anarpt2.xlsx>")
    sys.exit(0)
anarpt1 = sys.argv[1]
anarpt2 = sys.argv[2]
# get sample name
aName = os.path.basename(anarpt1).replace(".report.xlsx", "")
bName = os.path.basename(anarpt2).replace(".report.xlsx", "")
# read df
dfa = pd.read_excel(anarpt1, sheet_name="AllGeneExp")
dfb = pd.read_excel(anarpt2, sheet_name="AllGeneExp")
# intersection 
both = dfa.merge(dfb, how="inner", on="Gene", suffixes=(aName, bName))
both.to_csv("Gene.In.Both.csv", index=False)
# only a
onlya = dfa[~dfa["Gene"].isin(both["Gene"])]
onlya.to_csv(aName + ".csv", index=False)
# only b
onlyb = dfb[~dfb["Gene"].isin(both["Gene"])]
onlyb.to_csv(bName + ".csv", index=False)
# venn
gena = set(dfa.Gene)
genb = set(dfb.Gene)
plt.figure()
v = venn2([gena, genb], set_labels=(aName, bName))
plt.savefig("Venn.png", dpi=300)
