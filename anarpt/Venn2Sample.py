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
# stat
dfa = pd.read_excel(anarpt1, sheet_name="AllGeneExp")
dfb = pd.read_excel(anarpt2, sheet_name="AllGeneExp")
gena = set(dfa.Gene)
genb = set(dfb.Gene)
both = pd.DataFrame({"Gene": list(gena.intersection(genb))})
onlya = pd.DataFrame({"Gene": list(gena.difference(genb))})
onlyb = pd.DataFrame({"Gene": list(genb.difference(gena))})
both.to_csv("Gene.In.Both.csv", index=False)
aName = os.path.basename(anarpt1).replace(".report.xlsx", "")
bName = os.path.basename(anarpt2).replace(".report.xlsx", "")
onlya.to_csv(aName + ".csv", index=False)
onlyb.to_csv(bName + ".csv", index=False)
# plot
plt.figure()
v = venn2([gena, genb], set_labels=(aName, bName))
plt.savefig("Venn.png", dpi=300)
