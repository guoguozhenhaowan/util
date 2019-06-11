#!/usr/bin/env python3

import sys
import json
import matplotlib.pyplot as plt

if len(sys.argv) <= 2:
    print(sys.argv[0], "<bamqc.json>", "<out.png>")
    sys.exit(0)

jsnFile = sys.argv[1]
pngFile = sys.argv[2]

injson = json.load(open(jsnFile, "r"))
depList = injson["DupIncludeIntegrityDepth"]
totDepth = sum(depList)
normDepList = list(map(lambda x: x*100/totDepth, depList))
fig = plt.figure()
ax = plt.axes()
plt.plot(normDepList)
plt.xlabel("Relative Position in Gene(5'->3')")
plt.ylabel("Percentage of Reads(%)")
plt.subplots_adjust(bottom=.2, left=.2)
plt.savefig(pngFile, dpi=300, framon=True)
