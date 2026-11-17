#!/usr/bin/env python3
"""Compare FLUKA RAY segments with Geant4 point location in a GDML file.

For every FLUKA segment inside the clip box we test the midpoint and two
points just inside its ends.  Geant4 must report LV '<region>_lv'.  The
mother must be '<cell>_lv' for segments inside a lattice cell, and the
world otherwise.

usage: compare_rays.py fluka_rays.json file.gdml [--box x0 x1 y0 y1 z0 z1 (cm)]
"""
import argparse
import collections
import json
import os
import subprocess
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
IMAGE = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"
VIEW = "/cvmfs/sft.cern.ch/lcg/views/LCG_110a/x86_64-el9-gcc13-opt"

ap = argparse.ArgumentParser()
ap.add_argument("rays")
ap.add_argument("gdml")
ap.add_argument("--box", nargs=6, type=float,
                default=[-1250, 1250, -1500, 1500, -24950, 1500])
ap.add_argument("--out", default="compare_result.json")
a = ap.parse_args()
box = np.array(a.box).reshape(3, 2)

res = json.load(open(a.rays))
pts, meta = [], []
for rname, r in res.items():
    if "segments" not in r:
        print("FLUKA error for", rname, r)
        continue
    for s in r["segments"]:
        A, B = np.array(s["a"], float), np.array(s["b"], float)
        L = np.linalg.norm(B - A)
        if L < 1e-3:
            continue
        u = (B - A) / L
        eps = min(0.05, 0.25 * L)
        for tag, p in (("a", A + eps * u), ("m", 0.5 * (A + B)), ("b", B - eps * u)):
            if np.all(p >= box[:, 0]) and np.all(p <= box[:, 1]):
                pts.append(p * 10.0)  # cm -> mm
                meta.append((rname, tag, s, L))

with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, dir=HERE) as f:
    for p in pts:
        f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    ptsfile = f.name
cmd = (f"source {VIEW}/setup.sh >/dev/null 2>&1; "
       f"{HERE}/g4tools/gdmlprobe {os.path.abspath(a.gdml)} {ptsfile}")
out = subprocess.run(["apptainer", "exec", "-B", "/cvmfs", "-B", "/scratch", IMAGE,
                      "bash", "-c", cmd], capture_output=True, text=True).stdout
os.unlink(ptsfile)
g4 = {}
for ln in out.splitlines():
    f = ln.split()
    if len(f) == 5 and f[0].isdigit():
        g4[int(f[0])] = f[1:]

nok = 0
bad = []
regmis = collections.Counter()
for i, (rname, tag, s, L) in enumerate(meta):
    lv, mat, mo, depth = g4.get(i, ["<missing>", "-", "-", "-"])
    exp_lv = s["region"] + "_lv"
    exp_mo = (s["cell"] + "_lv") if s["cell"] else "wl"
    if lv == exp_lv and mo == exp_mo:
        nok += 1
    else:
        bad.append({"ray": rname, "tag": tag, "xyz_cm": list(pts[i] / 10),
                    "fluka": [s["region"], s["cell"], s["mat"]], "seglen_cm": L,
                    "g4": [lv, mo, mat]})
        regmis[(s["region"], s["cell"], lv, mo)] += 1

print(f"points tested: {len(meta)}   agree: {nok}   disagree: {len(bad)}")
fl_regions = collections.Counter(m[2]["region"] for m in meta)
print(f"distinct FLUKA regions sampled: {len(fl_regions)}; cells sampled: "
      f"{len(set(m[2]['cell'] for m in meta if m[2]['cell']))}")
for k, v in regmis.most_common(40):
    print(f"  {v:5d}  FLUKA {k[0]}@{k[1]}  ->  G4 {k[2]} in {k[3]}")
json.dump({"n": len(meta), "agree": nok, "bad": bad,
           "fluka_regions": dict(fl_regions)}, open(a.out, "w"), indent=0)
