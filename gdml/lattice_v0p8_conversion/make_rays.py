#!/usr/bin/env python3
"""Ray set for the FLUKA-vs-GDML comparison (FLUKA units: cm)."""
import json
import math

rays = []


def add(name, s, d):
    n = math.sqrt(sum(c * c for c in d))
    rays.append({"name": name, "start": list(map(float, s)), "dir": [c / n for c in d]})


# A. along the beamline, from both ends, at transverse offsets (cm)
offs = [(0, 0), (3, 0), (0, -3), (8, 0), (0, 12), (-15, 0), (0, 20), (25, 0),
        (18, 18), (-28, 5), (-40, 0), (0, 60), (100, 0), (-150, 30), (0, -120),
        (-51, 160), (-200, 5), (250, 0), (0, 190), (600, -300)]
for x, y in offs:
    add(f"zup_{x}_{y}", (x, y, -24890), (0, 0, 1))
    add(f"zdn_{x}_{y}", (x, y, 1495), (0, 0, -1))

# B. transverse through each lattice cell centre, plus drifts/tunnel
zc = {"IB201": -6570, "IB102": -5740, "IB203": -4910, "IQF204": -4280,
      "IQF2A05": -3650, "IQD107": -2870, "IQD1A09": -1940, "IQF1B11": -1360,
      "IQF1A13": -1080, "IQF115": -750, "IP100": 0, "IQF1d17": 750,
      "IQF1Ad19": 1080, "IQF1Bd21": 1360, "drift_-6255": -6255,
      "tun_-20000": -20000, "tun_-10000": -10000, "drift_-1470": -1470,
      "IP100_z300": 300, "IP100_z-450": -450, "drift_1470": 1470}
for c, z in zc.items():
    for y in (0, 4, 17):
        add(f"x_{c}_y{y}", (-1195, y, z + 1), (1, 0, 0))
    for x in (0, 6):
        add(f"y_{c}_x{x}", (x, -1495, z + 1), (0, 1, 0))
    add(f"xy_{c}", (-1195, -1195, z - 1), (1, 1, 0))
    add(f"xmy_{c}", (1195, -1195, z + 2), (-1, 1, 0))

# C. IP fan from the beam axis: through DET_CHAM, BP_VAC, liners, nozzles
for z0 in (0, 50, 150, 300, 450, 580, -100, -300, -560):
    for th in (2, 5, 10, 20, 45, 90, 135, 160, 175):
        for ph in (0, 90, 210):
            t, p = math.radians(th), math.radians(ph)
            add(f"ip_{z0}_{th}_{ph}", (0, 0, z0),
                (math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)))

json.dump(rays, open("rays.json", "w"), indent=0)
print(len(rays), "rays")
