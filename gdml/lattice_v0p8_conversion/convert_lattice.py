#!/usr/bin/env python3
"""Convert the MAIA lattice v0.8 FLUKA input to GDML with the LATTICE cells
instantiated.  Replaces the ad-hoc interactive pyg4ometry conversion of
2026-07-28 (mucoll_mdi_lattice_v0p8.gdml, 10 m world, no lattice copies).

What this does, and why (see GDML_LATTICE_REPAIR_HANDOFF.md):

1. World.  FLUKA's physical geometry is enclosed by BLCKHOLE regions
   (CENTRAL/OUTACCEL/CCENTRAL at |x| > ~12 m, AIR_000 at z < -249 m,
   AIR_031/AIR_032 at z > 14.77 m, INNERr/OUTERr outside body `cont`;
   `cont` already stops at y = -15 m).  We shrink body `cont` to the box
   CONT_MM below and intersect every beamline region with it, so all solids
   are finite.  Inside that box the FLUKA geometry is unchanged.  The only
   physics truncation is the rock above y = +15 m (FLUKA has it up to
   y = +300 m).  Leaving the box is equivalent to entering a FLUKA black
   hole: the box is convex and only vacuum lies outside it.  The Geant4
   world is centred on the FLUKA origin, so it is CONT_MM made symmetric
   about 0 (+1 mm).  FLUKA coordinates are kept unchanged.

1b. pyg4ometry's half-space filter and its mesh-based solid minimisation
   are disabled (see the comments at the imports): the first deletes
   bounding planes from thin regions, the second drops real regions.

2. LATTICE.  pyg4ometry 1.4.4 parses LATTICE, but its _convertLatticeCells
   places every prototype at the *cell* centre and rotation, ignoring each
   prototype's own offset.  It also gives all copies the same PV name.  The
   call is wrapped in `except UnboundLocalError: pass`.  We bypass it: the
   lattices are taken out of the registry before conversion, then for each
   cell we build a vacuum LV from the cell region (the cell body as FLUKA
   defines it, via $start_transform -CELL) and place it in the world.  The
   prototype LVs go inside it as daughters with
       M_daughter = M_cell^-1 . T^-1 . M_proto(parking)
   where T is the cell's ROT-DEFI chain (cell frame -> prototype frame,
   FLUKA manual LATTICE) as implemented by pyg4ometry's RotoTranslation.
   The prototypes in each cell are those whose parking AABB lies inside the
   cell's AABB after mapping it with T.  The parking copies and PARKr are
   removed.

3. Black holes.  BLCKHOLE regions that survive inside the world (DET_LAT,
   plus clipped slivers of the outer ones) get a dedicated vacuum material
   named BLCKHOLE.  Every LV gets GDML auxiliaries FlukaRegion and
   FlukaMaterial, and black holes also get FlukaBlackHole=1.  Cell LVs get
   FlukaLatticeCell=<cell name>.  The Geant4 application kills/scores by
   these.
"""
import argparse
import json
import re
import sys
import time

import numpy as np
import pyg4ometry
from pyg4ometry import geant4 as g4
from pyg4ometry import transformation as tr
from pyg4ometry.gdml.Defines import Auxiliary
from pyg4ometry.convert.fluka2Geant4 import fluka2Geant4
import pyg4ometry.fluka.body  # noqa: F401  (make sure the module is loaded)

# NB: `import pyg4ometry.convert.fluka2Geant4 as m` binds the *function* (the
# package re-exports it under the module's name), so patch via sys.modules.
_f2g = sys.modules["pyg4ometry.convert.fluka2Geant4"]

# pyg4ometry 1.4.4 _filterHalfSpaces drops any plane farther from the
# region AABB than its corner distance.  That is only valid when the other
# bodies already bound the region, and they don't for shared infinite bodies:
# e.g. BEAM_V5D (+T_IBP5 -C_BP4D +C_BP5D +bppla028 -bppla011, a 10 cm slice
# at z=14.7 m) lost -C_BP4D and -bppla011 and became a 32 m tube running
# through the IP and all IR quadrupole cells.  The old GDML put (0,0,0) in
# BEAM_V5D for this reason.  Keep all planes.
_f2g._filterHalfSpaces = lambda flukareg, regionZoneAABBs: flukareg

# pyg4ometry's minimiseSolids pre-pass meshes every zone to find null zones
# and to size the "infinite" bodies.  The arc tunnel is built from YCC
# cylinders of radius ~2.25 km (cPROTc02/cPROTp02); their polygon meshes are
# off by metres, so whole 70 cm tunnel slices (T00x_002, AIR_002, ...) were
# dropped as "null" and others were truncated.  The Geant4 solids themselves
# are exact (G4Tubs etc.), so we skip the pre-pass and bound infinite bodies
# at 1 km instead of 50 km: that still covers the 250 m world and keeps G4
# surface arithmetic far inside kCarTolerance.  (minimiseSolids=False itself
# crashes in 1.4.4 with UnboundLocalError, so we stub the steps instead: every
# zone is kept, no body gets a minimum AABB, and pyg4ometry's own lattice
# step -- which we replace anyway -- is a no-op.)
_fbody = sys.modules["pyg4ometry.fluka.body"]
_fbody.INFINITY = 1.0e6  # mm
_f2g._getRegionZoneAABBs = lambda flukareg, regions, q: {
    n: [True] * len(flukareg.regionDict[n].zones) for n in regions}
_f2g._makeBodyMinimumAABBMap = lambda flukareg, regionZoneAABBs, regions: None
_f2g._convertLatticeCells = lambda *a, **k: None

# cont clip box in cm (FLUKA units): x, y, z min/max.
CONT_CM = (-1250.0, 1250.0, -1500.0, 1500.0, -24950.0, 1500.0)
OMIT = ["PARKr", "INNERr", "OUTERr"]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def preprocess(src, dst):
    """Shrink body `cont` in the FLUKA text.  Everything else verbatim."""
    txt = open(src).read()
    pat = re.compile(r"^RPP cont .*$", re.M)
    assert len(pat.findall(txt)) == 1
    new = "RPP cont       {} {} {} {} {} {}".format(*CONT_CM)
    txt = pat.sub(new, txt)
    open(dst, "w").write(txt)


def active(pv):
    """Active (local -> mother) 4x4 matrix of a pyg4ometry PV.  pyg4ometry
    stores the passive TB angles: rot = reverse(tbxyz(R_active))."""
    ang = [float(v) for v in pv.rotation.eval()]
    R = tr.tbxyz2matrix(ang).T
    M = np.identity(4)
    M[:3, :3] = R
    M[:3, 3] = [float(v) for v in pv.position.eval()]
    if pv.scale is not None:
        s = [float(v) for v in pv.scale.eval()]
        assert np.allclose(s, 1), f"reflected placement {pv.name}"
    return M


def matrix(R, t):
    M = np.identity(4)
    M[:3, :3] = R
    M[:3, 3] = t
    return M


def aabb(lv, M):
    lo, hi = lv.mesh.getBoundingBox(M[:3, :3], list(M[:3, 3]))
    return np.array(lo, float), np.array(hi, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--report", default="convert_report.json")
    args = ap.parse_args()

    pre = args.out.replace(".gdml", "_cont.inp")
    preprocess(args.inp, pre)
    log("preprocessed ->", pre)

    fr = pyg4ometry.fluka.Reader(pre).getRegistry()
    lattices = dict(fr.latticeDict)
    fr.latticeDict.clear()
    log(len(fr.regionDict), "regions,", len(lattices), "lattices")

    cont = [b for b in fr.bodyDict.values() if b.name == "cont"][0]

    # Prototype regions: the regions FLUKA places in the parking area.  These
    # are the regions listed between PARKr and the first beamline region
    # (CENTRAL) in the input.
    names = list(fr.regionDict)
    proto = names[names.index("PARKr") + 1 : names.index("CENTRAL")]
    world_regions = [n for n in names if n not in proto and n not in OMIT]
    log(len(proto), "prototype regions,", len(world_regions), "world regions")

    # Clip every world region to the (shrunk) cont box.  This is a no-op
    # for regions that already carry +cont and a FLUKA no-op inside cont for
    # the rest (every non-OUTERr/INNERr/PARKr point lies inside cont).
    nclip = 0
    for n in world_regions:
        for z in fr.regionDict[n].zones:
            if cont not in z.bodies():
                z.addIntersection(cont)
                nclip += 1
    log("added +cont to", nclip, "zones")

    blackholes = sorted(n for n, m in fr.assignmas.items() if m[0] == "BLCKHOLE")
    matof = {n: fr.assignmas[n][0] for n in fr.assignmas}

    half = [max(abs(CONT_CM[0]), abs(CONT_CM[1])) * 10,
            max(abs(CONT_CM[2]), abs(CONT_CM[3])) * 10,
            max(abs(CONT_CM[4]), abs(CONT_CM[5])) * 10]
    worlddim = [2 * h + 1.0 for h in half]
    log("world full lengths (mm):", worlddim)

    t0 = time.time()
    greg = fluka2Geant4(fr, omitRegions=OMIT, omitBlackholeRegions=False,
                        worldDimensions=worlddim)
    log(f"fluka2Geant4 done in {time.time()-t0:.0f} s")
    assert _f2g._filterHalfSpaces.__name__ == "<lambda>" and _fbody.INFINITY == 1.0e6

    wlv = greg.getWorldVolume()
    pvs = {pv.name[:-3]: pv for pv in wlv.daughterVolumes}  # strip "_pv"
    converted = set(pvs)
    dropped = sorted(set(names) - converted - set(OMIT))
    log("regions dropped as null by pyg4ometry:", dropped)

    # ---- lattice cells ------------------------------------------------------
    protoM = {}
    protoBB = {}
    for n in proto:
        if n not in pvs:
            continue
        protoM[n] = active(pvs[n])
        protoBB[n] = aabb(pvs[n].logicalVolume, protoM[n])

    report = {"world_mm": worlddim, "cont_cm": CONT_CM, "omitted": OMIT,
              "dropped_null": dropped, "cells": {}}
    used = set()
    tol = 1.0  # mm
    for cname, lat in lattices.items():
        T = lat.getTransform().to4DMatrix()
        creg = lat.cellRegion
        Mc = matrix(creg.rotation(), list(creg.centre()))
        csolid = creg.geant4Solid(greg)
        clv = g4.LogicalVolume(csolid, "G4_Galactic", f"{cname}_lv", greg)
        clv.addAuxiliaryInfo(Auxiliary("FlukaLatticeCell", cname, greg, addRegistry=False))
        clv.addAuxiliaryInfo(Auxiliary("FlukaRegion", cname, greg, addRegistry=False))
        # cell AABB in the parking frame
        lo, hi = aabb(clv, T @ Mc)
        members = [n for n, (plo, phi) in protoBB.items()
                   if np.all(plo >= lo - tol) and np.all(phi <= hi + tol)]
        Minv = np.linalg.inv(Mc) @ np.linalg.inv(T)
        for n in members:
            Md = Minv @ protoM[n]
            rot = tr.matrix2tbxyz(Md[:3, :3].T)
            g4.PhysicalVolume(list(rot), list(Md[:3, 3]), pvs[n].logicalVolume,
                              f"{n}_in_{cname}_pv", clv, greg)
        used.update(members)
        g4.PhysicalVolume(list(tr.matrix2tbxyz(Mc[:3, :3].T)), list(Mc[:3, 3]),
                          clv, f"{cname}_pv", wlv, greg)
        wlo, whi = aabb(clv, Mc)
        report["cells"][cname] = {"members": members,
                                  "aabb_world_mm": [wlo.tolist(), whi.tolist()],
                                  "T": T.tolist()}
        log(f"cell {cname}: {len(members)} prototypes, world AABB z "
            f"[{wlo[2]/1000:.3f}, {whi[2]/1000:.3f}] m")

    unused = sorted(set(protoM) - used)
    report["unused_prototypes"] = unused
    log("prototype regions in no cell:", unused)

    # remove parking placements from the world
    for n in proto:
        if n in pvs:
            pv = pvs[n]
            wlv.daughterVolumes.remove(pv)
            del greg.physicalVolumeDict[pv.name]
            wlv._daughterVolumesDict.pop(pv.name, None)
    for n in unused:
        lv = pvs[n].logicalVolume
        greg.logicalVolumeDict.pop(lv.name, None)

    # ---- materials / tags ---------------------------------------------------
    bh = g4.MaterialSingleElement("BLCKHOLE", 1, 1.008, 1e-25, greg)
    for lvname, lv in list(greg.logicalVolumeDict.items()):
        reg = lvname[:-3] if lvname.endswith("_lv") else lvname
        if reg not in matof or lv is wlv:
            continue
        lv.addAuxiliaryInfo(Auxiliary("FlukaRegion", reg, greg, addRegistry=False))
        lv.addAuxiliaryInfo(Auxiliary("FlukaMaterial", matof[reg], greg, addRegistry=False))
        if matof[reg] == "BLCKHOLE" and reg not in lattices:
            lv.material = bh
            lv.addAuxiliaryInfo(Auxiliary("FlukaBlackHole", "1", greg, addRegistry=False))
    report["blackholes_kept"] = sorted(
        r for r in blackholes if f"{r}_lv" in greg.logicalVolumeDict)
    log("black holes kept:", report["blackholes_kept"])

    w = pyg4ometry.gdml.Writer()
    w.addDetector(greg)
    w.write(args.out)
    log("wrote", args.out)
    json.dump(report, open(args.report, "w"), indent=1)


if __name__ == "__main__":
    main()
