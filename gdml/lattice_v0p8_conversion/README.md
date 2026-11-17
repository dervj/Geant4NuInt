# MAIA lattice v0.8: FLUKA → GDML conversion and validation (2026-11-17)

Output: `../mucoll_mdi_lattice_v0p8.gdml` (1.28 MB, sha256 0daf49b1d7e7…), from
`lattice_v_0p8_EU_BIB_.inp` (mucoll-mdi-simulation-lattice-0.8, `fluka_inputs/`) with `python3 convert_lattice.py lattice_v_0p8_EU_BIB_.inp mucoll_mdi_lattice_v0p8.gdml`
(pyg4ometry 1.4.4, 90 s).  It replaces the 2026-07-28 conversion (commit 110076c).

## What the file contains
- World: centred box, full lengths 25.0 × 30.0 × 499.0 m.  FLUKA body `cont` is clipped to
  x ±12.5 m, y ±15 m, z −249.5…+15 m.  Outside that box the FLUKA geometry is black hole,
  except the rock above y = +15 m, which is dropped.  Coordinates are FLUKA's, in mm.
- All 343 FLUKA regions minus PARKr/INNERr/OUTERr.  pyg4ometry dropped none of them.
- 14 LATTICE cells as real volumes (`<cell>_pv`, vacuum LV `<cell>_lv`), filled with
  their prototype regions (`<region>_in_<cell>_pv`): 16 per magnet, 7 in IP100 (LW/RW nozzles,
  BORETH liners, Be pipe DET_CHAM, BP_VAC, DET_LAT).  There are no parking copies.
- Black holes inside the world: material `BLCKHOLE` plus aux `FlukaBlackHole=1`.  These are
  AIR_000, AIR_031, AIR_032, CENTRAL, CCENTRAL, DET_LAT and OUTACCEL.
- Every LV carries the aux `FlukaRegion` and `FlukaMaterial`; cells also carry `FlukaLatticeCell`.
- 275 materials, FLUKA names kept (INERM180, BORETH, NbTi*, LHe_*, ...).

## pyg4ometry 1.4.4 problems worked around (see convert_lattice.py)
1. `_convertLatticeCells` puts every prototype at the cell centre and gives all copies one name,
   and its errors are silently swallowed.  We place the copies ourselves with M_cell⁻¹·T⁻¹·M_proto.
2. `_filterHalfSpaces` deletes bounding planes of thin regions.  BEAM_V5D, a 10 cm slice,
   became a 32 m tube through the IP.
3. The mesh-based "minimise solids" pre-pass gets the 2.25 km arc-tunnel cylinders wrong by
   metres.  It dropped 35 real regions as null (tunnel slices `_002`, AIR_002, BEAM_*1/2/C, ...)
   and truncated others.  We skip it and bound infinite bodies at 1 km (default 50 km).
4. Patches must go through `sys.modules[...]`: `import pyg4ometry.convert.fluka2Geant4 as m`
   binds the function, not the module.

## Validation
- **FLUKA ground truth.** We traced 430 RAY pseudo-particles through the *original* input with
  fluka4-5.2 (`fluka_rays.py`): beam-axis-parallel rays, transverse rays through every cell, and an
  IP fan.  Each FLUKA segment was sampled at three points with the Geant4 navigator
  (`compare_rays.py`, `g4tools/gdmlprobe`), which must return the same region and the same
  lattice cell (or the world).
  - 14 982 points across 240 regions and 14 cells: **14 923 agree (99.6 %)**.
  - All 59 others lie exactly on a FLUKA surface: the rays run along the r = 3/12 cm pipe
    surfaces, the y = −120 cm plane, the z = ±100 cm nozzle planes, or tangent to a plane tilted
    1.45 mrad.  Every one agrees after a ≤0.1 mm displacement.
  - Before fix 3 the agreement was 96.3 %.
- **Overlaps and gaps** (`g4tools/gdmlinside`, Inside() sampling, 2×10⁶ points per box over the
  IR core, the hall and the whole world): 0 overlapping pairs, 0 points inside a cell but in no
  prototype, 0 points in the world but in no region.  We did not use G4 CheckOverlaps:
  GetPointOnSurface fails on most of these booleans.
- **Structure** (`g4tools/gdmlcheck`): 203 world daughters, none outside the world.  Two
  bounding boxes are flagged "protruding" (CENTRAL, CCENTRAL: loose G4MultiUnion bounds from the
  arc cylinders).  Sampling 8×10⁵ points just outside all world faces found no solid there.
- **Handoff probes.** (0,0,±1 m) and (0, 1 m, 1 m) sit exactly on the nozzle planes pLn/pRn.
  ±0.5 mm away they give BP_VAC on the axis, LW/RW (INERM180) at r = 2 and 10 cm, and DET_LAT at
  y = 1 m, as in FLUKA.

## Caveats for users
- Kill tracks on entering any `FlukaBlackHole` volume, as FLUKA does.
- pyg4ometry's length safety leaves ~1e-6 mm slivers of mother vacuum at internal zone
  boundaries.  They are harmless for tracking, but a point sitting exactly on such a plane
  locates to the mother.
- Magnetic fields (the MGNFIELD/`*_2` ROT-DEFIs) are not in GDML.

## Reproducing
- Converter: `convert_lattice.py` needs pyg4ometry 1.4.4. It monkeypatches internals of that
  version, so pin the version.
- FLUKA rays: `make_rays.py` writes the ray list, and `fluka_rays.py orig.inp rays.json out.json`
  traces it. FLUKA is set by `FLUPRO` at the top of the script; fluka4-5.2 was used.
- Comparison: `compare_rays.py out.json ../mucoll_mdi_lattice_v0p8.gdml`.
- Geant4 tools: `g4tools/build.sh` builds `gdmlcheck`, `gdmlprobe` and `gdmlinside` against
  Geant4 11.4.2 (LCG_110a via apptainer). Usage:
  - `gdmlcheck file.gdml`
  - `gdmlprobe file.gdml points_mm.txt`
  - `gdmlinside file.gdml N x0 x1 y0 y1 z0 z1 [seed]`
