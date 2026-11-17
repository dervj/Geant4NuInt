#!/usr/bin/env python3
"""Trace RAY pseudoparticles (FLUKA manual ch. 14) through the ORIGINAL
lattice v0.8 FLUKA input and record every boundary crossing:
(position, region name, lattice cell name, material name).

This is the FLUKA ground truth for the GDML point/ray checks.  It uses the
unmodified geometry, material and ASSIGNMA cards, and FLUKA's own LATTICE
treatment.

usage: fluka_rays.py <orig.inp> <rays.json> <out.json>
rays.json: [{"name":..., "start":[x,y,z] cm, "dir":[u,v,w]}, ...]
"""
import json
import os
import re
import struct
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

FLUPRO = os.environ.get("FLUPRO", "/scratch/jdervan/flukastudy/fluka4-5.2")


def fmt10(v):
    v = float(v)
    if abs(v) < 1e-9:
        v = 0.0
    # most precise representation that fits a 10-column FLUKA field
    for p in (8, 7, 6, 5, 4, 3):
        for t in (format(v, f".{p}g"),
                  format(v, f".{p}e").replace("e-0", "e-").replace("e+0", "e+")):
            if len(t) <= 10:
                return t
    raise ValueError(v)


def card(kw, w=(), sdum=""):
    s = f"{kw:<10}"
    for i in range(6):
        v = w[i] if i < len(w) else None
        s += f"{'' if v is None else fmt10(v):>10}"
    return s + sdum


def make_input(orig, start, d):
    txt = open(orig).read().splitlines()
    out = []
    for ln in txt:
        if re.match(r"^#define (SOURCE_DECAY|USRDUMP|PHYSICS_DECAY)\b", ln):
            ln = "!" + ln
        elif ln.startswith("BEAM ") and "MUON" in ln:
            ln = card("BEAM", (-0.001,), "RAY")
        elif re.match(r"^BEAMPOS\s*$", ln):
            u, v, w = d
            # FLUKA derives w from u, v: the formatted u, v must satisfy
            # u^2+v^2 <= 1 or it re-normalises and crashes.
            while float(fmt10(u)) ** 2 + float(fmt10(v)) ** 2 > 1:
                u *= 1 - 1e-7
                v *= 1 - 1e-7
            ln = card("BEAMPOS", (start[0], start[1], start[2], u, v),
                      "NEGATIVE" if w < 0 else "")
        elif ln.startswith("START"):
            ln = card("START", (1.0,))
        out.append(ln)
    return "\n".join(out) + "\n"


def records(path):
    data = open(path, "rb").read()
    i = 0
    while i < len(data):
        (n,) = struct.unpack_from("<i", data, i)
        yield data[i + 4 : i + 4 + n]
        i += 8 + n


def region_names(outfile):
    """Region and lattice numbering from FLUKA's .out geometry listing."""
    txt = open(outfile, errors="replace").read()
    regs = {}
    for m in re.finditer(r"^\s*Region n\.\s*(\d+)\s+(\S+)", txt, re.M):
        regs[int(m.group(1))] = m.group(2)
    if not regs:
        for m in re.finditer(r"^\s+(\d+)\s+(\S+)\s+\d+\s+\d+\s*$", txt, re.M):
            regs.setdefault(int(m.group(1)), m.group(2))
    # lattice cell numbers = WHAT(4) of the LATTICE cards, echoed in .out
    lats = {}
    for m in re.finditer(r"NEXT CONTROL CARD \*+\s+LATTICE\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)", txt):
        lats[int(float(m.group(4)))] = m.group(5)
    return regs, lats


def parse_rays(path):
    recs = list(records(path))
    pts = []
    k = 0
    # start record: NRAYRN, MREG, MLATTC(i8), MMAT, EKIN
    nray, mreg = struct.unpack_from("<ii", recs[k], 0)
    mlat = struct.unpack_from("<q", recs[k], 8)[0]
    mmat = struct.unpack_from("<i", recs[k], 16)[0]
    xyz = struct.unpack_from("<3f", recs[k + 1], 0)
    pts.append({"xyz": xyz, "mreg": mreg, "mlat": mlat, "mmat": mmat, "mat": None})
    k += 3
    while k + 3 < len(recs) + 1 and k < len(recs):
        r = recs[k]
        mreg, = struct.unpack_from("<i", r, 0)
        mlat = struct.unpack_from("<q", r, 4)[0]
        mmat, = struct.unpack_from("<i", r, 12)
        mat = r[16:24].decode(errors="replace").strip()
        idisc, = struct.unpack_from("<i", r, 24)
        xyz = struct.unpack_from("<3f", recs[k + 1], 0)
        pts.append({"xyz": xyz, "mreg": mreg, "mlat": mlat, "mmat": mmat,
                    "mat": mat, "idisc": idisc})
        k += 4
        if idisc:
            break
    return pts


def run_ray(orig, ray, workdir):
    d = os.path.join(workdir, ray["name"])
    os.makedirs(d, exist_ok=True)
    inp = os.path.join(d, "ray.inp")
    open(inp, "w").write(make_input(orig, ray["start"], ray["dir"]))
    subprocess.run([f"{FLUPRO}/bin/rfluka", "-N0", "-M1", "ray"], cwd=d,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    f10 = os.path.join(d, "ray001_fort.10")
    out = os.path.join(d, "ray001.out")
    if not os.path.exists(f10):
        return ray["name"], {"error": "no fort.10", "dir": d}
    regs, lats = region_names(out)
    pts = parse_rays(f10)
    for p in pts:
        p["region"] = regs.get(p["mreg"], str(p["mreg"]))
        p["cell"] = lats.get(p["mlat"], None)
    # record i = region traversed from point i-1 to point i (start record =
    # starting region)
    segs = [{"a": pts[i - 1]["xyz"], "b": pts[i]["xyz"], "region": pts[i]["region"],
             "cell": pts[i]["cell"], "mat": pts[i]["mat"]} for i in range(1, len(pts))]
    return ray["name"], {"ray": ray, "segments": segs, "end_blackhole": pts[-1].get("idisc", 0) != 0}


def main():
    orig, raysf, outf = sys.argv[1:4]
    rays = json.load(open(raysf))
    work = os.path.abspath("fluka_ray_runs")
    with ThreadPoolExecutor(24) as ex:
        res = dict(ex.map(lambda r: run_ray(os.path.abspath(orig), r, work), rays))
    json.dump(res, open(outf, "w"), indent=0)
    for n, r in res.items():
        print(n, "error" if "error" in r else len(r["segments"]))


if __name__ == "__main__":
    main()
