#!/usr/bin/env python3
"""Open-body (cutaway) version of the turbojet: static shells are cut in half
(z > 0 removed) so the compressor, combustor, turbine, shaft etc. are visible.
Rotating parts stay whole.   python build_cutaway.py [outdir]"""
import sys
import cadquery as cq
import build_turbojet as B

CUT = ["intake_shell", "casing", "combustor_outer", "combustor_dome", "diffuser", "ngv"]

def build_cutaway():
    parts = B.build()
    cutter = cq.Workplane("XY").box(400, 300, 150, centered=(False, True, False)).translate((-50, 0, 0))
    out = {}
    for n, (wp, rgb) in parts.items():
        out[n] = (wp.cut(cutter), rgb) if n in CUT else (wp, rgb)
    return out

if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    parts = build_cutaway()
    for n, (wp, _) in parts.items():
        assert all(s.isValid() for s in wp.vals()), n
    B.export(parts, f"{d}/turbojet_cutaway.step", f"{d}/turbojet_cutaway.stl")
    print("exported cutaway")
