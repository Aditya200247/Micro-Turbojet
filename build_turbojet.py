#!/usr/bin/env python3
"""
Micro turbojet, KJ66 class (~110 mm OD, ~234 mm long, ~100 N thrust class)
Parametric CadQuery model.  Units: mm.  Engine axis = +X.
Inlet lip at x = 0, nozzle exit at x = 234.

    python build_turbojet.py        ->  turbojet.step (assembly, colours + names)
                                        turbojet.stl  (single mesh)

Layout along the axis
    0-49     intake bellmouth + compressor shroud
    23-58    centrifugal compressor wheel (8 main + 8 splitter blades), O.D. 66
    49-56    radial vaned diffuser (15 vanes)
    59-181   bearing tunnel + 2 bearings + shaft (O.D. 8)
    70-96    fuel ring, feed stems, 12 vaporiser tubes
    80-172   annular combustor (outer/inner liner, dome, primary + dilution holes)
    172-180  nozzle guide vanes (13)
    182-196  axial turbine wheel (12 twisted blades), O.D. 67.2
    196-234  exhaust cone, casing converges into the nozzle
"""
import math
import sys

import numpy as np
import cadquery as cq

# --------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------
N_COMP_MAIN = 8
N_COMP_SPLIT = 8
N_DIFFUSER = 15
N_NGV = 13
N_TURBINE = 12
N_VAPORIZER = 12
SHAFT_R = 4.0

AX0, AX1 = (0, 0, 0), (1, 0, 0)
HALF_PI = math.pi / 2


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def profile(start, *segs):
    """Closed (x, r) profile in the XY plane: ('l', pt) lines, ('s', [pts]) splines."""
    wp = cq.Workplane("XY").moveTo(*start)
    for kind, data in segs:
        wp = wp.lineTo(*data) if kind == "l" else wp.spline(data, includeCurrent=True)
    return wp.close()


def revolve(wp):
    return wp.revolve(360, AX0, AX1)


def poly_rev(pts):
    return revolve(cq.Workplane("XY").polyline(pts).close())


def rot(wp, deg):
    return wp.rotate(AX0, AX1, deg)


def offset_left(pts, d):
    """Offset an (x, r) polyline by d to its left (material side of the flow path)."""
    p = np.array(pts, float)
    t = np.gradient(p, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    n = np.column_stack([-t[:, 1], t[:, 0]])
    return (p + d * n).tolist()


def ell(t, x0, r0, dx, dr):
    """Quarter-ellipse meridional contour: horizontal at t=0, vertical at t=pi/2."""
    return (x0 + dx * math.sin(t), r0 + dr * (1 - math.cos(t)))


def radial_holes(solid, x, n, d, rmax, phase=0.0):
    cyls = []
    for i in range(n):
        c = cq.Solid.makeCylinder(d / 2, rmax, cq.Vector(x, 0, 0), cq.Vector(0, 1, 0))
        cyls.append(c.rotate(cq.Vector(*AX0), cq.Vector(*AX1), phase + i * 360.0 / n))
    return solid.cut(cq.Compound.makeCompound(cyls))


# --------------------------------------------------------------------------
# parts
# --------------------------------------------------------------------------
def compressor_wheel():
    th = np.linspace(0, HALF_PI, 10)[1:]
    hub = revolve(
        profile(
            (23, 4),
            ("s", [(25.5, 5.6), (28, 7.1), (30, 8)]),          # spinner / nose
            ("s", [ell(t, 30, 8, 24, 25) for t in th]),        # hub contour
            ("l", (58, 33)),
            ("l", (58, 4)),                                      # back face, bore r=4
        )
    )

    def blade(t0, thick):
        ts = np.linspace(t0, HALF_PI, 20)
        H = np.array([ell(t, 30, 8, 24, 25) for t in ts])
        d = np.array([(24 * math.cos(t), 25 * math.sin(t)) for t in ts])
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        H = H + 0.4 * np.column_stack([d[:, 1], -d[:, 0]])      # root sinks into hub
        # blade tip = shroud contour pulled 0.5 mm toward the hub, along the true normal
        S = np.array([(30 + 19 * math.sin(t), 21.5 + 11.5 * (1 - math.cos(t))) for t in ts])
        nrm = np.array([(-11.5 * math.sin(t), 19 * math.cos(t)) for t in ts])
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
        T = S - 0.5 * nrm
        pts = [tuple(p) for p in H] + [tuple(p) for p in T[::-1]]
        return cq.Workplane("XY").polyline(pts).close().extrude(thick / 2, both=True)

    wheel = hub
    for i in range(N_COMP_MAIN):
        wheel = wheel.union(rot(blade(0.0, 1.0), i * 360.0 / N_COMP_MAIN))
    for i in range(N_COMP_SPLIT):
        wheel = wheel.union(rot(blade(0.55, 0.9), (i + 0.5) * 360.0 / N_COMP_SPLIT))
    return wheel


def intake_shell():
    ph = np.linspace(0, HALF_PI, 14)
    bell = [(30 * (1 - math.cos(p)), 21.5 + 16.5 * (1 - math.sin(p))) for p in ph]
    th = np.linspace(0, HALF_PI, 14)[1:]
    shroud = [(30 + 19 * math.sin(t), 21.5 + 11.5 * (1 - math.cos(t))) for t in th]
    inner = bell + shroud
    outer = offset_left(inner, 2.0)
    wp = profile(
        inner[0],
        ("s", inner[1:]),
        ("l", (49, 53)),
        ("l", (47, 53)),
        ("l", tuple(outer[-1])),
        ("s", outer[::-1][1:]),
    )
    return revolve(wp)


def diffuser():
    plate = poly_rev([(54, 35), (56, 35), (56, 52), (54, 52)])
    L, t, a = 15.0, 1.2, math.radians(22)
    pts = []
    for u, v in [(-L / 2, -t / 2), (L / 2, -t / 2), (L / 2, t / 2), (-L / 2, t / 2)]:
        pts.append((42.5 + u * math.cos(a) - v * math.sin(a), u * math.sin(a) + v * math.cos(a)))
    vane = cq.Workplane("YZ").workplane(offset=49.0).polyline(pts).close().extrude(5.5)
    out = plate
    for i in range(N_DIFFUSER):
        out = out.union(rot(vane, i * 360.0 / N_DIFFUSER))
    return out


def casing():
    return poly_rev(
        [(47, 55), (170, 55), (182, 36.2), (196, 36.2), (232, 31.5),
         (232, 30), (196, 34.2), (180, 34.2), (168, 53), (47, 53)]
    )


def outer_liner():
    s = poly_rev([(80, 44), (80, 45), (158, 45), (172, 34), (172, 33), (158, 44)])
    s = radial_holes(s.val(), 95, 18, 4.0, 47)             # primary
    s = radial_holes(s, 135, 18, 5.0, 47, phase=10)        # dilution
    return cq.Workplane("XY").add(s)


def inner_liner():
    s = poly_rev([(80, 26), (80, 27), (158, 27), (172, 22), (172, 21), (158, 26)])
    s = radial_holes(s.val(), 95, 18, 3.5, 29)
    s = radial_holes(s, 135, 18, 4.0, 29, phase=10)
    return cq.Workplane("XY").add(s)


def dome():
    d = poly_rev([(80, 27), (81.2, 27), (81.2, 44), (80, 44)]).val()
    holes = []
    for i in range(N_VAPORIZER):
        a = math.radians(i * 360.0 / N_VAPORIZER)
        holes.append(cq.Solid.makeCylinder(3.0, 4, cq.Vector(79, 35.5 * math.cos(a), 35.5 * math.sin(a)),
                                           cq.Vector(1, 0, 0)))
    return cq.Workplane("XY").add(d.cut(cq.Compound.makeCompound(holes)))


def fuel_system():
    ring = cq.Solid.makeTorus(49, 1.5, cq.Vector(70, 0, 0), cq.Vector(1, 0, 0))
    angs = [math.radians(i * 360.0 / N_VAPORIZER) for i in range(N_VAPORIZER)]
    pts = [(35.5 * math.cos(a), 35.5 * math.sin(a)) for a in angs]
    tubes = cq.Workplane("YZ").workplane(offset=70).pushPoints(pts).circle(2.5).circle(1.5).extrude(26)
    fuel = cq.Workplane("XY").add(ring)
    for a in angs:
        stem = cq.Solid.makeCylinder(1.2, 13.5, cq.Vector(70, 35.5 * math.cos(a), 35.5 * math.sin(a)),
                                     cq.Vector(0, math.cos(a), math.sin(a)))
        fuel = fuel.union(stem)
    return fuel.union(tubes)


def ngv():
    rings = poly_rev([(172, 33), (180, 33), (180, 34), (172, 34)]).union(
        poly_rev([(172, 21), (180, 21), (180, 22), (172, 22)]))
    vane = (cq.Workplane("YZ").workplane(offset=172)
            .pushPoints([(27.5, 0)]).rect(12, 1.2).twistExtrude(8, 28))
    out = rings
    for i in range(N_NGV):
        out = out.union(rot(vane, i * 360.0 / N_NGV))
    return out


def turbine_wheel():
    hub = poly_rev([(182, 4), (182, 15), (183.5, 22), (194.5, 22), (196, 15), (196, 4)])
    blade = (cq.Workplane("YZ").workplane(offset=183)
             .pushPoints([(27.55, 0)]).rect(12.1, 1.4).twistExtrude(12, -35))
    out = hub
    for i in range(N_TURBINE):
        out = out.union(rot(blade, i * 360.0 / N_TURBINE))
    return out


def shaft():
    return cq.Workplane("YZ").workplane(offset=23).circle(SHAFT_R).extrude(173)


def bearing_tunnel():
    return poly_rev([(59, 8.2), (59, 14), (72, 14), (72, 11), (181, 11), (181, 8.2)])


def bearings():
    a = poly_rev([(61, 4), (61, 8), (66, 8), (66, 4)])
    b = poly_rev([(175, 4), (175, 8), (180, 8), (180, 4)])
    return a, b


def spider():
    strut = cq.Workplane("XY").box(3, 15.0, 3).translate((151.5, 18.5, 0))   # r = 11 .. 26, seats on tunnel + liner
    out = strut
    for i in (1, 2):
        out = out.union(rot(strut, i * 120))
    return out


def exhaust_cone():
    return revolve(
        profile((196, 0), ("l", (196, 20)),
                ("s", [(205, 17), (215, 12), (225, 7), (234, 2.5)]),
                ("l", (234, 0)))
    )


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
def build():
    b_front, b_rear = bearings()
    return {
        "intake_shell":       (intake_shell(),     (0.78, 0.80, 0.83)),
        "compressor_wheel":   (compressor_wheel(), (0.86, 0.87, 0.89)),
        "diffuser":           (diffuser(),         (0.55, 0.57, 0.60)),
        "casing":             (casing(),           (0.40, 0.42, 0.45)),
        "combustor_outer":    (outer_liner(),      (0.62, 0.50, 0.38)),
        "combustor_inner":    (inner_liner(),      (0.62, 0.50, 0.38)),
        "combustor_dome":     (dome(),             (0.58, 0.46, 0.35)),
        "fuel_system":        (fuel_system(),      (0.80, 0.62, 0.25)),
        "ngv":                (ngv(),              (0.35, 0.30, 0.28)),
        "turbine_wheel":      (turbine_wheel(),    (0.30, 0.27, 0.27)),
        "shaft":              (shaft(),            (0.90, 0.90, 0.92)),
        "bearing_tunnel":     (bearing_tunnel(),   (0.50, 0.52, 0.56)),
        "bearing_front":      (b_front,            (0.20, 0.35, 0.60)),
        "bearing_rear":       (b_rear,             (0.20, 0.35, 0.60)),
        "support_spider":     (spider(),           (0.50, 0.52, 0.56)),
        "exhaust_cone":       (exhaust_cone(),     (0.45, 0.40, 0.38)),
    }


def export(parts, step="turbojet.step", stl="turbojet.stl"):
    assy = cq.Assembly(name="KJ66_class_turbojet")
    for name, (wp, rgb) in parts.items():
        assy.add(wp, name=name, color=cq.Color(*rgb, 1.0))
    assy.save(step)
    comp = cq.Compound.makeCompound([wp.val() if len(wp.vals()) == 1
                                     else cq.Compound.makeCompound(wp.vals())
                                     for wp, _ in parts.values()])
    cq.exporters.export(cq.Workplane("XY").add(comp), stl, tolerance=0.05, angularTolerance=0.1)


if __name__ == "__main__":
    parts = build()
    for n, (wp, _) in parts.items():
        v = wp.val()
        bb = v.BoundingBox()
        print(f"{n:18s} valid={v.isValid()!s:5s} vol={v.Volume():10.1f} mm3  "
              f"x[{bb.xmin:6.1f},{bb.xmax:6.1f}] rmax={max(bb.ymax, bb.zmax):5.1f}")
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    export(parts, f"{out}/turbojet.step", f"{out}/turbojet.stl")
    print("exported")
