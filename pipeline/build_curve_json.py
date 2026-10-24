#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
build_curve_json.py
===================
.fseries  ->  compact JSON point array for the interactive 3D web viewer.

Third modality next to the PD pencil drawing (2D SVG) and the Cycles JPG (photoreal):
the bare space curve (centerline) as a rotatable chrome tube in the browser.

Bevel equivalence with the 4K renders (Werner's directive 2026-05-24:
"one bevel, equivalent to the 4k renders ... where nobody fits, nobody fits"):
the tube radius is computed with THE SAME auto_bevel_depth algorithm as
render_knot_blender.py -- min(1.5% bbox, MCA/2*safety). MCA = smallest distance
between non-adjacent centerline points. Where strands come close, the
tube stays thin enough that they do not merge ("where nobody fits, nobody fits").

Reuse instead of rebuilding: reconstruct_curve + auto_bevel_depth are ported 1:1 from
render_knot_blender.py (pure numpy, no cupy/seal import needed).
Defaults mirrored: n_samples=1024, bevel_frac=0.015, bevel_safety=0.4.

Output per knot: curves/<name>.json with
  { name, group, points:[[x,y,z],...]  (normalised to bbox_max=2, centred),
    radius (in the same normalised space), n_points, meta:{...} }
Plus curve_manifest.json (list of all knots + groups).
"""

import argparse
import site_config
import json
import math
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Ported from render_knot_blender.py (bit-identical mathematics).
# ---------------------------------------------------------------------------

def parse_fseries(path):
    """offset(3,), a(n,3), b(n,3) from a .fseries file."""
    coeffs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            coeffs.append([float(v) for v in line.split()])
    if len(coeffs[0]) == 3:
        offset = np.array(coeffs[0], dtype=np.float64)
        harmonics = coeffs[1:]
    else:
        offset = np.zeros(3)
        harmonics = coeffs
    arr = np.array(harmonics, dtype=np.float64)          # (n, 6): ax bx ay by az bz
    a = np.stack([arr[:, 0], arr[:, 2], arr[:, 4]], axis=1)   # (n,3) cos coefficients
    b = np.stack([arr[:, 1], arr[:, 3], arr[:, 5]], axis=1)   # (n,3) sin coefficients
    return offset, a, b


def reconstruct_curve(offset, a, b, n_samples=1024):
    """3D curve from Fourier coefficients, closed (s in [0, 2pi))."""
    n_terms = a.shape[0]
    s = np.linspace(0.0, 2 * np.pi, n_samples, endpoint=False)
    out = np.tile(offset.astype(np.float64), (n_samples, 1))
    for j in range(1, n_terms + 1):
        out += np.outer(np.cos(j * s), a[j - 1]) + np.outer(np.sin(j * s), b[j - 1])
    return out


# Segment-exact self-approach (numpy only, so that it runs in the base env that
# publish_site.sh uses; scipy is not installed there). No silent fallback to the old
# point measurement: if the tool is missing, the build aborts LOUDLY, because an overestimated
# clearance looks like a valid number (Disc #746/#747).
sys.path.insert(0, site_config.TOOLS_DIR)
from clearance_numpy import clearance as _segment_clearance  # noqa: E402


def auto_bevel_depth(points, max_bevel, safety=0.4):
    """Tube radius chosen so that tubes of non-adjacent strands keep a gap.
    Returns (bevel, mca, ring_skip).

    CHANGED 2026-09-04 (Werner: "you decide 3"): fixed at the root instead of tuning the
    bevel number. Until then this function measured the distance between the SAMPLED
    POINTS. That overestimates in principle (the closest point of the curve almost never
    lies on a sample point) and converges only linearly, so EVERY number here was an
    upper bound and the radius computed from it too large: for kf329 by a factor of 5.1
    at n=4096. Werner's eye caught it in the viewer, not the measurement (Disc #746/#747).
    In addition the old version built the full n^2 distance matrix and died at n=32768
    with 24 GB, at 65536 with 96 GB, so measuring more finely was not possible at all.

    Now: distance between the SEGMENTS (lower bound, quadratic convergence,
    spatial grid instead of a matrix). Measured on kf329:
        n= 4096  point 0.003264  segment 0.000638   (0.7 s)
        n= 8192  point 0.002244  segment 0.000634   (2.7 s)
        n=16384  point 0.001050  segment 0.000633   (10.6 s)
    The segment number is already stable at n=4096; the point number was still off
    by a factor of 1.66 at 16384. The radius therefore never gets larger, only more correct.
    Calibrated in clearance_numpy.calibrate() (analytic + brute force + scipy reference)."""
    n = len(points)
    diffs_seq = np.linalg.norm(np.diff(points, axis=0, append=points[:1]), axis=1)
    arc_length_total = float(diffs_seq.sum())
    ds = arc_length_total / n
    min_arc_dist = max(0.05 * arc_length_total, 2.0 * max_bevel)
    ring_skip = max(3, int(np.ceil(min_arc_dist / max(ds, 1e-9))))
    mca, _n_pairs = _segment_clearance(points, ring_skip)
    bevel = min(max_bevel, mca * 0.5 * safety)
    return bevel, mca, ring_skip


# ---------------------------------------------------------------------------

# Sculptures modelled in isometric orientation -> stand them up by +90° about X in the viewer
# (2026-06-03, Werner). ONLY these own_knots sculptures; machine knots never. Must be bit-identical to
# claude_setup_v2/tools/serve_sculpture_orientation_pd.py. y'=-z, z'=y.
SCULPTURE_ROTATIONS = {
    'trainstation', 'xyz-free-view', 'double-3', 'ping-pong', 'double-6',
    'celtic-for-playing', 'celtic-trefoil', 'knots-of-love', 'mirror',
    'sundown', 'sundowner'}
_RX90 = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])


# Hairline centerline instead of a tube (Werner 2026-06-04: "with any bevel at all
# the geometry intersects itself"). Knots with tight self-passages whose true
# clearance only becomes visible at n>=4096 (r3_phantom_diag, Disc #251 family):
# high measuring resolution, many viewer control points, sharp safety -> thin
# centerline without intersection. Survives rebuilds (the override lives in the generator).
CENTERLINE_HAIRLINE = {
    # safety 0.75 @8192: radius ~0.0005 (viewer space) = thickest intersection-free
    # hairline (true mca@8192=0.00118, the old @1024 measurement saw 0.0033 -> the tube stuck through).
    "S_45_d49_D9_chiral__m01": {"n_ctrl": 2400, "n_bevel": 8192, "safety": 0.75},
    # tw_k5: c<=60, vol 85, 256 harmonics, a lot of length packed in (Werner 2026-06-22).
    # The default of 600 n_ctrl was barely above Nyquist (2*256=512) -> coarse rendering.
    # n_ctrl/n_bevel 8192 (~16 points per highest harmonic), safety 0.75 = thin
    # intersection-free hairline for the densely packed knot.
    "tw_k5": {"n_ctrl": 65536, "n_bevel": 8192, "safety": 0.75},
    # kf329: Werner saw it in the viewer before any measurement showed it: "the tube is too
    # thick, it touches itself, but no selfintersection" (29.08.2026). Re-measured,
    # and the default n_bevel=1024 overestimates the clearance by a factor of 4:
    #     n=1024/2048  mca 0.009334  ->  radius 0.00204   (what the viewer used)
    #     n=4096       mca 0.003264  ->  radius 0.00071
    #     n=8192       mca 0.002244  ->  radius 0.00049
    # ⚠️ At 8192 the measurement has NOT converged (4096->8192 still drops by 31 %),
    #    so 0.002244 is itself still an upper bound. Measuring higher costs
    #    quadratically (199 s at 8192, ~6.4 GB at 16384), hence 8192 as for the
    #    siblings above, and safety stays at the default 0.4 instead of 0.75, because
    #    the margin is what is needed here, not the visibility.
    # 🔴 ADDENDUM 29.08. afternoon: n_bevel 8192 was NOT enough. The segment-exact
    #    measurement (`claude_setup_v2/tools/clearance_segment.py`, converges at n=16384)
    #    says the true clearance is **0.000633**; the point sampling sees instead:
    #        n= 8192  0.002244 (factor 3.5 too large)   n=65536  0.000652 (3 % too large)
    #    Recomputed in curve space: radius_view 0.00049 ⟹ curve-space radius 0.000448
    #    ⟹ two tubes need an axis distance of 0.000896, available is 0.000633
    #    ⟹ **intersection by a factor of 1.42.** (The original 0.00204: factor 5.90.)
    #    Werner had approved 0.00049 in the viewer ("the tube fits like that"); the approval
    #    was for a number I had presented to him as corrected. It was only less wrong.
    #    🔴 n_bevel 65536 is NOT feasible: this generator builds the full distance matrix
    #       (65536^2*3*8 bytes) and died with "Unable to allocate 96.0 GiB". I had inferred
    #       feasibility from my KD-tree tool (324 s there) and carried it over to the
    #       brute-force path here without checking it; memory grows
    #       quadratically, run time was the wrong metric.
    #    WORKED THROUGH, what actually works here (against the true clearance 0.000633):
    #        n= 4096  0.4 GB  safety 0.4  radius 0.000714  🔴 intersection 2.06x
    #        n= 8192  1.5 GB  safety 0.4  radius 0.000491  🔴 intersection 1.42x  ← Werner's look
    #        n=16384  6.0 GB  safety 0.4  radius 0.000230  ✅ free, margin 1.51x
    #        n=32768 24.0 GB                                🔴 exceeds memory
    #    ⟹ 16384/0.4 is the ONLY feasible value that holds. The line thereby becomes
    #       thinner once more than what Werner saw in the viewer.
    #    📌 WERNER DECIDES ON RETURN whether that is fine with him, or whether the generator
    #       is switched to segment-exact measurement (`tools/clearance_segment.py`,
    #       converges at n=16384 in 5.6 s instead of a 6 GB matrix; that would be the root fix).
    #    ✅ DONE 2026-09-04 (Werner: "you decide 3"), root fix applied:
    #    auto_bevel_depth now measures segment-exact. The true clearance is thus
    #    already stable at n_bevel=4096 (0.000638 vs 0.000633 at 16384, 0.8 % apart),
    #    in 0.7 s instead of a 6 GB matrix. The 16384 override was a crutch against the
    #    overestimating point measurement and is no longer needed; n_ctrl stays at
    #    2400 because that concerns the RENDERING, not the measurement.
    "kf329_polys08_sh0_tw-1": {"n_ctrl": 2400, "n_bevel": 4096, "safety": 0.4},
}


def process_one(fseries_path, n_ctrl, n_bevel, bevel_frac, safety):
    """-> dict for curves/<name>.json (or None on error)."""
    _name = os.path.splitext(os.path.basename(fseries_path))[0]
    if _name in CENTERLINE_HAIRLINE:
        ov = CENTERLINE_HAIRLINE[_name]
        n_ctrl = ov.get("n_ctrl", n_ctrl)
        n_bevel = ov.get("n_bevel", n_bevel)
        safety = ov.get("safety", safety)
    offset, a, b = parse_fseries(fseries_path)
    if a.shape[0] == 0:
        return None
    # Compute bevel/MCA at full render resolution (= 4K-equivalent).
    pts_hi = reconstruct_curve(offset, a, b, n_samples=n_bevel)
    ptp = np.ptp(pts_hi, axis=0)  # numpy 2.x: ndarray.ptp() removed -> np.ptp
    bbox_max = float(max(ptp[0], ptp[1], ptp[2]))
    if bbox_max <= 0:
        return None
    max_bevel = bbox_max * bevel_frac
    bevel, mca, ring_skip = auto_bevel_depth(pts_hi, max_bevel=max_bevel, safety=safety)

    # Control points for the viewer (CatmullRom interpolates smoothly in between).
    pts = reconstruct_curve(offset, a, b, n_samples=n_ctrl)

    # Centre + normalise to bbox_max = 2 (one camera setting for all knots).
    centroid = pts_hi.mean(axis=0)
    scale = 2.0 / bbox_max
    pts_n = (pts - centroid) * scale
    radius_n = bevel * scale

    name = os.path.splitext(os.path.basename(fseries_path))[0]
    # Stand isometric sculptures up by +90° about X (survives rebuilds; only own_knots,
    # machines never). 2026-06-03 Werner; see serve_sculpture_orientation_pd.py.
    rotated = name in SCULPTURE_ROTATIONS
    if rotated:
        pts_n = pts_n @ _RX90.T
    pts_list = [[round(float(v), 4) for v in p] for p in pts_n]
    meta = {
        "n_harmonics": int(a.shape[0]),
        "bbox_max": round(bbox_max, 6),
        "bevel": round(bevel, 6),
        "bevel_frac_effective": round(bevel / bbox_max, 6),
        "mca": round(mca, 6),
        # source_fseries removed (Werner 2026-07-16, variant 2): the CWD-relative path was
        # non-deterministic (phantom delta on every build) + unused by the viewer +
        # a HOW leak (internal fseries path). No active reader -> removed without replacement.
    }
    if rotated:
        meta["viewer_rotation"] = "Rx+90"
    if _name in CENTERLINE_HAIRLINE:
        meta["centerline_hairline"] = True
        meta["bevel_note"] = ("hairline centerline — tight self-passages, any normal tube "
                              "self-intersects (Werner 2026-06-04, r3_phantom_diag clearance)")
    return {
        "name": name,
        "license": "CC-BY-NC-4.0",
        "attribution": "Werner Alois Stanggassinger & Claude Code (Anthropic)",
        "points": pts_list,
        "radius": round(radius_n, 5),
        "n_points": len(pts_list),
        "closed": True,
        "meta": meta,
    }


def group_for(fseries_path, root):
    """Group label from the path (clusters/<C>, sculptures, singletons_novelty)."""
    rel = os.path.relpath(fseries_path, root)
    parts = rel.split(os.sep)
    if parts[0] == "clusters" and len(parts) >= 2:
        return parts[1]
    return parts[0]


def main():
    ap = argparse.ArgumentParser(description=".fseries -> 3D-Web-Viewer JSON")
    ap.add_argument("--root", default=site_config.WEB_EXPORT_DIR,
                    help="root folder containing clusters/ sculptures/ singletons_novelty/")
    ap.add_argument("--out", default=None,
                    help="output folder (default: <root>/curve_viewer)")
    ap.add_argument("--n-ctrl", type=int, default=600,
                    help="control points per curve in the JSON (default 600)")
    ap.add_argument("--n-bevel", type=int, default=4096,
                    help="samples for the bevel/MCA computation (= render default 1024)")
    ap.add_argument("--bevel-frac", type=float, default=0.015)
    ap.add_argument("--bevel-safety", type=float, default=0.4)
    args = ap.parse_args()

    out = args.out or os.path.join(args.root, "curve_viewer")
    curves_dir = os.path.join(out, "curves")
    os.makedirs(curves_dir, exist_ok=True)

    # Display-name layer (Werner 2026-06-21 "browse all shows singletons without an S designation"):
    # manifest label = S name/KnotInfo/artist name via naming_registry.json (keyed on distinct_knot_id;
    # for singletons == stem). Fallback = raw stem. Only the DISPLAY label; the file key stays raw (URL).
    # Consistent with build_wp_pages.disp_canonical. This second data path was never updated when the
    # naming layer was introduced (2026-06-16), which is why "browse all" stayed raw.
    _HERE = os.path.dirname(os.path.abspath(__file__))
    try:
        _NAMES = json.load(open(os.path.join(_HERE, "naming_registry.json"))).get("names", {})
    except Exception:
        _NAMES = {}
    # The display name follows the KNOT identity (distinct_knot_id), not the realisation stem:
    # realisations (e.g. st118_k2 = realisation of cm_k2) should show the S name of THEIR knot,
    # not their raw build stem. stem->distinct_knot_id from corpus_db (fix 2026-06-21, round 2).
    try:
        _cdb = json.load(open(os.path.join(_HERE, "corpus_db.json")))
        _crecs = _cdb.get("knots", _cdb)
        _crecs = list(_crecs.values()) if isinstance(_crecs, dict) else _crecs
        _STEM2DID = {r.get("stem"): r.get("distinct_knot_id") for r in _crecs if r.get("stem")}
    except Exception:
        _STEM2DID = {}
    def disp(stem):
        did = _STEM2DID.get(stem, stem)
        e = _NAMES.get(did) or _NAMES.get(stem)
        return e["display_name"] if e and e.get("display_name") else stem

    fseries = []
    for dirpath, _dirs, files in os.walk(args.root):
        if os.path.join(out, "") in os.path.join(dirpath, ""):
            continue  # do not read our own output
        for fn in files:
            if fn.endswith(".fseries"):
                fseries.append(os.path.join(dirpath, fn))
    fseries.sort()

    manifest = []
    used_keys = set()
    ok = fail = 0
    for fp in fseries:
        try:
            rec = process_one(fp, args.n_ctrl, args.n_bevel, args.bevel_frac, args.bevel_safety)
        except Exception as e:                       # noqa: BLE001
            print(f"[FAIL] {fp}: {e}", file=sys.stderr)
            fail += 1
            continue
        if rec is None:
            print(f"[SKIP] {fp}: empty/degenerate curve", file=sys.stderr)
            fail += 1
            continue
        rec["group"] = group_for(fp, args.root)
        # Unique file key: on a basename collision prepend the folder
        # (two knots with the same base name in different folders). name stays the display label.
        key = rec["name"]
        if key in used_keys:
            parent = os.path.basename(os.path.dirname(fp))
            key = f"{parent}__{rec['name']}"
        while key in used_keys:
            key += "_x"
        used_keys.add(key)
        rec["file"] = key
        with open(os.path.join(curves_dir, key + ".json"), "w") as f:
            json.dump(rec, f, separators=(",", ":"))
        # Sculptures display their ARTIST name (the stem), NOT the knot-identity name.
        # Else a sculpture realized as a known knot shows the math designation (double-3 ->
        # 12a_1019), and the 3 sculptures sharing one knot collapse to 3x "12a_1019" / 2x
        # "K8_19" in the list (Werner 2026-06-24, screenshot). disp() (knot-identity naming via
        # distinct_knot_id, the 2026-06-21 fix) stays for machines/singletons — there the build
        # stem is meaningless and the knot identity IS the right label.
        disp_name = rec["name"] if rec["group"] == "sculptures" else disp(rec["name"])
        manifest.append({
            "name": disp_name, "file": key, "group": rec["group"],
            "n_harmonics": rec["meta"]["n_harmonics"],
            "radius": rec["radius"],
        })
        ok += 1

    manifest.sort(key=lambda m: (m["group"], m["name"]))
    groups = sorted({m["group"] for m in manifest})
    with open(os.path.join(out, "curve_manifest.json"), "w") as f:
        json.dump({"count": len(manifest), "groups": groups,
                   "license": "CC-BY-NC-4.0",
                   "attribution": "Werner Alois Stanggassinger & Claude Code (Anthropic)",
                   "license_url": "https://creativecommons.org/licenses/by-nc/4.0/",
                   "knots": manifest}, f, indent=1)

    print(f"[curve-json] {ok} OK, {fail} fail -> {curves_dir}")
    print(f"[curve-json] {len(groups)} groups, manifest -> {out}/curve_manifest.json")


if __name__ == "__main__":
    main()
