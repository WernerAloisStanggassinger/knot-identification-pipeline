# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""render_knot_blender.py — headless Blender render of a knot from a .fseries.

Call (from outside Blender):
    /home/squeezy/.local/bin/blender -b \\
        /media/squeezy/kauffman_storage/freeze_manifest/render_template.blend \\
        --python /media/squeezy/kauffman_storage/freeze_manifest/pipeline_scripts/render_knot_blender.py \\
        -- --fseries /path/to/knot.fseries \\
           [--png /path/out.png] [--jpg /path/out.jpg] [--glb /path/out.glb] [--blend /path/out.blend] \\
           [--n-samples 1024] [--bevel-depth 0.003] [--padding 1.1] \\
           [--res 1920x1920] [--samples 4096]

Reads the .fseries (Fourier coefficients, 6 values per line for x/y/z
cosine/sine, optional 3-value offset first line per Discovery #69),
reconstructs the 3D curve, builds a closed poly spline with a
bevel tube, uses the stainless-steel material + studio lighting + camera
from the template .blend, fits cam.ortho_scale to the object, renders
PNG or JPG and/or exports GLB/OBJ.

Default style: kf32_radial_outward.blend (prepared by the user 2026-05).
Discovery #146 (family cleanup) — folder layout: erzeugte_knoten/<family>/<knot>/<knot>.fseries.

No external packages except Blender bpy + numpy (included in Blender).
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import bpy
import numpy as np


# --------------------------------------------------------------- fseries Parser -

def parse_fseries(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Robust against a 3-value offset first line (Discovery #69).

    Returns:
        offset: (3,) vector; (0,0,0) if there is no offset line.
        a, b:   (n_terms, 3) cosine/sine coefficients.
    """
    coeffs = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('%'):
            continue
        vals = [float(t) for t in s.split()]
        coeffs.append(vals)

    if coeffs and len(coeffs[0]) == 3:
        offset = np.array(coeffs[0], dtype=float)
        harm = coeffs[1:]
    else:
        offset = np.zeros(3)
        harm = coeffs

    if not harm:
        raise ValueError(f'.fseries contains no Fourier lines: {path}')
    if any(len(row) != 6 for row in harm):
        raise ValueError(f'.fseries lines must have 6 values (ax bx ay by az bz)')

    h = np.array(harm, dtype=float)
    a = np.stack([h[:, 0], h[:, 2], h[:, 4]], axis=1)
    b = np.stack([h[:, 1], h[:, 3], h[:, 5]], axis=1)
    return offset, a, b


def reconstruct_curve(offset: np.ndarray, a: np.ndarray, b: np.ndarray,
                      n_samples: int = 1024) -> np.ndarray:
    """3D curve from Fourier coefficients, closed (s in [0, 2pi))."""
    n_terms = a.shape[0]
    s = np.linspace(0.0, 2 * np.pi, n_samples, endpoint=False)
    out = np.tile(offset, (n_samples, 1))
    for j in range(1, n_terms + 1):
        out += np.outer(np.cos(j * s), a[j - 1]) + np.outer(np.sin(j * s), b[j - 1])
    return out


def projection_pps_arrays(points: np.ndarray, resolution_x: int,
                           padding: float) -> tuple[np.ndarray, float]:
    """Return (pps_per_sample, world_per_pixel) for the top-view render.

    pps_per_sample[i] = pixel distance between sample i and i+1 in the render
    (cyclically closed). world_per_pixel is constant from ortho_scale/res.

    Werner's diagnosis 2026-05-20 ~18:00, weak spot #2 follow-up: a global mean
    pps is not enough, because the 2D speed of the curve varies strongly locally
    at each crossing (steep Z stretches have dense 2D samples, flat
    XY stretches wide ones). Per-crossing splice needs locally measured pps."""
    pts_2d = points[:, :2]
    diffs_2d = np.linalg.norm(np.diff(pts_2d, axis=0, append=pts_2d[:1]), axis=1)
    w_img = float(max(np.ptp(points[:, 0]), np.ptp(points[:, 1])))
    ortho_scale_est = w_img * padding
    world_per_pixel = ortho_scale_est / resolution_x
    pps_arr = diffs_2d / max(world_per_pixel, 1e-12)
    return pps_arr, world_per_pixel


def estimate_pixels_per_sample(points: np.ndarray, resolution_x: int,
                                padding: float) -> float:
    """Mean pixel distance (backward-compatible wrapper around
    projection_pps_arrays)."""
    pps_arr, _ = projection_pps_arrays(points, resolution_x, padding)
    return float(pps_arr.mean())


def local_pps_at_index(pps_arr: np.ndarray, idx: int, win: int = 2) -> float:
    """Smoothed local pps value around sample index idx (cyclic). Window
    win=2 -> 5 samples (idx-2..idx+2)."""
    n = len(pps_arr)
    idxs = [(idx + d) % n for d in range(-win, win + 1)]
    return float(pps_arr[idxs].mean())


def compute_subsample_cuts(points: np.ndarray, resolution_x: int, padding: float,
                            crossings: list[dict], target_px: float
                            ) -> tuple[list[tuple[int, float, int, float]], float]:
    """Compute one cut [s_low, s_high] per crossing with s = sample_index + fraction,
    such that the 2D arc distance between s_low and s_high in the
    projection plane is exactly target_px (Werner's directive 2026-05-20:
    'tackle the sub-sample splice' for true 30 px from the projection).

    Before (integer sample): n_erase steps in pps jumps, the achieved gap
    varies between min/max depending on rounding. Now: walk the arc from the exact
    crossing point, cut at a fractional sample position.

    Returns (cuts, world_per_pixel) where cuts = list of
    (low_idx, low_frac, high_idx, high_frac).
    """
    n = len(points)
    pts_2d = points[:, :2]
    diffs_2d = np.linalg.norm(np.diff(pts_2d, axis=0, append=pts_2d[:1]), axis=1)
    # cumarc[i] = 2D arc distance from sample 0 to sample i (cumarc[n] = total length)
    cumarc = np.concatenate([[0.0], np.cumsum(diffs_2d)])
    L = float(cumarc[-1])

    w_img = float(max(np.ptp(points[:, 0]), np.ptp(points[:, 1])))
    ortho_scale_est = w_img * padding
    world_per_pixel = ortho_scale_est / resolution_x
    half_world = (target_px / 2.0) * world_per_pixel

    def locate(arc_value: float) -> tuple[int, float]:
        """Convert 2D-arc-position to (sample_idx, fraction)."""
        arc_value = arc_value % L
        idx = int(np.searchsorted(cumarc, arc_value, side='right') - 1)
        idx = max(0, min(idx, n - 1))
        seg_len = diffs_2d[idx] if idx < n else diffs_2d[-1]
        frac = (arc_value - cumarc[idx]) / max(seg_len, 1e-12)
        frac = max(0.0, min(frac, 1.0))
        return idx, float(frac)

    cuts = []
    for c in crossings:
        ui = c['under_idx']
        tu = c['t_under']
        arc_c = cumarc[ui] + tu * diffs_2d[ui]
        low_idx, low_frac = locate(arc_c - half_world)
        high_idx, high_frac = locate(arc_c + half_world)
        cuts.append((low_idx, low_frac, high_idx, high_frac))
    return cuts, world_per_pixel


def autobump_n_samples_for_pixel_gap(offset: np.ndarray, a: np.ndarray, b: np.ndarray,
                                      n_samples_initial: int, target_px: int,
                                      resolution_x: int, padding: float,
                                      cap: int = 8192
                                      ) -> tuple[int, np.ndarray, float, float]:
    """Auto-increase n_samples until max(pps_per_sample) ≤ target_px (Werner's directive
    2026-05-20, weak spots #2 + per-crossing follow-up: splice gap exactly 30 px
    IN THE PROJECTION PLANE, also at crossings where the curve locally runs fast through
    the XY plane).

    Before: pps=mean -> n_erase=1 suffices only for the average region; on steep
    XY stretches (pps_local > target) achieved_gap is far above target.
    Now: pps=max -> n_erase=1 suffices EVERYWHERE <= target, the per-crossing splice
    computes per under_idx individually.

    Cap at 8192: detect_top_view_crossings() is O(n^2) -- 16K samples would
    build 1GB arrays. 8192 -> 32M pairs ~256MB per array, survivable.

    Returns: (n_samples_eff, points, pps_max, pps_mean)
    """
    pts = reconstruct_curve(offset, a, b, n_samples=n_samples_initial)
    pps_arr, _ = projection_pps_arrays(pts, resolution_x, padding)
    pps_max = float(pps_arr.max())
    pps_mean = float(pps_arr.mean())
    if pps_max <= target_px * 1.05:
        return n_samples_initial, pts, pps_max, pps_mean
    required = int(np.ceil(n_samples_initial * pps_max / max(target_px, 1)))
    new_n = min(required, cap)
    if new_n == n_samples_initial:
        return n_samples_initial, pts, pps_max, pps_mean
    pts_new = reconstruct_curve(offset, a, b, n_samples=new_n)
    pps_arr_new, _ = projection_pps_arrays(pts_new, resolution_x, padding)
    return new_n, pts_new, float(pps_arr_new.max()), float(pps_arr_new.mean())


def detect_top_view_crossings(points: np.ndarray, ring_skip: int = 20
                              ) -> list[dict]:
    """Detect 2D-crossings of the closed curve under top-view (XY projection).

    Vectorized: pairwise segment-vs-segment intersection in XY for indices
    (i, j) with i<j and ring_skip <= j-i <= n-ring_skip (skip self-adjacency
    and wrap-around-close). Returns list of dicts with under/over indices
    determined by world-Z at the intersection point.
    """
    n = len(points)
    # Generate valid (i, j) pairs: i<j, ring_skip <= (j-i) <= n-ring_skip
    I, J = np.triu_indices(n, k=ring_skip)
    mask = (J - I) <= (n - ring_skip)
    I = I[mask]; J = J[mask]
    Ip1 = (I + 1) % n
    Jp1 = (J + 1) % n

    x1, y1 = points[I, 0],   points[I, 1]
    x2, y2 = points[Ip1, 0], points[Ip1, 1]
    x3, y3 = points[J, 0],   points[J, 1]
    x4, y4 = points[Jp1, 0], points[Jp1, 1]

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    safe = np.where(np.abs(denom) > 1e-12, denom, 1.0)
    t1 = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / safe
    t2 = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / safe

    valid = (np.abs(denom) > 1e-12) & (t1 >= 0.0) & (t1 <= 1.0) & (t2 >= 0.0) & (t2 <= 1.0)
    I = I[valid]; J = J[valid]; t1 = t1[valid]; t2 = t2[valid]
    Ip1 = (I + 1) % n; Jp1 = (J + 1) % n

    zi = points[I, 2] + t1 * (points[Ip1, 2] - points[I, 2])
    zj = points[J, 2] + t2 * (points[Jp1, 2] - points[J, 2])

    out = []
    for k in range(len(I)):
        if zi[k] > zj[k]:
            over_idx, under_idx, t_under = int(I[k]), int(J[k]), float(t2[k])
        else:
            over_idx, under_idx, t_under = int(J[k]), int(I[k]), float(t1[k])
        out.append({
            'i': int(I[k]), 'j': int(J[k]),
            'ti': float(t1[k]), 'tj': float(t2[k]),
            'zi': float(zi[k]), 'zj': float(zj[k]),
            'dz': float(abs(zi[k] - zj[k])),
            'over_idx': over_idx, 'under_idx': under_idx,
            't_under': t_under,
        })
    return out


def auto_bevel_depth(points: np.ndarray, max_bevel: float,
                     safety: float = 0.4) -> tuple[float, float, int]:
    """Choose the bevel tube radius such that non-adjacent strand tubes always
    keep a visible gap (user directive 2026-05-20: otherwise Blender shows
    intersections when centerlines pass close to each other).

    Algorithm:
      1) Compute MCA = minimum closest approach between non-adjacent
         centerline points (ring_skip adaptive: 2*max_bevel in arc length,
         excludes the tube's own neighbourhood, NOT real strand crossings).
      2) safe_bevel = MCA / 2 * safety  (safety<1.0 leaves a gap between tubes)
      3) cap at max_bevel (style limit from --bevel-frac).

    Returns: (bevel, mca, ring_skip).
    """
    n = len(points)
    # arc-length per sample (closed curve)
    diffs_seq = np.linalg.norm(np.diff(points, axis=0, append=points[:1]), axis=1)
    arc_length_total = float(diffs_seq.sum())
    ds = arc_length_total / n
    # ring_skip must satisfy two conditions:
    #   (a) arc distance > 5% of the curve length (= "another loop", not self-curvature)
    #   (b) arc distance > 2 * max_bevel (= outside the tube's own reach)
    # Without (a) MCA measures tight-bend proximity instead of real strand crossings.
    min_arc_dist = max(0.05 * arc_length_total, 2.0 * max_bevel)
    ring_skip = max(3, int(np.ceil(min_arc_dist / max(ds, 1e-9))))

    # pairwise dist matrix
    diff = points[:, None, :] - points[None, :, :]        # (n, n, 3)
    dist = np.linalg.norm(diff, axis=2)                   # (n, n)
    idx = np.arange(n)
    ring_dist = np.abs(idx[:, None] - idx[None, :])
    ring_dist = np.minimum(ring_dist, n - ring_dist)      # closed-curve ring distance
    dist[ring_dist < ring_skip] = np.inf
    mca = float(dist.min())

    safe_bevel = mca * 0.5 * safety
    bevel = min(max_bevel, safe_bevel)
    return bevel, mca, ring_skip


# --------------------------------------------------------------- Scene-Helpers -

def remove_existing_curves() -> int:
    """Delete all curve objects in the scene (the template has none, but a previous
    render run could have left some behind)."""
    n = 0
    for obj in list(bpy.data.objects):
        if obj.type == 'CURVE':
            bpy.data.objects.remove(obj, do_unlink=True)
            n += 1
    for cu in list(bpy.data.curves):
        if cu.users == 0:
            bpy.data.curves.remove(cu)
    return n


def remove_existing_lights() -> int:
    """Delete all light objects in the scene (user directive 2026-05-20: lamps
    out of the template, they produce odd artefacts). The world background
    remains as the only light source."""
    n = 0
    for obj in list(bpy.data.objects):
        if obj.type == 'LIGHT':
            bpy.data.objects.remove(obj, do_unlink=True)
            n += 1
    for lt in list(bpy.data.lights):
        if lt.users == 0:
            bpy.data.lights.remove(lt)
    return n


def set_world_white(strength: float = 2.0) -> tuple[tuple, float] | None:
    """Turn the world background to pure white + selectable strength (user directive
    2026-05-20: 'turn the world lighting further into pure white'). Since the lamps
    were removed, the world is the only light source, hence also the main lighting
    for the metal. Returns (old_colour, old_strength) for logging."""
    world = bpy.context.scene.world
    if world is None:
        return None
    world.use_nodes = True
    bg = next((n for n in world.node_tree.nodes if n.type == 'BACKGROUND'), None)
    if bg is None:
        return None
    old_color = tuple(bg.inputs['Color'].default_value)
    old_strength = float(bg.inputs['Strength'].default_value)
    bg.inputs['Color'].default_value = (1.0, 1.0, 1.0, 1.0)
    bg.inputs['Strength'].default_value = strength
    return (old_color, old_strength)


def set_standard_view_transform() -> dict:
    """Switch Blender's view transform to 'Standard' + look='None' (Werner's
    directive 2026-05-20, weak spot #1: with Filmic/AgX the white background came out only as
    R=238, not 255). Filmic/AgX apply a tone-mapping highlight roll-
    off that gamma-darkens a pure 1.0 world strength. Standard = linear sRGB
    pass-through, exposure=0, gamma=1.0 -> world 1.0 becomes a clean 255.

    Only takes effect for line art (flat black lines on pure white); for the
    stainless-steel tube the caller leaves the default Filmic curve in place, because
    metal reflections would clip otherwise.

    Returns a dict with the old values for the restore log.
    """
    sc = bpy.context.scene
    old = {
        'view_transform': sc.view_settings.view_transform,
        'look': sc.view_settings.look,
        'exposure': float(sc.view_settings.exposure),
        'gamma': float(sc.view_settings.gamma),
    }
    sc.view_settings.view_transform = 'Standard'
    sc.view_settings.look = 'None'
    sc.view_settings.exposure = 0.0
    sc.view_settings.gamma = 1.0
    return old


def get_or_create_line_art_material() -> bpy.types.Material:
    """Use the pre-baked LineArt material from render_template_lineart.blend if
    present, otherwise create it on the fly (fallback for render_template.blend
    + --line-art without the special blend)."""
    if 'LineArt' in bpy.data.materials:
        return bpy.data.materials['LineArt']
    return _create_line_art_material()


def _create_line_art_material() -> bpy.types.Material:
    """Flat-dark emission material for line-art mode — no reflection, no
    highlights, constant line look whatever the curvature. On a white
    world background, the classic knot-diagram look."""
    mat = bpy.data.materials.new(name='LineArt')
    mat.use_nodes = True
    nt = mat.node_tree
    for n in list(nt.nodes):
        nt.nodes.remove(n)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Color'].default_value = (0.0, 0.0, 0.0, 1.0)  # pure ink-black
    em.inputs['Strength'].default_value = 1.0
    nt.links.new(em.outputs['Emission'], out.inputs['Surface'])
    return mat


def find_steel_material() -> bpy.types.Material | None:
    """Find the stainless-steel material from the template (name 'Edelstahl' preferred,
    otherwise the first metallic material)."""
    if 'Edelstahl' in bpy.data.materials:
        return bpy.data.materials['Edelstahl']
    for m in bpy.data.materials:
        if m.use_nodes:
            for n in m.node_tree.nodes:
                if n.type == 'BSDF_PRINCIPLED' and n.inputs['Metallic'].default_value > 0.5:
                    return m
    return None


def _interp_point(points: np.ndarray, idx: int, frac: float) -> np.ndarray:
    """Linear interpolation between sample idx and (idx+1)%n at parameter frac."""
    n = len(points)
    return points[idx] + frac * (points[(idx + 1) % n] - points[idx])


def _merge_subsample_cuts(cuts: list[tuple[int, float, int, float]], n: int
                          ) -> list[tuple[int, float, int, float]]:
    """Sort cuts + merge overlapping compound crossings.

    Cuts span parameter intervals (li+lf, hi+hf) on the cyclic curve
    (s ∈ [0, n)). Wrap cuts (hi+hf < li+lf) are split into two intervals.
    """
    intervals = []  # (s_low, s_high, low_idx, low_frac, high_idx, high_frac)
    for li, lf, hi, hf in cuts:
        s_low = li + lf
        s_high = hi + hf
        if s_high < s_low:
            # wrap cut: split into two parts
            intervals.append((s_low, float(n), li, lf, n - 1, 1.0))
            intervals.append((0.0, s_high, 0, 0.0, hi, hf))
        else:
            intervals.append((s_low, s_high, li, lf, hi, hf))
    intervals.sort(key=lambda x: x[0])

    merged = []
    for iv in intervals:
        if not merged:
            merged.append(iv)
            continue
        p_lo, p_hi, p_li, p_lf, p_hi_idx, p_hf = merged[-1]
        s_low, s_high, li, lf, hi_idx, hf = iv
        if s_low <= p_hi:  # overlapping or touching
            if s_high > p_hi:
                merged[-1] = (p_lo, s_high, p_li, p_lf, hi_idx, hf)
        else:
            merged.append(iv)

    # wrap merge: last interval reaches the n boundary, first starts at 0
    if len(merged) >= 2:
        last = merged[-1]
        first = merged[0]
        if last[1] >= float(n) - 1e-9 and first[0] <= 1e-9:
            # fuse: low = last.low, high = first.high
            merged[0] = (last[0] - n, first[1], last[2], last[3], first[4], first[5])
            merged.pop()

    return [(iv[2], iv[3], iv[4], iv[5]) for iv in merged]


def _build_subsample_splines(cu, points: np.ndarray,
                              cuts: list[tuple[int, float, int, float]]) -> None:
    """Build splines for sub-sample cuts. Each keep region becomes a non-
    cyclic spline with interpolated start/end point + inner samples."""
    n = len(points)
    merged = _merge_subsample_cuts(cuts, n)
    K = len(merged)
    if K == 0:
        # should not happen if cuts is not empty
        return

    n_segments = 0
    for k in range(K):
        cut_k = merged[k]
        cut_next = merged[(k + 1) % K]
        # keep segment: from the end of cut_k to the start of cut_(k+1)
        s_start_idx, s_start_frac = cut_k[2], cut_k[3]
        s_end_idx, s_end_frac = cut_next[0], cut_next[1]

        # zero-length skip
        s_start = s_start_idx + s_start_frac
        s_end = s_end_idx + s_end_frac
        if s_end < s_start:
            s_end += n  # cyclic
        if (s_end - s_start) < 1e-9:
            continue

        pts_seg = [_interp_point(points, s_start_idx, s_start_frac)]
        # inner samples: walk cyclically from (s_start_idx + 1) to s_end_idx inclusive
        idx = (s_start_idx + 1) % n
        steps = 0
        if not (s_start_idx == s_end_idx and s_start_frac < s_end_frac):
            # more than one sample segment in the keep region
            while steps <= n:
                pts_seg.append(points[idx])
                if idx == s_end_idx:
                    break
                idx = (idx + 1) % n
                steps += 1
        end_pt = _interp_point(points, s_end_idx, s_end_frac)
        # double end point only if the end point differs from the last inner sample
        if not pts_seg or np.linalg.norm(np.asarray(pts_seg[-1]) - end_pt) > 1e-9:
            pts_seg.append(end_pt)

        if len(pts_seg) < 2:
            continue
        sp = cu.splines.new('POLY')
        sp.points.add(len(pts_seg) - 1)
        for i, p in enumerate(pts_seg):
            sp.points[i].co = (float(p[0]), float(p[1]), float(p[2]), 1.0)
        sp.use_cyclic_u = False
        # do NOT set order_u: POLY is linear; in Blender 5.x order_u=3 could
        # confuse the bevel pipeline's end-point behaviour (NURBS taper look).
        n_segments += 1

    print(f'[render] Sub-sample splice: {n_segments} spline pieces '
          f'(merged {K} cuts from {len(cuts)} crossings)')


def build_curve_object(points: np.ndarray, name: str, bevel_depth: float,
                       material: bpy.types.Material | None,
                       erase_mask: np.ndarray | None = None,
                       cuts: list[tuple[int, float, int, float]] | None = None,
                       line_art: bool = False
                       ) -> bpy.types.Object:
    """Build curve. Three modes:

    1) erase_mask=None, cuts=None -> single cyclic spline (no splice).
    2) cuts=list -> SUB-SAMPLE-precision splice (Werner 2026-05-20: each cut
       has an exact 2D arc length in the projection plane; end points are
       INTERPOLATED points between samples, no integer-sample quantisation
       error any more). Merge of overlapping cuts (compound crossings).
    3) erase_mask=bool array -> old integer-sample path (backward-compatible
       for --splice-radius or a splice without a --splice-pixels target).
    """
    n = len(points)
    cu = bpy.data.curves.new(name=name, type='CURVE')
    cu.dimensions = '3D'
    cu.resolution_u = 64
    cu.bevel_depth = bevel_depth
    cu.bevel_resolution = 8
    # tube mode: fill_caps close the cyclic spline harmlessly. Line art:
    # on a thin tube fill_caps make conical arrow tips at every splice end
    # ("spikes"), see kf270_v3_subsample.jpg → off in line art.
    cu.use_fill_caps = not line_art
    cu.twist_mode = 'MINIMUM'

    if cuts is not None and cuts:
        _build_subsample_splines(cu, points, cuts)
    elif erase_mask is None or not erase_mask.any():
        # standard: one cyclic spline
        sp = cu.splines.new('POLY')
        sp.points.add(n - 1)
        for i, p in enumerate(points):
            sp.points[i].co = (float(p[0]), float(p[1]), float(p[2]), 1.0)
        sp.use_cyclic_u = True
        sp.order_u = 3
    else:
        # Splice-mode: build piecewise non-cyclic splines
        keep = ~erase_mask
        # Walk cyclically — start at the first KEEP after an ERASE (= run start)
        # special case: all keep -> would be standard (handled above)
        # special case: all erased -> no spline (warning)
        if not keep.any():
            print('[render] WARN: erase_mask covers all indices -- no curve built', file=sys.stderr)
        else:
            # find the first index k with keep[k-1] False and keep[k] True (run start)
            start = 0
            for k in range(n):
                if keep[k] and not keep[(k - 1) % n]:
                    start = k
                    break
            else:
                # no ERASE found despite erase_mask.any() — should not happen
                start = 0

            runs = []
            cur = []
            for off in range(n + 1):  # +1 to close the last run
                idx = (start + off) % n
                if off < n and keep[idx]:
                    cur.append(idx)
                else:
                    if cur:
                        runs.append(cur)
                        cur = []
                    if off >= n:
                        break

            for run in runs:
                if len(run) < 2:
                    continue
                sp = cu.splines.new('POLY')
                sp.points.add(len(run) - 1)
                for i, idx in enumerate(run):
                    p = points[idx]
                    sp.points[i].co = (float(p[0]), float(p[1]), float(p[2]), 1.0)
                sp.use_cyclic_u = False
                sp.order_u = 3
            print(f'[render] Splice: {len(runs)} spline pieces built '
                  f'({n - keep.sum()}/{n} indices erased)')

    obj = bpy.data.objects.new(name, cu)
    if material is not None:
        cu.materials.append(material)
    bpy.context.collection.objects.link(obj)
    return obj


def autofit_camera(obj: bpy.types.Object, padding: float,
                   view_theta: float | None = None,
                   view_phi: float | None = None) -> dict:
    """Fits cam.ortho_scale to the object bbox.

    Default: top view, looking along -Z.
    With (theta, phi) set: spherical view from direction
        (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
    (pipeline convention from kinematic_result.json chosen.theta/phi).
    """
    import mathutils

    cam_obj = bpy.context.scene.camera
    if cam_obj is None:
        for o in bpy.data.objects:
            if o.type == 'CAMERA':
                cam_obj = o
                bpy.context.scene.camera = cam_obj
                break
    if cam_obj is None:
        raise RuntimeError('No camera in the scene (template defective)')

    # spline points in world space
    mw = obj.matrix_world
    pts = []
    for sp in obj.data.splines:
        for p in sp.points:
            v = mw @ p.co.xyz
            pts.append((v.x, v.y, v.z))
    arr = np.array(pts)
    bmin = arr.min(axis=0)
    bmax = arr.max(axis=0)

    cam_data = cam_obj.data
    cam_data.type = 'ORTHO'

    if view_theta is None or view_phi is None:
        # original top view, unchanged
        xmin, ymin, zmin = bmin
        xmax, ymax, zmax = bmax
        width = max(xmax - xmin, ymax - ymin)
        cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
        cam_z = float(zmax + max(width, 1.0))
        cam_obj.location = (float(cx), float(cy), cam_z)
        cam_obj.rotation_euler = (0.0, 0.0, 0.0)
        cam_data.ortho_scale = float(width * padding)
        cam_data.clip_start = 0.001
        cam_data.clip_end = float(max(100.0, cam_z * 10.0))
        return {
            'view_mode': 'top',
            'bbox_x': float(xmax - xmin), 'bbox_y': float(ymax - ymin), 'bbox_z': float(zmax - zmin),
            'ortho_scale': cam_data.ortho_scale, 'cam_z': cam_z,
        }

    # spherical view
    st, ct = float(np.sin(view_theta)), float(np.cos(view_theta))
    sp, cp = float(np.sin(view_phi)), float(np.cos(view_phi))
    view_dir = np.array([st*cp, st*sp, ct], dtype=float)  # FROM origin TOWARD camera

    center = 0.5 * (bmin + bmax)
    # choose a world up vector not parallel to view_dir
    if abs(float(np.dot(view_dir, [0., 0., 1.]))) < 0.95:
        up_world = np.array([0., 0., 1.])
    else:
        up_world = np.array([0., 1., 0.])

    # camera basis: the Blender camera looks along local -Z.
    # forward = look direction (from camera to object) = -view_dir
    forward = -view_dir
    right = np.cross(forward, up_world); right /= np.linalg.norm(right)
    up_cam = np.cross(right, forward);   up_cam /= np.linalg.norm(up_cam)

    # projection of all points onto the camera plane for bbox/ortho scale
    pts_c = arr - center
    proj_x = pts_c @ right
    proj_y = pts_c @ up_cam
    proj_z = pts_c @ view_dir
    width_x = float(proj_x.max() - proj_x.min())
    width_y = float(proj_y.max() - proj_y.min())
    width = max(width_x, width_y)

    cam_distance = float(proj_z.max() + max(width, 1.0))
    cam_pos = center + view_dir * cam_distance

    # rotation matrix: columns = local axes in world (X=right, Y=up_cam, Z=view_dir)
    R_world = mathutils.Matrix((
        (float(right[0]),  float(up_cam[0]),  float(view_dir[0])),
        (float(right[1]),  float(up_cam[1]),  float(view_dir[1])),
        (float(right[2]),  float(up_cam[2]),  float(view_dir[2])),
    ))

    cam_obj.location = (float(cam_pos[0]), float(cam_pos[1]), float(cam_pos[2]))
    cam_obj.rotation_euler = R_world.to_euler()
    cam_data.ortho_scale = float(width * padding)
    cam_data.clip_start = 0.001
    cam_data.clip_end = float(max(100.0, cam_distance * 10.0))

    return {
        'view_mode': 'spherical',
        'theta': float(view_theta), 'phi': float(view_phi),
        'view_dir': tuple(float(v) for v in view_dir),
        'cam_pos': tuple(float(v) for v in cam_pos),
        'bbox_x': float(bmax[0]-bmin[0]), 'bbox_y': float(bmax[1]-bmin[1]), 'bbox_z': float(bmax[2]-bmin[2]),
        'proj_width_x': width_x, 'proj_width_y': width_y,
        'ortho_scale': cam_data.ortho_scale, 'cam_distance': cam_distance,
    }


def setup_render(res: tuple[int, int], samples: int | None,
                 file_format: str = 'PNG', jpeg_quality: int = 90,
                 engine: str = 'cycles'):
    sc = bpy.context.scene
    engine_lower = engine.lower()

    # engine ID mapping (Blender 5.x): Eevee Next is the standard from 4.2 on
    if engine_lower == 'eevee':
        # try Eevee Next first (Blender 4.2+), fall back to legacy
        for eid in ('BLENDER_EEVEE_NEXT', 'BLENDER_EEVEE'):
            try:
                sc.render.engine = eid
                actual_engine = eid
                break
            except (TypeError, KeyError):
                continue
        else:
            raise RuntimeError('No Eevee engine available')
        # TAA samples (Eevee Next + legacy both)
        if samples is not None:
            try:
                sc.eevee.taa_render_samples = samples
            except AttributeError:
                pass
        sys.stderr.write(f'[render] Engine: {actual_engine} (TAA samples={samples})\n')
    elif engine_lower == 'workbench':
        sc.render.engine = 'BLENDER_WORKBENCH'
        sys.stderr.write(f'[render] Engine: BLENDER_WORKBENCH\n')
    else:
        sc.render.engine = 'CYCLES'
        if samples is not None:
            sc.cycles.samples = samples
        # GPU only for Cycles
        try:
            sc.cycles.device = 'GPU'
            prefs = bpy.context.preferences.addons['cycles'].preferences
            # OPTIX uses hardware RT cores on RTX cards -> massively faster for tube geometry.
            # Fall back to CUDA if OPTIX is not available.
            # (Texts written via sys.stderr.write are not message calls for the gate and stay as in the source.)
            for dt in ('OPTIX', 'CUDA'):
                prefs.compute_device_type = dt
                prefs.get_devices()
                gpu_devs = [d for d in prefs.devices if d.type == dt]
                if gpu_devs:
                    for d in prefs.devices:
                        d.use = (d.type == dt)
                    sys.stderr.write(f'[render] Cycles GPU: {dt} ({len(gpu_devs)} device(s))\n')
                    break
        except Exception as e:
            sys.stderr.write(f'[render] GPU-Setup-Warnung: {e}\n')

    sc.render.resolution_x = res[0]
    sc.render.resolution_y = res[1]
    sc.render.resolution_percentage = 100
    sc.render.image_settings.file_format = file_format
    if file_format == 'JPEG':
        sc.render.image_settings.color_mode = 'RGB'
        sc.render.image_settings.quality = jpeg_quality
    else:
        sc.render.image_settings.color_mode = 'RGBA'


def export_glb(obj: bpy.types.Object, out_path: Path):
    """Convert curve→mesh (duplicate) and export as GLB."""
    bpy.ops.object.select_all(action='DESELECT')
    dup = obj.copy()
    dup.data = obj.data.copy()
    dup.name = obj.name + '_glb_mesh'
    bpy.context.collection.objects.link(dup)
    bpy.context.view_layer.objects.active = dup
    dup.select_set(True)
    bpy.ops.object.convert(target='MESH')

    out_path.parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.export_scene.gltf(
        filepath=str(out_path),
        export_format='GLB',
        use_selection=True,
        export_apply=True,
        export_yup=True,
    )

    bpy.data.objects.remove(dup, do_unlink=True)


def export_obj(obj: bpy.types.Object, out_path: Path):
    """Convert curve→mesh (duplicate) and export as Wavefront OBJ.

    Universal format for 3D software (Meshlab, OpenSCAD, Mathematica,
    Geomview). Text-based, inspectable. Y-up like the GLB path.
    """
    bpy.ops.object.select_all(action='DESELECT')
    dup = obj.copy()
    dup.data = obj.data.copy()
    dup.name = obj.name + '_obj_mesh'
    bpy.context.collection.objects.link(dup)
    bpy.context.view_layer.objects.active = dup
    dup.select_set(True)
    bpy.ops.object.convert(target='MESH')

    out_path.parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.obj_export(
        filepath=str(out_path),
        export_selected_objects=True,
        forward_axis='NEGATIVE_Z',
        up_axis='Y',
        apply_modifiers=True,
    )

    bpy.data.objects.remove(dup, do_unlink=True)


# --------------------------------------------------------------- main -

def parse_args():
    argv = sys.argv
    user_args = argv[argv.index('--') + 1:] if '--' in argv else []
    p = argparse.ArgumentParser(prog='render_knot_blender')
    p.add_argument('--fseries', type=Path, required=True)
    p.add_argument('--png',     type=Path, default=None, help='PNG render path (optional)')
    p.add_argument('--jpg',     type=Path, default=None, help='JPEG render path (optional, RGB without alpha)')
    p.add_argument('--jpeg-quality', type=int, default=90, help='JPEG quality 1-100 (default 90)')
    p.add_argument('--glb',     type=Path, default=None, help='GLB export path (optional)')
    p.add_argument('--obj',     type=Path, default=None, help='Wavefront OBJ export path (optional)')
    p.add_argument('--blend',   type=Path, default=None, help='.blend save path (optional)')
    p.add_argument('--n-samples',   type=int,   default=1024)
    p.add_argument('--bevel-depth', type=float, default=None,
                   help='Absolute bevel tube radius (override, disables auto bevel).')
    p.add_argument('--bevel-frac',  type=float, default=0.015,
                   help='Maximum tube radius as a fraction of bbox_max (style cap for auto bevel; default 0.015)')
    p.add_argument('--bevel-safety', type=float, default=0.4,
                   help='Auto bevel: safety factor (tube radius = MCA/2 * safety). '
                        'Default 0.4 = tubes take 40%% of half the distance, a 20%% gap remains.')
    p.add_argument('--padding',     type=float, default=1.15,
                   help='ortho_scale = max(bbox_xy) * padding')
    p.add_argument('--res',         type=str,   default='1920x1920')
    p.add_argument('--samples',     type=int,   default=None,
                   help='Cycles samples (default: from the template = 4096)')
    p.add_argument('--world-strength', type=float, default=5.0,
                   help='World background strength (pure white); default 5.0, the template default was 3.9 with a slightly bluish colour')
    p.add_argument('--view-theta',  type=float, default=None,
                   help='Spherical-view polar angle (rad). With --view-phi: the camera is turned to this view (pipeline kinematic_result chosen.theta/phi). Without: the top view stays.')
    p.add_argument('--view-phi',    type=float, default=None,
                   help='Spherical-view azimuth angle (rad). See --view-theta.')
    p.add_argument('--splice-undercrossings', action='store_true', default=False,
                   help='Knot-diagram style: detect top-view crossings, erase the under-strand ±splice-radius samples around the crossing point. Only sensible for the top view (no --view-theta).')
    p.add_argument('--splice-radius', type=int, default=4,
                   help='[deprecated if --splice-pixels is set] splice half-width in sample indices (default 4 -> 9 samples gap per crossing).')
    p.add_argument('--splice-pixels', type=int, default=None,
                   help='Splice gap in pixels (Werner directive 2026-05-20: 30 px). '
                        'If set: overrides --splice-radius, samples are computed from the render resolution + 2D arc density.')
    p.add_argument('--splice-n-cap', type=int, default=8192,
                   help='Maximum n_samples in the auto bump for --splice-pixels '
                        '(default 8192). Higher = closer gap approximation, but '
                        'detect_top_view_crossings is O(n^2) -- 16K would build GB arrays.')
    p.add_argument('--splice-ring-skip', type=int, default=20,
                   help='Lower bound of the ring distance for crossing detection (default 20 -> ignore local self-curvature).')
    p.add_argument('--splice-dz-factor', type=float, default=3.0,
                   help='Splice only if |Z_over - Z_under| < dz_factor * bevel_depth. '
                        'Default 3.0 -> only tight crossings where tubes would touch. '
                        'Widely spread crossings are left to the natural Z-buffer render logic.')
    p.add_argument('--line-art', action='store_true', default=False,
                   help='Line-art mode: thin tube + flat-dark material. Knot-diagram look '
                        'for knots where the tube render becomes catastrophic (too dense crossings). '
                        'Forces a splice for ALL crossings (no dz filter), overrides --bevel-depth '
                        'and the stainless-steel material.')
    p.add_argument('--line-bevel-frac', type=float, default=0.0015,
                   help='Line-art tube thickness as a fraction of bbox_max (default 0.0015 ≈ 4-5px at 4K).')
    p.add_argument('--engine', choices=['cycles', 'eevee', 'workbench'], default=None,
                   help='Render engine. Default: cycles (3D tube) or eevee (--line-art). '
                        'Werner 2026-05-20: line-art special cases use Eevee for speed.')
    return p.parse_args(user_args)


def main():
    args = parse_args()
    if not args.fseries.exists():
        print(f'ERROR: fseries not found: {args.fseries}', file=sys.stderr)
        sys.exit(1)
    if args.png is None and args.jpg is None and args.glb is None and args.obj is None and args.blend is None:
        print('ERROR: give at least --png, --jpg, --glb, --obj or --blend', file=sys.stderr)
        sys.exit(1)

    w, h = [int(x) for x in args.res.lower().split('x')]
    name = args.fseries.stem

    print(f'[render] Knot: {name}')
    print(f'[render] Template objects: {[o.name for o in bpy.data.objects]}')

    # The template comes already loaded (blender -b template.blend --python ...).
    # If someone starts the script without a template argument, there are no
    # lights/camera/material -> abort.
    if not any(o.type == 'CAMERA' for o in bpy.data.objects):
        print('ERROR: no camera in the scene -- please load render_template.blend', file=sys.stderr)
        sys.exit(2)

    remove_existing_curves()
    n_lights = remove_existing_lights()
    if n_lights:
        print(f'[render] {n_lights} light(s) removed (user directive 2026-05-20)')
    world_change = set_world_white(strength=args.world_strength)
    if world_change is not None:
        old_color, old_strength = world_change
        print(f'[render] World->pure white strength={args.world_strength:.2f} '
              f'(before: color={tuple(round(c,3) for c in old_color)} strength={old_strength:.2f})')
    if args.line_art:
        old_vt = set_standard_view_transform()
        print(f'[render] View-Transform: Standard look=None exposure=0 gamma=1.0 '
              f'(before: {old_vt["view_transform"]}/{old_vt["look"]} exp={old_vt["exposure"]:.2f} '
              f'gamma={old_vt["gamma"]:.2f}) — white-background fix, weak spot #1')
    steel = find_steel_material()
    if steel is None:
        print('WARNING: no stainless-steel material found in the template', file=sys.stderr)

    offset, a, b = parse_fseries(args.fseries)
    points = reconstruct_curve(offset, a, b, n_samples=args.n_samples)
    n_samples_eff = args.n_samples

    # Auto-bump n_samples for line art + a splice_pixels target (weak spot #2).
    # Makes the splice gap pixel-exact in the projection plane -- computed locally
    # per crossing, therefore the MAX pps must lie below target, not the mean.
    # At crossings where the curve runs fast through XY locally, a
    # mean-based bump has achieved_gap > target (Werner's diagnosis 2026-05-20).
    if args.line_art and args.splice_pixels is not None:
        n_new, points_new, pps_max, pps_mean = autobump_n_samples_for_pixel_gap(
            offset, a, b, args.n_samples, args.splice_pixels, w, args.padding,
            cap=args.splice_n_cap)
        if n_new != args.n_samples:
            print(f'[render] Auto-Bump n_samples: {args.n_samples} -> {n_new} '
                  f'(splice_pixels target={args.splice_pixels}, pps_max={pps_max:.2f}, pps_mean={pps_mean:.2f})')
            if pps_max > args.splice_pixels * 1.05:
                print(f'[render] WARN: pps_max={pps_max:.2f} > target={args.splice_pixels} '
                      f'even after the bump (cap={args.splice_n_cap}). At the steepest '
                      f'crossings the achieved gap will be ~{pps_max:.1f} px instead of {args.splice_pixels} px. '
                      f'Raise --splice-n-cap for a closer projection approximation.', file=sys.stderr)
            points = points_new
            n_samples_eff = n_new

    bbox_max = float(max(np.ptp(points[:, 0]), np.ptp(points[:, 1]), np.ptp(points[:, 2])))
    if args.line_art:
        bevel = bbox_max * args.line_bevel_frac
        print(f'[render] bbox_max={bbox_max:.3f}  bevel_depth={bevel:.5f}  '
              f'(line-art: line_bevel_frac={args.line_bevel_frac})')
    elif args.bevel_depth is not None:
        bevel = args.bevel_depth
        print(f'[render] bbox_max={bbox_max:.3f}  bevel_depth={bevel:.4f}  (explicit override)')
    else:
        max_bevel = bbox_max * args.bevel_frac
        bevel, mca, ring_skip = auto_bevel_depth(points, max_bevel=max_bevel,
                                                  safety=args.bevel_safety)
        cap_hit = "CAP" if bevel >= max_bevel - 1e-12 else "MCA-limited"
        print(f'[render] bbox_max={bbox_max:.3f}  bevel_depth={bevel:.4f}  '
              f'(auto: {cap_hit}, MCA={mca:.4f}, cap={max_bevel:.4f}, safety={args.bevel_safety}, ring_skip={ring_skip})')

    erase_mask = None
    subsample_cuts = None
    # Line art: forces a splice for ALL crossings (knot-diagram convention,
    # no dz filter — a line has no Z depth, so every under-strand must be broken).
    do_splice = args.splice_undercrossings or args.line_art
    if do_splice:
        if args.view_theta is not None or args.view_phi is not None:
            print('[render] WARN: splice + spherical view not combined -- splice is designed for the top view only', file=sys.stderr)
        crossings = detect_top_view_crossings(points, ring_skip=args.splice_ring_skip)
        if args.line_art:
            tight = crossings  # ALL crossings spliced
            dz_threshold = None
        else:
            dz_threshold = args.splice_dz_factor * bevel
            tight = [c for c in crossings if c['dz'] < dz_threshold]
        if crossings:
            dz_min = min(c['dz'] for c in crossings)
            dz_med = float(np.median([c['dz'] for c in crossings]))
            print(f'[render] Crossings detected: {len(crossings)} total, '
                  f'dz_min={dz_min:.5f} dz_median={dz_med:.5f}')

        # splice computation:
        # - with --splice-pixels: sub-sample cuts (exact 2D arc distance in the
        #   projection plane, end points interpolated between samples)
        # - otherwise: integer-sample erase (old path, --splice-radius)
        if args.splice_pixels is not None:
            if tight:
                subsample_cuts, world_per_pixel = compute_subsample_cuts(
                    points, w, args.padding, tight, args.splice_pixels)
                # diagnostic: each cut spans exactly target_px in 2D arc -- but
                # after _merge_subsample_cuts overlapping compound crossings are
                # merged, so a merged cut width can be larger.
                # Compute here the planned single-cut widths (target_px everywhere)
                # + the effective merge statistics after _merge_subsample_cuts.
                merged_for_stat = _merge_subsample_cuts(subsample_cuts, len(points))
                pts_2d = points[:, :2]
                diffs_2d = np.linalg.norm(np.diff(pts_2d, axis=0, append=pts_2d[:1]), axis=1)
                merged_widths_world = []
                for li, lf, hi, hf in merged_for_stat:
                    # Total 2D arc inside the merged interval
                    s_lo = li + lf
                    s_hi = hi + hf
                    if s_hi < s_lo:
                        s_hi += len(points)
                    # walk integer samples
                    w_total = 0.0
                    cur_idx = li
                    cur_frac = lf
                    target_param = s_hi
                    cur_param = s_lo
                    while cur_param < target_param - 1e-12:
                        next_param = float(cur_idx + 1)
                        step = min(next_param, target_param) - cur_param
                        w_total += step * diffs_2d[cur_idx % len(points)]
                        cur_param += step
                        cur_idx = int(cur_param) % len(points)
                    merged_widths_world.append(w_total)
                merged_widths_px = [w / max(world_per_pixel, 1e-12) for w in merged_widths_world]
                pps_arr, _ = projection_pps_arrays(points, w, args.padding)
                pps_mean = float(pps_arr.mean())
                if merged_widths_px:
                    mw_min = min(merged_widths_px); mw_max = max(merged_widths_px)
                    mw_med = float(np.median(merged_widths_px))
                    print(f'[render] Sub-sample splice: target={args.splice_pixels}px exactly per '
                          f'crossing, {len(tight)} crossings -> {len(merged_for_stat)} merged cuts '
                          f'(compound crossings merged). Effective cut widths '
                          f'min/med/max = {mw_min:.1f}/{mw_med:.1f}/{mw_max:.1f} px '
                          f'(n_samples={n_samples_eff}, pps_mean={pps_mean:.2f})')
            else:
                print(f'[render] Splice pixel-mode: no crossing cluster to splice')
        else:
            # backward-compat: splice_radius gives total span = 2*radius + 2
            n_erase_global = 2 * args.splice_radius + 2
            print(f'[render] Splice sample-mode: radius={args.splice_radius} -> {n_erase_global} samples erased per crossing')
            if tight:
                erase_mask = np.zeros(len(points), dtype=bool)
                for c in tight:
                    ui = c['under_idx']
                    half_lo = (n_erase_global - 1) // 2
                    half_hi = n_erase_global // 2
                    for d in range(-half_lo, half_hi + 1):
                        erase_mask[(ui + d) % len(points)] = True
                mode = "ALL (line-art)" if args.line_art else f"tight (dz < {dz_threshold:.4f} = {args.splice_dz_factor}x bevel)"
                print(f'[render] Splice applied: {len(tight)} {mode}')

        if not tight and args.splice_pixels is None:
            print(f'[render] Splice: 0 crossings to splice -- the Z-buffer handles all over/under')

    mat_to_use = get_or_create_line_art_material() if args.line_art else steel
    if args.line_art:
        print('[render] Material: LineArt (flat near-black emission, no reflection)')
    obj = build_curve_object(points, name=name, bevel_depth=bevel,
                             material=mat_to_use, erase_mask=erase_mask,
                             cuts=subsample_cuts, line_art=args.line_art)

    info = autofit_camera(obj, padding=args.padding,
                          view_theta=args.view_theta, view_phi=args.view_phi)
    print(f'[render] bbox(x,y,z)=({info["bbox_x"]:.3f}, {info["bbox_y"]:.3f}, {info["bbox_z"]:.3f})')
    if info['view_mode'] == 'top':
        print(f'[render] view=top  ortho_scale={info["ortho_scale"]:.3f}  cam_z={info["cam_z"]:.3f}')
    else:
        print(f'[render] view=spherical theta={info["theta"]:.4f} phi={info["phi"]:.4f}  '
              f'view_dir=({info["view_dir"][0]:.3f}, {info["view_dir"][1]:.3f}, {info["view_dir"][2]:.3f})')
        print(f'[render] ortho_scale={info["ortho_scale"]:.3f}  cam_distance={info["cam_distance"]:.3f}  '
              f'proj_bbox=({info["proj_width_x"]:.3f}, {info["proj_width_y"]:.3f})')

    # engine choice: explicit --engine, otherwise eevee for line art, otherwise cycles
    engine = args.engine if args.engine else ('eevee' if args.line_art else 'cycles')
    print(f'[render] Engine: {engine}')

    if args.png is not None:
        setup_render((w, h), args.samples, file_format='PNG', engine=engine)
        args.png.parent.mkdir(parents=True, exist_ok=True)
        bpy.context.scene.render.filepath = str(args.png)
        eff_samples = bpy.context.scene.cycles.samples
        print(f'[render] PNG -> {args.png}  ({w}x{h}, samples={eff_samples})')
        bpy.ops.render.render(write_still=True)
        print(f'[render] PNG OK: {args.png}')

    if args.jpg is not None:
        setup_render((w, h), args.samples, file_format='JPEG', jpeg_quality=args.jpeg_quality, engine=engine)
        args.jpg.parent.mkdir(parents=True, exist_ok=True)
        bpy.context.scene.render.filepath = str(args.jpg)
        eff_samples = bpy.context.scene.cycles.samples
        print(f'[render] JPG -> {args.jpg}  ({w}x{h}, samples={eff_samples}, q={args.jpeg_quality})')
        bpy.ops.render.render(write_still=True)
        print(f'[render] JPG OK: {args.jpg}')

    if args.glb is not None:
        print(f'[glb]    GLB -> {args.glb}')
        export_glb(obj, args.glb)
        print(f'[glb]    GLB OK: {args.glb}')

    if args.obj is not None:
        print(f'[obj]    OBJ -> {args.obj}')
        export_obj(obj, args.obj)
        print(f'[obj]    OBJ OK: {args.obj}')

    if args.blend is not None:
        args.blend.parent.mkdir(parents=True, exist_ok=True)
        bpy.ops.wm.save_as_mainfile(filepath=str(args.blend.resolve()), compress=True)
        print(f'[blend]  Saved: {args.blend}')


if __name__ == '__main__':
    main()
