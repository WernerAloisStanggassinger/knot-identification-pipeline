#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# ============================================================================
# canonical_flow.py — 6D machine BUILD STAGE 2, SONO v3 (2026-06-06)
#
# History (all three attempts bit-diagnosed, see NOW.md 2026-06-05):
#   v1 force flow (spring+repulsion, fixed dt):   200ms/step + dt creep.
#   v2 force flow (adaptive): stiff system, shape crumples instead of relaxing.
#   SONO v1: trefoil collapsed to a CIRCLE — excl window grew with D
#     without bound, from excl>=n/2 no checkable pairs -> "ok" always true ->
#     D grew freely, smoothing+resample pulled the knot through itself.
#   CORE LESSON: preservation of isotopy depends EXCLUSIVELY on the active
#   overlap wall — smoothing without a wall = tunnelling.
#
# SONO v2 — fix recipe (NOW.md safety wrap-up 2026-06-05 ~20:00):
#  (a) excl cap:  excl = min(int(1.5*D/h)+2, n//8). Math receipt: the cap
#      only binds from D > 1/12 (L=1) <=> RL=2/D < 24 — BELOW the trefoil
#      lower bound 32.74. For real knots it NEVER binds, it is the pure
#      anti-tunnel wall against the circle collapse.
#  (b) D growth ONLY if real pairs were checked (pairs_seen > 0);
#      an empty candidate set due to excl saturation != success.
#  (c) Convergence = RL plateau (best_D not improved for STALL_STOP
#      iterations) + K_tot plausibility: Fary-Milnor says the total curvature
#      of every real knot is >= 4pi <=> K_tot >= 2.0. K_tot < 1.85 during
#      the flow => tunnelled => ABORT + CIRCLE_COLLAPSE (unless the
#      record is a declared unknot, singletons_unknots — there the
#      circle is the CORRECT normal form).
#  (d) Calibration anchor: ideal trefoil RL ~ 32.74 (literature norm L/r; D is
#      the tube DIAMETER, so RL = 2/D at L=1), |Wr| ~ 3.42.
#  (e) Couple n to RL: after the plateau enforce n >= 4*RL (upsample +
#      continued flow); in addition start-density probe on the 16n polyline
#      (densify doctrine): bead spacing h <= t0/2, otherwise the wall is porous.
#
# Wall-integrity hardenings (from the core lesson, beyond the recipe):
#  (f) Step cap in the overlap fix: displacement per bead per pass
#      <= 0.30*D — large accumulated pushes can tunnel by themselves.
#  (g) SONO reaction to unresolvable overlap: loosen D slightly (*0.998)
#      instead of leaving beads permanently inside the wall.
#
# SONO v3 — window-crumple fix (calibration FAIL #4, NOW.md 2026-06-05 ~20:30):
#   v2 held the GLOBAL wall, but the flow found the nearest blind
#   spot: the curve crumpled LOCALLY into the excl window (ring distance
#   3..excl is invisible to the wall; kf22: thickness 1.6e-5 at D=0.007,
#   K_tot 248 = zigzag on bead scale). Mechanism: point pushes up to 0.3*D
#   = at D >> h several bead spacings per pass -> kinks; smoothing 0.02
#   too weak against it. The system ALWAYS finds the wall's blind spot.
#  (1) CURVATURE CONTROL = second half of the wall (the missing SONO piece,
#      Pieranski has it for a reason): per iteration enforce bending radius >= D/2
#      via a kink-angle cap per bead theta_max = 2*asin(h/D).
#      Math receipt crumple tightness: curve with r >= D/2 has in the window
#      (arc length <= 1.5*D <= pi*r) chord length 2r*sin(s/2r) monotonically
#      growing, at s=1.5D: 0.997*D — the curve CANNOT in the window come back
#      to a distance << D. Wall gap-free: outside overlap pushes,
#      inside curvature. D growth only if BOTH walls stand; otherwise
#      loosen (rope too thick to bend = same reaction as (g)).
#  (2) Spread pushes over neighbouring beads (Gaussian falloff 5 taps,
#      sigma 1.2) instead of a point yank — kinks do not arise in the first place.
#  (3) Step cap additionally in units of h: <= min(0.3*D, 1.0*h).
#  (4) Honest metrology: v2 wall_ratio used thickness(ring>2) — at
#      D >> h the min distance of ring-3 beads is ~3h, the measure could NEVER
#      reach 1.0. v3 measures thickness_wall at the final excl window
#      (wall_ratio ~1.0 = the wall really held) + r_min (min bending radius,
#      r_min_over_halfD ~>1.0 = the curvature wall held).
#   Calibration expectation: trefoils |Wr|~3.42 / RL~32.7 / K_tot~2.1-2.5,
#   unknot -> circle without alarm.
#
# SONO v3.2 — double-chokehold fix (calibration run 2, 2026-06-06, 7/7 ok but
#   6/7 JAM on intermediate plateaus RL 78-170; only celtic hit the anchors.
#   Even the UNKNOT did not find the circle (RL=77.7, Wr=-5.24 instead of 6.3/0).
#   Both findings MEASURED on the kf22 patience run — STALL_STOP=2000 delivered
#   bit-the-same plateau, hence a jam, not impatience):
#  (A) PINCH-CANCEL: 16236 residual overlap pairs on the 'converged' curve,
#      1755/1768 contact beads in >=2 pairs = pressure from both sides, pushes
#      cancel when summed, overlap structurally NEVER becomes zero ->
#      zero-tolerance gate blocked D forever. Fix: tolerance gate — D
#      grows as long as residual DEPTH <= TOL_DEPTH*D (measured pinch: 2.4%).
#      Pieranski tolerates residual overlap for exactly this reason.
#  (B) KERNEL INJECTION: the sigma=1.2 kernel injected kinks per wall pass of
#      ~|w''|*push/h ~ 0.7 rad >> theta_max -> curv_resid oscillated above
#      the 2% threshold (measured 0.029/0.173/0.030) -> the curvature gate
#      additionally choked D. Fix: sigma=5, +-10 taps -> injection
#      ~1/sigma^2 = 0.04 rad ~ theta_max(D/h=32). Downsampling keeps
#      D/h <= 32, so the one sigma value fits everywhere.
#  (C) POLISH PHASE before the fingerprint: measurement happened in the middle
#      of the push-iron cycle (celtic r/D2=0.04, K=9.6 instead of ~2.2) -> at
#      fixed D=0.98*plateau relax wall+curvature, THEN measure.
#  (D) Downsample halving: the condition n_t < n//2 required 4x jumps,
#      halving never triggered (kf59/kf51 dragged n=4096 along for nothing).
#
# SONO v3.3 — writhe-transport fix (calibration run 3, 2026-06-06: FIRST
#   NORMAL-FORM COLLAPSE celtic/kf85/kf59 onto ONE point RL 33.7-33.9 /
#   Wr -3.49..-3.51; remaining stuck ones = exactly the writhe-distant embeddings
#   kf20/22 (+1.0 instead of +3.42), kf51 (+0.14), unknot (-5.96 instead of 0) —
#   Loop flipping is the slowest mode. Pressure-hold experiment
#   (kf22 plateau, D fixed 1.02*plateau, 3000 iter, tmp_durable/
#   diag_v3_pressure_hold.*) delivered both lessons bit-exact):
#  (E) THE SHRINK REFLEX (g) WAS A PRESSURE-RELIEF VALVE: every failed iteration
#      loosened D at once -> the sustained pressure needed for loop flipping
#      was reset every time (Sisyphus). Fix: slow valve —
#      shrink only after VALVE_EVERY failed iterations in a row, in between
#      HOLD D (sustained pressure).
#  (F) PRESSURE ONLY PRODUCTIVE ON THE SMOOTH MANIFOLD: in the hold
#      experiment curv_resid stood at ~0.6 (8 ironing passes in a row
#      hopeless against 3000-pair pressure fields) and the 'creep' was
#      pure kink shuffling (mean_move stagnated/fell, K~34). Fix:
#      INTERLEAVING — after each push pass a mini-iron (4 passes),
#      after the pass block iron-until-clean (32 passes, O(n), free);
#      the pre-block iron is dropped (redundant).
#
# SONO v3.4 — DEEP-CYCLE / basin hopping (calibration run 4 + anneal experiments
#   2026-06-06: v3.3 anchor trio fixed (33.6-33.7/-3.46..-3.51), unknot
#   92.9->19.4, BUT gauss class hard clamped at RL 95-96 (Wr +0.93
#   instead of +3.42) — shallow anneal (0.6*D, 300 iter) bounced off 5/5.
#   DIFFUSION MATH: writhe sits in global modes k=1..10, smoothing
#   relaxes mode k ~ S*iters*(2pi k/n)^2 -> unreachable at n=1024,
#   at n=256 + S=0.1 x5 substeps ~O(1) in 3000 iter. DEEP-ANNEAL
#   EXPERIMENT (tmp_durable/diag_v3_deep_anneal.*): kf22 clamp after
#   deep rounds at Wr=+3.341/ACN=4.03/K=2.81 (ideal shape BEFORE the
#   reinflate!), end RL=38.1/Wr=+4.03 — BASIN HOPPED, topology-
#   safe (wall stayed active, clearance ~ D_plateau >> 0.3*D_plateau).
#   INTEGRATION: on plateau stall -> deep_round() (coarse n~4/best_D, deflate
#   0.3*best_D, DEEP_ROUND_ITERS strong smoothing + wall) -> reinflate in the
#   main flow (D=0.3*best_D, stall reset), max DEEP_CYCLES attempts;
#   snapshot restore: if a cycle makes things worse, the best D counts.
# ============================================================================
import argparse
import site_config
import json
import math
import os
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from sixd_scan import parse_fseries, detect_period_pi, eval_derivatives

ROOT = site_config.ROOT
CORPUS_DB = ROOT + "/freeze_manifest/pipeline_scripts/corpus_db.json"
OUT_DIR = ROOT + "/sixd_machine/out2"
CURVE_DIR = ROOT + "/sixd_machine/canonical_curves"

N_BEADS = 512
N_BEADS_DENSE = 1024      # H > 100
N_MAX = 4096              # memory cap (O(n^2) pair matrix)
MAX_ITERS = 20000
MAX_PASSES = 8            # overlap-resolution passes per iteration
CHECK_EVERY = 500         # log interval
STALL_STOP = 400          # best_D not improved for this long -> plateau
GROW_CLEAN = 1.01         # no overlap found -> inflate fast
GROW_WORK = 1.002         # overlap resolved -> inflate carefully
SHRINK = 0.998            # (g) unresolvable -> loosen the rope
SMOOTH_S = 0.02
DISP_CAP = 0.30           # (f) max displacement per bead/pass in units of D
KTOT_CIRCLE = 1.85        # (c) Fary-Milnor alarm (real knot: >= 2.0)
RL_TREFOIL = 32.74        # (d) literature anchor, |Wr| ~ 3.42
CURV_RESID_OK = 0.02      # v3(1): max fraction of beads above kink cap for D-grow
TOL_DEPTH = 0.05          # v3.2(A): pinch tolerance — residual overlap depth <= 5%*D
POLISH_ITERS = 300        # v3.2(C): max relax iterations before the fingerprint
VALVE_EVERY = 10          # v3.3(E): shrink only after this many fails in a row
CURV_MINI = 4             # v3.3(F): mini-iron passes after each push pass
CURV_CLEAN = 32           # v3.3(F): iron-until-clean after the pass block
# v3.4c ADAPTIVE CYCLES (run-5b finding: 5/6 trefoils at 33.6-33.7/|3.47|,
# but kf20 stayed in the 37.8 basin and the unknot stopped at RL 16.2 — 2 fixed
# cycles are sometimes too few; basin-hopping outcome depends on the state
# before rounding): keep hopping until DEEP_FAILS_STOP cycles in a row brought
# < 2% best_D improvement (cap DEEP_CYCLES_MAX); rounding budget
# escalates per cycle (x cycle number = deterministic diversity).
DEEP_CYCLES_MAX = 8       # 2026-06-12 REVERTED 3->8: the cap=3 was WRONG —
                          # finding (out2 vs out2_robust): STIFF records (trefoils!)
                          # need >3 deep cycles to collapse; cap=3 left them
                          # un-collapsed (K3_1 RL 33.6->138, 11n_133 54->193). My
                          # "stiff stops early anyway via DEEP_FAILS_STOP" was empirically
                          # wrong. cap=8 (v3.5b) collapses correctly (Disc #259).
                          # Floppy waste remains open -> adaptive floppy skip
                          # (RL stays high -> non-collapsing -> abort) = TODO.
DEEP_FAILS_STOP = 3       # v3.4c: stop after non-improvers in a row (2->3 likewise)
# v3.6 (2026-06-12, Werner's go): adaptive FLOPPY SKIP. Problem: floppy records
# all grind through DEEP_CYCLES_MAX cycles because they drift marginally (>DEEP_IMPROVE),
# without ever collapsing — time wasters, and according to the 6D scope finding unreliable on plateaus
# anyway (W4). FIX: abort deep cycles when a record does NOT collapse
# = RL does not fall meaningfully from the deep start. TRAJECTORY-based (NOT absolute
# RL — otherwise one cuts off mid-collapse trefoils, cap=3 lesson 2026-06-12). Stiff
# (RL falls steeply) keeps its cycles; floppy (RL plateau) stops early.
FLOPPY_SKIP_ENABLED = False  # 2026-06-12: implemented, but UNVERIFIED (smoke
                          # timed out, no confirmation that trefoil still collapses).
                          # cap=3 lesson: unverified flow change NOT live.
                          # Next session: smoke (trefoil→33.6 + floppy skips) → only
                          # then to True. Until then: calibration runs the verified cap=8.
FLOPPY_PROBE_CYCLES = 4   # give AT LEAST this many deep cycles (chance to collapse)
FLOPPY_KEEP_FRAC = 0.85   # RL must have fallen >15% from the deep start = "collapsing"
FLOPPY_RL_FLOOR = 90      # target only high-RL time wasters (stiff ones <90 stop adaptively anyway)
DEEP_IMPROVE = 1.02       # v3.4c: a cycle counts as improvement from +2% best_D
DEEP_ROUND_ITERS = 3000   # v3.4e (Disc #260): CONSTANT — linear (x cycle
                          # number) was 56% idle; replay: const. budget =
                          # 35% cost with 329/329 deliveries kept
DEEP_SMOOTH = 0.1         # v3.4: strong smoothing during deep rounding
DEEP_SUBSTEPS = 5         # v3.4: smoothing substeps per rounding iteration
DEEP_DEFLATE = 0.3        # v3.4: deflation factor on best_D
# --- v3.5 (Disc #262): TUNNEL GUARD + DENSE DEEP WALL ----------------
# Finding Disc #261: deep_round at n_deep=max(256, 2*RL) has h ~ Dd/0.6 =
# porous wall -> strand passage for all normal forms with RL >~ 100
# (v14 tunnel, RL~118-127 collapse class). The objective (D growth)
# REWARDS the passage -> silent CONVERGED on the wrong knot. Third
# instance of the meta-lesson "the flow finds the blind spot of the wall".
TUNNEL_GUARD = 1.25       # v3.5: best_D jump over cycle mark after
                          # deep cycle > +25% (= RL crash > -20%) = strand
                          # passage; legitimate cycle gains are a few %
                          # (measured: kf237 plateaus 136->130.7->128.3 = -4/-2%)
# v3.5b (Werner stopped 85-min runs, Disc #260 pattern repeated!): the
# v3.5 "dense deep wall" (h<=Dd/8 -> n_deep up to 4096) was a budget
# buster WITHOUT delivery capability — diffusion scales n^-2, the dense
# round CANNOT hop, it only burns O(n^2)x3000 per cycle.
# MEASUREMENT (re-run-35 abort): 3/35 after 85 min, the landed ones show
# restore point == end result (kf40: RL 125.0). Consequence: after the
# first tunnel finding the deep machine is SHUT DOWN (deep_dead)
# — the restore-plateau state is the honest result (upper bound).
# Basin hopping and large normal forms are mutually exclusive in the SONO frame
# : coarse=porous (tunnels), fine=impotent (does not hop).
# v3.2(B): Gaussian spread kernel for pushes, sigma=5, +-10 taps, PEAK-normalised
# (contact bead gets the full push). Kink injection ~ 1/sigma^2 = 0.04 rad ~
# theta_max at D/h=32 (downsample cap) — before, sigma=1.2: 0.7 rad,
# 10x above cap, fed the curv_resid choker.
KER_HALF = 10
_k = np.exp(-0.5 * (np.arange(-KER_HALF, KER_HALF + 1) / 5.0) ** 2)
PUSH_KERNEL = _k / _k.max()
# v3 speed: 13 half-neighbourhood offsets for the cell grid (do not count
# pairs twice; same cell handled separately)
_HALF_OFFSETS = [(0, 0, 1), (0, 1, -1), (0, 1, 0), (0, 1, 1),
                 (1, -1, -1), (1, -1, 0), (1, -1, 1),
                 (1, 0, -1), (1, 0, 0), (1, 0, 1),
                 (1, 1, -1), (1, 1, 0), (1, 1, 1)]

_RING_CACHE = {}          # ring-distance matrix per n (constant, NEVER per step)
_SEEN_CACHE = {}          # number of checkable pairs per (n, excl)


def ring_dist(n):
    if n not in _RING_CACHE:
        idx = np.arange(n)
        b = np.abs(idx[:, None] - idx[None, :])
        _RING_CACHE[n] = np.minimum(b, n - b).astype(np.int32)
    return _RING_CACHE[n]


def pairs_seen(n, excl):
    """(b): how many pairs does the wall really see at this excl?
    Closed formula (ring): pairs with ring distance <= excl are n*excl
    (excl < n/2), the rest of n*(n-1)/2 is checkable — no O(n^2) cache any more."""
    if excl >= n // 2:
        return 0
    return n * (n - 1) // 2 - n * excl


def overlap_pairs(x, n, excl, D):
    """v3 speed: candidate pairs (dist < D, ring distance > excl) via a
    cell grid with edge length D (Pieranski's 'skipping'): same cell
    + 13 half-neighbourhood offsets, O(n*k) instead of the O(n^2) full matrix
    (calibration-run-1 finding 2026-06-06: 4096^2 pair_dist = ~470 MB temporaries
    and ~2 s PER PASS -> 4-8 h per record; the grid makes ms of it)."""
    c = np.floor(x / D).astype(np.int64)
    c -= c.min(0)
    dims = c.max(0) + 2
    key = (c[:, 0] * dims[1] + c[:, 1]) * dims[2] + c[:, 2]
    order = np.argsort(key, kind="stable")
    ks = key[order]
    cand_i, cand_j = [], []

    # same cell: all pairs within each group (groups are small)
    starts = np.flatnonzero(np.concatenate([[True], ks[1:] != ks[:-1]]))
    ends = np.concatenate([starts[1:], [len(ks)]])
    sizes = ends - starts
    for size in np.unique(sizes):
        if size < 2:
            continue
        gs = starts[sizes == size]
        loc_i, loc_j = np.triu_indices(int(size), 1)
        cand_i.append(order[(gs[:, None] + loc_i[None, :]).ravel()])
        cand_j.append(order[(gs[:, None] + loc_j[None, :]).ravel()])

    # neighbour cells: per offset, hit ranges via searchsorted, expanded
    for ox, oy, oz in _HALF_OFFSETS:
        nk = key + (ox * dims[1] + oy) * dims[2] + oz
        lo = np.searchsorted(ks, nk, side="left")
        hi = np.searchsorted(ks, nk, side="right")
        cnt = hi - lo
        src = np.flatnonzero(cnt)
        if len(src) == 0:
            continue
        m = cnt[src]
        rep_src = np.repeat(src, m)
        off = np.arange(int(m.sum())) - np.repeat(np.cumsum(m) - m, m)
        rep_dst = order[np.repeat(lo[src], m) + off]
        cand_i.append(rep_src)
        cand_j.append(rep_dst)

    if not cand_i:
        return (np.empty(0, np.int64), np.empty(0, np.int64),
                np.empty(0, np.float64))
    ii = np.concatenate(cand_i)
    jj = np.concatenate(cand_j)
    rd = np.abs(ii - jj)
    rd = np.minimum(rd, n - rd)
    keep = rd > excl
    ii, jj = ii[keep], jj[keep]
    d = x[ii] - x[jj]
    dist = np.sqrt((d * d).sum(1))
    keep = dist < D
    return ii[keep], jj[keep], dist[keep]


def pair_dist(x, n, excl):
    pdm = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    pdm[ring_dist(n) <= excl] = np.inf
    return pdm


def min_nonneighbor_dist(x, excl, chunk=2048):
    """min non-neighbour distance, chunked (also for the 16n density probe)."""
    n = len(x)
    best = np.inf
    idx = np.arange(n)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        d = np.linalg.norm(x[a:b, None, :] - x[None, :, :], axis=2)
        band = np.abs(idx[a:b, None] - idx[None, :])
        band = np.minimum(band, n - band)
        d[band <= excl] = np.inf
        best = min(best, float(d.min()))
    return best


# ------------------------------------------------------------- polyline toolkit
def resample_uniform(x, n):
    """closed polyline -> n points uniform in arc length, L=1."""
    seg = np.linalg.norm(np.roll(x, -1, 0) - x, axis=1)
    L = seg.sum()
    s = np.concatenate([[0.0], np.cumsum(seg)])
    xs = np.vstack([x, x[:1]])
    su = np.linspace(0.0, L, n, endpoint=False)
    out = np.empty((n, 3))
    for k in range(3):
        out[:, k] = np.interp(su, s, xs[:, k])
    return out / L


def quick_ktot(x):
    e = np.roll(x, -1, 0) - x
    T = e / np.linalg.norm(e, axis=1)[:, None]
    cosang = np.clip(np.einsum('ij,ij->i', T, np.roll(T, 1, 0)), -1, 1)
    return float(np.arccos(cosang).sum() / (2 * math.pi))


def gauss_wr_acn(mid, T, ls, chunk=1024):
    """Gauss double sum Wr/ACN, chunked (otherwise n=4096 gives ~GB temporaries)."""
    n = len(mid)
    wr = acn = 0.0
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        d = mid[a:b, None, :] - mid[None, :, :]
        dist = np.linalg.norm(d, axis=2)
        rows = np.arange(a, b)
        dist[rows - a, rows] = np.inf
        TxT = np.cross(T[a:b, None, :], T[None, :, :])
        num = np.einsum('ijk,ijk->ij', TxT, d)
        w = ls[a:b, None] * ls[None, :]
        integ = num / dist ** 3 * w
        wr += float(integ.sum())
        acn += float(np.abs(integ).sum())
    return wr / (4 * math.pi), acn / (4 * math.pi)


def polyline_fingerprint(x):
    n = len(x)
    e = np.roll(x, -1, 0) - x                    # segments
    ls = np.linalg.norm(e, axis=1)
    L = ls.sum()
    T = e / ls[:, None]
    mid = x + 0.5 * e

    wr, acn = gauss_wr_acn(mid, T, ls)

    # turning angle (discrete curvature) + dihedral (discrete torsion)
    cosang = np.clip(np.einsum('ij,ij->i', T, np.roll(T, 1, 0)), -1, 1)
    turn = np.arccos(cosang)
    K_tot = turn.sum() / (2 * math.pi)
    B = np.cross(np.roll(T, 1, 0), T)
    nB = np.linalg.norm(B, axis=1)
    okB = nB > 1e-12
    Bu = np.zeros_like(B)
    Bu[okB] = B[okB] / nB[okB, None]
    b1, b2 = Bu, np.roll(Bu, -1, 0)
    cosd = np.clip(np.einsum('ij,ij->i', b1, b2), -1, 1)
    sind = np.einsum('ij,ij->i', np.cross(b1, b2), T)
    dih = np.arctan2(sind, cosd)
    T_sig = dih.sum() / (2 * math.pi)

    thickness = min_nonneighbor_dist(x, 2)

    com = x.mean(0)
    rgyr = float(np.sqrt(((x - com) ** 2).sum(1).mean()))

    qs = [0.05, 0.25, 0.50, 0.75, 0.95]
    turn_n = turn * n / (2 * math.pi)            # dimensionless kink
    kq = [float(np.quantile(turn_n, q)) for q in qs]

    return dict(Wr=float(wr), ACN=float(acn), K_tot=float(K_tot),
                T_sig=float(T_sig), SL=float(wr + T_sig),
                thickness=float(thickness),
                inv_thickness=float(1.0 / thickness),
                R_gyr=rgyr, kappa_q=kq, L=float(L))


def turn_angles(x):
    """Kink angle at each bead (between segment i-1 and segment i)."""
    e = np.roll(x, -1, 0) - x
    T = e / np.maximum(np.linalg.norm(e, axis=1), 1e-12)[:, None]
    cosang = np.clip(np.einsum('ij,ij->i', T, np.roll(T, 1, 0)), -1, 1)
    return np.arccos(cosang)


def curvature_passes(x, n, h, D, max_passes=MAX_PASSES):
    """v3(1): enforce bending radius >= D/2 — kink-angle cap per bead.
    Discrete: r = h/(2*sin(theta/2)) >= D/2  <=>  theta <= 2*asin(h/D).
    Pull violating beads towards the neighbour midpoint (at the midpoint
    collinear = angle 0; lam = 1 - theta_max/theta to first order, the passes
    fetch the rest). Step cap like v3(3): <= min(0.3*D, h) — the
    curvature correction must not itself tunnel.
    Returns: x, residual fraction of beads above the cap (wall-2 residual)."""
    th_max = 2.0 * math.asin(min(1.0, h / D))
    cap = min(DISP_CAP * D, 1.0 * h)
    for _pass in range(max_passes):
        turn = turn_angles(x)
        viol = turn > th_max
        if not viol.any():
            return x, 0.0
        lam = np.zeros(n)
        lam[viol] = np.clip(1.0 - th_max / turn[viol], 0.0, 1.0)
        mid = 0.5 * (np.roll(x, 1, 0) + np.roll(x, -1, 0))
        move = lam[:, None] * (mid - x)
        dn = np.linalg.norm(move, axis=1)
        scale = np.minimum(1.0, cap / np.maximum(dn, 1e-12))
        x = x + move * scale[:, None]
    return x, float((turn_angles(x) > th_max * 1.05).mean())


def spread_disp(disp):
    """v3(2)/v3.2(B): point yank -> Gaussian-spread displacement over
    +-KER_HALF neighbours (circular convolution, peak-normalised) — kinks
    do not arise in the first place."""
    out = PUSH_KERNEL[KER_HALF] * disp
    for k in range(1, KER_HALF + 1):
        out += PUSH_KERNEL[KER_HALF - k] * np.roll(disp, -k, 0)
        out += PUSH_KERNEL[KER_HALF + k] * np.roll(disp, k, 0)
    return out


def min_bend_radius(x):
    """v3(4): minimal discrete bending radius r = h_loc/(2*sin(theta/2))."""
    e = np.roll(x, -1, 0) - x
    ls = np.linalg.norm(e, axis=1)
    th = turn_angles(x)
    hloc = 0.5 * (ls + np.roll(ls, 1))
    r = hloc / np.maximum(2.0 * np.sin(th / 2.0), 1e-12)
    return float(r.min())


def excl_cap(n, circle):
    """v3.4d: excl cap. n//8 was anti-tunnel insurance against the
    v1 circle collapse, but binds ONLY in the RL<24 regime (no real knot
    has RL<24, trefoil minimum 32.74) and blocked the unknot end state
    at RL~16 (circle chord over 1/8 arc = 0.122 -> wall read its own
    curvature neighbourhood as overlap; measured 16.0/16.2 bit-exact).
    allow_circle records get n//2-1: the wall sees exactly the
    opposite pairs -> D_max=2R -> RL=2pi. The anti-tunnel guarantee is carried by
    (b) seen>0 + the curvature wall."""
    return (n // 2 - 1) if circle else (n // 8)


def wall_block(x, n, h, D, circle=False):
    """One complete wall block (v3.3 mechanics): push<->mini-iron
    interleave + iron-until-clean. Returns x, ok_geom (no overlap
    or residual depth <= TOL_DEPTH*D, AND curvature clean), worked."""
    excl = min(int(1.5 * D / h) + 2, excl_cap(n, circle))
    cap = min(DISP_CAP * D, 1.0 * h)
    ok, worked = False, False
    for _pass in range(MAX_PASSES + 1):
        ii, jj, dist = overlap_pairs(x, n, excl, D)
        if len(ii) == 0:
            ok = True
            break
        if _pass == MAX_PASSES:
            ok = float(D - dist.min()) <= TOL_DEPTH * D  # v3.2(A)
            break
        worked = True
        d = x[ii] - x[jj]
        u = d / np.maximum(dist, 1e-12)[:, None]
        push = 0.55 * (D - dist)[:, None] * u
        disp = np.zeros_like(x)
        np.add.at(disp, ii, push)
        np.add.at(disp, jj, -push)
        disp = spread_disp(disp)                         # v3(2)
        dn = np.linalg.norm(disp, axis=1)
        scale = np.minimum(1.0, cap / np.maximum(dn, 1e-12))
        x = x + disp * scale[:, None]
        x, _ = curvature_passes(x, n, h, D, max_passes=CURV_MINI)  # v3.3(F)
    x, curv_resid = curvature_passes(x, n, h, D, max_passes=CURV_CLEAN)
    return x, (ok and curv_resid <= CURV_RESID_OK), worked, curv_resid


def deep_round(x, n, best_D, log_lines, iters=DEEP_ROUND_ITERS,
               circle=False):
    """v3.4: deep rounding = basin-hopping core. Coarse (diffusion of global
    modes ~n^-2 -> writhe transport becomes reachable), deflate, strong
    smoothing. CAUTION (Disc #261): the old assumption "wall stays active =
    topology-safe" is WRONG for RL >~ 100 (h ~ Dd/0.6 = porous) —
    hence v3.5b: tunnel finding shuts down the deep machine for the record
    (deep_dead in the flow), NO dense rounds (budget buster
    without hop capability). Return value resampled to the main n."""
    Dd = DEEP_DEFLATE * best_D
    n_deep = max(256, 1 << math.ceil(math.log2(max(4.0 / best_D, 256.0))))
    n_deep = min(n, n_deep)
    h = 1.0 / n_deep
    xr = resample_uniform(x, n_deep)
    for _it in range(iters):
        for _s in range(DEEP_SUBSTEPS):
            xr = (1 - DEEP_SMOOTH) * xr + DEEP_SMOOTH * 0.5 * (
                np.roll(xr, 1, 0) + np.roll(xr, -1, 0))
        xr = resample_uniform(xr, n_deep)
        xr, _, _, _ = wall_block(xr, n_deep, h, Dd, circle)
    log_lines.append(f"deep_round: n_deep={n_deep} Dd={Dd:.5f} "
                     f"({iters} iter) done")
    return resample_uniform(xr, n)


# ------------------------------------------------------------------- flow
def flow(x, n, allow_circle, log_lines):
    """SONO v3 (Pieranski style): let the tube diameter D grow, push overlap
    pairs apart geometrically. Fixes (a)-(g) + v3(1)-(3), see header.
    Tube = TWO walls: outside overlap pushes, inside bending radius >= D/2.
    Convergence = RL=2/D plateau + Fary-Milnor guard."""
    h = 1.0 / n
    x = resample_uniform(x, n)
    t_start = min_nonneighbor_dist(x, 2)
    D = max(min(0.9 * t_start, 2.5 * h), 1.2 * h)
    best_D, stall = D, 0
    cur_best = D     # v3.4 fix: resettable watermark — stall measures progress
    #                  of the CURRENT phase; measured against the all-time best_D
    #                  every reinflate after a deflate (needs ~600 iter
    #                  at GROW_WORK to get back above best_D) would be cut off by the
    #                  400 stall window (run-5 finding: celtic
    #                  999+2x400, kf22 1865+2x400 steps bit-exact)
    x_best, n_best = x.copy(), n                         # v3.4 snapshot
    curv_resid = 0.0
    fails = 0                                            # v3.3(E) valve counter
    deep_used = 0                                        # v3.4 cycle counter
    deep_spent = 0                                       # v3.4e Disc #260: reporting
    deep_fails = 0                                       # v3.4c consecutive duds
    best_at_cycle = 0.0                                  # v3.4c comparison mark
    rl_deep_start = 0.0                                  # v3.6 RL at the 1st deep cycle (floppy skip)
    x_cycle, n_cycle, D_cycle = None, None, 0.0          # v3.5 cycle snapshot
    deep_dead = False                                    # v3.5b deep shut down
    tunnel_hits = 0                                      # v3.5 guard counter
    status = "RUNNING"
    it = 0
    while it < MAX_ITERS:
        it += 1
        # smoothing (= shrinks) + renormalisation L=1 (= relative D growth)
        x = (1 - SMOOTH_S) * x + SMOOTH_S * 0.5 * (np.roll(x, 1, 0) +
                                                   np.roll(x, -1, 0))
        x = resample_uniform(x, n)
        seen = pairs_seen(n, min(int(1.5 * D / h) + 2,
                                 excl_cap(n, allow_circle)))      # (b)
        x, ok_geom, worked, curv_resid = wall_block(
            x, n, h, D, allow_circle)                             # v3.3
        ok = ok_geom and seen > 0      # (b) empty candidate set != success
        if ok:
            D *= GROW_WORK if worked else GROW_CLEAN
            fails = 0
        else:
            fails += 1
            if fails >= VALVE_EVERY:                     # v3.3(E) slow
                D = max(D * SHRINK, 1.2 * h)             # valve instead of reflex
                fails = 0
            # otherwise: HOLD D = sustained pressure for loop flipping
        if D > cur_best * 1.0005:                        # phase watermark
            cur_best, stall = D, 0
        else:
            stall += 1
        if D > best_D * 1.0005:                          # all-time best mark
            best_D = D
            x_best, n_best = x.copy(), n                 # v3.4 snapshot
        if it % 50 == 0 and not allow_circle and quick_ktot(x) < KTOT_CIRCLE:
            # v3.5b: after a deep cycle the passage is healable — restore
            # to the pre-cycle best state + SHUT DOWN deep machine
            # (heals the bit-stable collapse alarms, kf40 type; dense
            # rounds were budget busters without delivery capability)
            if deep_used > 0 and x_cycle is not None and not deep_dead:
                tunnel_hits += 1
                log_lines.append(
                    f"it {it}: K_tot < {KTOT_CIRCLE} after deep cycle "
                    f"-> tunnelled (guard {tunnel_hits}) -> restore "
                    f"x_cycle (n={n_cycle}, RL={2.0/D_cycle:.1f}), "
                    f"deep machine shut down")
                x, n = x_cycle.copy(), n_cycle
                h = 1.0 / n
                best_D = D_cycle
                x_best, n_best = x_cycle.copy(), n_cycle
                D = D_cycle
                cur_best = D
                stall, fails = 0, 0
                deep_dead = True
                continue
            status = "CIRCLE_COLLAPSE"                   # (c) circle alarm
            log_lines.append(f"it {it}: K_tot < {KTOT_CIRCLE} -> "
                             f"Fary-Milnor violated -> tunnelled, ABORT")
            break
        if it % 50 == 0 and n > N_BEADS:
            # (e) in reverse: couple n to D — the start n protects the THIN
            # raw curve (h <= t0/2); once the tube is inflated, the
            # wall at h <= D/16 holds just as tight, and the cell grid
            # stays scale-free (~const. beads/cell). Only real
            # halvings (no flicker); moves stay <= h << D.
            need_n = max(16.0 / D, float(N_BEADS))   # clamp BEFORE the shift
            n_t = min(N_MAX, 1 << math.ceil(math.log2(need_n)))
            # v3.2(D): n_t < n allows halvings (n_t is pow2, so
            # n_t < n <=> n_t <= n/2); the old n//2 condition required
            # 4x jumps and never triggered at halving
            if n_t < n:
                log_lines.append(f"it {it}: D={D:.5f} -> Downsample "
                                 f"n {n} -> {n_t} (h <= D/16 kept)")
                x = resample_uniform(x, n_t)
                n = n_t
                h = 1.0 / n
                continue
        if it % CHECK_EVERY == 0:
            log_lines.append(f"it {it}: D={D:.5f} RL={2.0/D:.1f} "
                             f"stall={stall} "
                             f"curv_resid={curv_resid:.3f} n={n}")
        if stall >= STALL_STOP:
            # v3.5b TUNNEL GUARD: the objective (D growth) rewards
            # strand passages — a deep cycle after which best_D jumps by
            # more than TUNNEL_GUARD (RL crash) has changed
            # the knot. Restore + shut down the deep machine instead of
            # a silent CONVERGED on the wrong knot (kf237->trefoil type).
            # Only as long as deep lives — after deep_dead best_D only grows
            # via the main loop (legitimate), no false positive.
            if (deep_used > 0 and not deep_dead and best_at_cycle > 0.0
                    and best_D > TUNNEL_GUARD * best_at_cycle):
                tunnel_hits += 1
                log_lines.append(
                    f"it {it}: TUNNEL GUARD {tunnel_hits}: best_D "
                    f"{best_D:.5f} > {TUNNEL_GUARD}x cycle mark "
                    f"{best_at_cycle:.5f} (RL {2.0/best_D:.1f} << "
                    f"{2.0/best_at_cycle:.1f}) -> strand passage -> "
                    f"restore x_cycle (n={n_cycle}), deep shut down")
                x, n = x_cycle.copy(), n_cycle
                h = 1.0 / n
                best_D = D_cycle
                x_best, n_best = x_cycle.copy(), n_cycle
                D = D_cycle
                cur_best = D
                stall, fails = 0, 0
                deep_dead = True
                continue
            RL = 2.0 / best_D
            need = 4.0 * RL                              # (e) n >= 4*RL
            if n < need and n < N_MAX:
                n_new = min(N_MAX, 1 << math.ceil(math.log2(need)))
                log_lines.append(f"it {it}: Plateau RL={RL:.1f}, "
                                 f"n {n} -> {n_new} (recipe e: n>=4*RL)")
                x = resample_uniform(x, n_new)
                n = n_new
                h = 1.0 / n
                stall = 0
                continue
            # v3.4c ADAPTIVE: evaluate the previous cycle, then possibly the next one
            if deep_used > 0:
                if best_D > DEEP_IMPROVE * best_at_cycle:
                    deep_fails = 0                       # cycle delivered
                else:
                    deep_fails += 1                      # dud
            # v3.6 FLOPPY SKIP: abort non-collapsing records (RL has not
            # fallen meaningfully since the deep start). Accept the plateau (W4: floppy
            # anyway unreliable). Trajectory-based -> does NOT cut collapsing ones.
            rl_now = 2.0 / best_D
            if (FLOPPY_SKIP_ENABLED and deep_used >= FLOPPY_PROBE_CYCLES
                    and rl_now > FLOPPY_RL_FLOOR
                    and rl_now > FLOPPY_KEEP_FRAC * rl_deep_start
                    and not allow_circle):
                status = "CONVERGED"
                log_lines.append(
                    f"it {it}: FLOPPY-SKIP after {deep_used} deep cycles — RL={rl_now:.1f} "
                    f"> {FLOPPY_KEEP_FRAC}*deep start({rl_deep_start:.1f}) = not "
                    f"collapsing, stopped (plateau accepted)")
                break
            if (deep_used < DEEP_CYCLES_MAX and deep_fails < DEEP_FAILS_STOP
                    and not deep_dead):
                # v3.4 DEEP-CYCLE: plateau != proof of ideality — basin hopping
                # (deep rounding + reinflate); the snapshot keeps the best
                best_at_cycle = best_D
                # v3.5: pre-cycle snapshot — rollback anchor of the guard
                x_cycle, n_cycle, D_cycle = x_best.copy(), n_best, best_D
                if deep_used == 0:
                    rl_deep_start = 2.0 / best_D         # v3.6 RL at the 1st deep cycle
                deep_used += 1
                log_lines.append(f"it {it}: Plateau RL={RL:.1f} -> "
                                 f"DEEP-CYCLE {deep_used}/{DEEP_CYCLES_MAX} "
                                 f"(fails={deep_fails}, "
                                 f"budget={DEEP_ROUND_ITERS})")
                deep_spent += DEEP_ROUND_ITERS           # v3.4e constant
                x = deep_round(x, n, best_D, log_lines,
                               DEEP_ROUND_ITERS,
                               allow_circle)
                D = DEEP_DEFLATE * best_D
                cur_best = D           # set the watermark to the phase start
                stall, fails = 0, 0
                continue
            status = "CONVERGED"
            log_lines.append(f"converged @ it {it}: D plateau {best_D:.5f}, "
                             f"RL={RL:.1f}, n={n}, curv_resid={curv_resid:.3f}")
            break
    if status == "RUNNING":
        status = "MAX_ITERS"
        log_lines.append(f"MAX_ITERS {MAX_ITERS} reached, D={D:.5f} "
                         f"RL={2.0/D:.1f}")
    if D < 0.97 * best_D:
        # v3.4 snapshot restore: end state worse than the best plateau
        # (e.g. a failed last deep cycle) -> the best D counts
        log_lines.append(f"snapshot restore: final D={D:.5f} < best "
                         f"{best_D:.5f} -> x_best (n={n_best})")
        x, n = x_best, n_best
    return (x, n, best_D, it, status, curv_resid, deep_spent,
            tunnel_hits, deep_dead)


def polish(x, n, best_D, log_lines, circle=False):
    """v3.2(C): snapshot hygiene — the flow ends in the middle of the push-iron
    cycle (wall passes leave kinks that only the NEXT iteration
    would iron out). At fixed D = 0.98*plateau, without smoothing/growth,
    relax curvature+wall, ONLY THEN measure."""
    h = 1.0 / n
    D = 0.98 * best_D
    excl = min(int(1.5 * D / h) + 2, excl_cap(n, circle))
    cap = min(DISP_CAP * D, 1.0 * h)
    curv_resid, n_pairs = 1.0, -1
    for it in range(POLISH_ITERS):
        x, curv_resid = curvature_passes(x, n, h, D, max_passes=CURV_CLEAN)
        ii, jj, dist = overlap_pairs(x, n, excl, D)
        n_pairs = len(ii)
        deep = n_pairs and float(D - dist.min()) > TOL_DEPTH * D
        if curv_resid == 0.0 and not deep:
            log_lines.append(f"polish: clean after {it} iter "
                             f"(D=0.98*plateau, remaining pairs {n_pairs} "
                             f"all <= {TOL_DEPTH:.0%} depth)")
            return x
        if n_pairs:
            d = x[ii] - x[jj]
            u = d / np.maximum(dist, 1e-12)[:, None]
            push = 0.55 * (D - dist)[:, None] * u
            disp = np.zeros_like(x)
            np.add.at(disp, ii, push)
            np.add.at(disp, jj, -push)
            disp = spread_disp(disp)
            dn = np.linalg.norm(disp, axis=1)
            scale = np.minimum(1.0, cap / np.maximum(dn, 1e-12))
            x = x + disp * scale[:, None]
    log_lines.append(f"polish: remainder after {POLISH_ITERS} iter "
                     f"(curv_resid={curv_resid:.3f}, pairs={n_pairs})")
    return x


# ------------------------------------------------------------------ record
def run_record(job):
    stem, fpath = job
    try:
        offset, A, B = parse_fseries(fpath)
        Hn = A.shape[0]
        period = math.pi if detect_period_pi(A, B) else 2 * math.pi
        allow_circle = "singletons_unknots" in fpath
        n0 = N_BEADS_DENSE if Hn > 100 else N_BEADS
        # start curve: dense Fourier evaluation (densify doctrine)
        t = np.linspace(0.0, period, 16 * n0, endpoint=False)
        r, _, _, _ = eval_derivatives(offset, A, B, t)
        r16 = resample_uniform(r, 16 * n0)
        # (e) start-density probe on the 16n polyline: h <= t0/2 needed
        t0 = min_nonneighbor_dist(r16, 32)
        n_req = int(math.ceil(2.0 / max(t0, 1e-9)))
        n = max(n0, 1 << math.ceil(math.log2(max(n_req, 2))))
        n = min(n, N_MAX)
        wall_porous = (1.0 / n) > (t0 / 2.0)
        x = resample_uniform(r16, n)
        fp0 = polyline_fingerprint(x)
        log_lines = []
        if n != n0:
            log_lines.append(f"start density t0={t0:.5f} -> n={n} "
                             f"(porous={wall_porous})")
        (x, n_fin, best_D, steps, status, curv_resid, deep_spent,
         tunnel_hits, deep_dead) = flow(x, n, allow_circle, log_lines)
        x = resample_uniform(x, n_fin)
        x = polish(x, n_fin, best_D, log_lines, allow_circle)  # v3.2(C)
        fp = polyline_fingerprint(x)
        fp["ropelength"] = float(2.0 / best_D)           # (d) literature norm L/r
        fp["D_final"] = float(best_D)
        # residual overlap energy (wall residual after the final resample)
        excl_fin = min(int(1.5 * best_D * n_fin) + 2,
                       excl_cap(n_fin, allow_circle))
        pdm = pair_dist(x, n_fin, excl_fin)
        viol = best_D - pdm[pdm < best_D]
        fp["E_rep"] = float((viol ** 2).sum())
        # (c) end guard (catches a collapse between two 50-iteration checks)
        if not allow_circle and fp["K_tot"] < KTOT_CIRCLE \
                and status == "CONVERGED":
            status = "CIRCLE_COLLAPSE"
        conv = status == "CONVERGED" or (allow_circle and status != "MAX_ITERS")
        degenerate = fp["thickness"] < 0.25 / n_fin      # pull-tight guard
        # v3(4) honest wall metrology: outside measured at the excl window
        # (~1.0 = the overlap wall held), inside via bending radius (>=1.0 = the
        # curvature wall held). v2 wall_ratio (thickness ring>2 / D) was
        # structurally blind at D >> h (min distance ring-3 ~ 3h, never ~D).
        thickness_wall = min_nonneighbor_dist(x, excl_fin)
        wall_ratio = thickness_wall / best_D
        r_min = min_bend_radius(x)
        fp["thickness_wall"] = float(thickness_wall)
        fp["r_min"] = float(r_min)
        fp["r_min_over_halfD"] = float(r_min / (0.5 * best_D))
        out = {
            "stem": stem, "fseries": fpath,
            "n_beads": n, "n_beads_final": int(n_fin),
            "n_harmonics": int(Hn),
            "period": "pi" if period == math.pi else "2pi",
            "steps": steps, "deep_iters": int(deep_spent),
            "status": status, "converged": bool(conv),
            "tunnel_hits": int(tunnel_hits),             # v3.5 guard
            "deep_dead": bool(deep_dead),                # v3.5b deep shut down
            "allow_circle": bool(allow_circle),
            "wall_porous_at_start": bool(wall_porous),
            "t0_dense": float(t0), "wall_ratio": float(wall_ratio),
            "curv_resid_final": float(curv_resid),
            "degenerate_thickness": bool(degenerate),
            "fingerprint_start": fp0, "fingerprint": fp,
            "flow_log": log_lines,
            "machine": "canonical_flow SONO-v3.5b 2026-06-07 (Disc #262: tunnel guard + deep_dead after tunnel; #260: constant deep budget)",
        }
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(os.path.join(OUT_DIR, stem + "_sixd2.json"), "w") as f:
            json.dump(out, f, indent=1)
        np.save(os.path.join(CURVE_DIR, stem + "_canonical.npy"), x)
        flag = status + " "
        if degenerate:
            flag = "DEGEN " + flag
        if tunnel_hits:
            flag += f"tun={tunnel_hits} "
        return (f"{flag}{stem} steps={steps} deep={deep_spent} RL={fp['ropelength']:.1f} "
                f"Wr={fp['Wr']:+.3f} ACN={fp['ACN']:.2f} "
                f"K={fp['K_tot']:.3f} wall={wall_ratio:.2f} "
                f"r/D2={fp['r_min_over_halfD']:.2f} n={n_fin}")
    except Exception:
        return f"FAIL {stem}\n{traceback.format_exc()}"


def collect_jobs(only=None):
    db = json.load(open(CORPUS_DB))
    jobs = []
    pool = list(db["knots"]) + list(db.get("validation_knots", []))
    for rec in pool:
        stem = rec["stem"]
        if only and stem not in only:
            continue
        folder = os.path.dirname(rec["pflicht_path"])
        fpath = os.path.join(folder, stem + ".fseries")
        if not os.path.isfile(fpath):
            cands = [f for f in os.listdir(folder)
                     if f.endswith(".fseries")
                     and not f.startswith("_sculpture_")] \
                if os.path.isdir(folder) else []
            if len(cands) == 1:
                fpath = os.path.join(folder, cands[0])
            else:
                continue
        jobs.append((stem, fpath))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--skip-done", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(CURVE_DIR, exist_ok=True)
    jobs = collect_jobs(set(args.only) if args.only else None)
    if args.skip_done:
        done = {f[:-len("_sixd2.json")] for f in os.listdir(OUT_DIR)
                if f.endswith("_sixd2.json")}
        jobs = [j for j in jobs if j[0] not in done]
    print(f"[FLOW] SONO v3.5b: {len(jobs)} records, {args.workers} workers, "
          f"beads {N_BEADS}/{N_BEADS_DENSE}/max {N_MAX}, "
          f"max {MAX_ITERS} iter, anchor trefoil RL~{RL_TREFOIL}")
    n_ok = n_fail = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_record, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            line = fu.result()
            print(f"[{i}/{len(jobs)}] {line}", flush=True)
            bad = line.startswith("FAIL") or "CIRCLE_COLLAPSE" in line
            n_fail += bad
            n_ok += not bad
    print(f"[FLOW] DONE: {n_ok} ok, {n_fail} fail/collapse")


if __name__ == "__main__":
    main()
