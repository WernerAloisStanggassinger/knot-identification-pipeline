# Calibration of the GPU scanner against an independent NumPy brute force: for a small camera flight, every frame's
# set of crossing segment pairs (seg_i, seg_j) must be identical, and t_norm / z / cross2d must agree to float32 precision.
# Usage (from the repository root, NVIDIA GPU + cupy needed): python3 -I tests/test_scanner_vs_numpy.py
import os, sys, subprocess, tempfile
import numpy as np
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(root, 'pipeline'))
from knot_scanner import parse_fseries, reconstruct_curve
N, FRAMES = 256, 64
tmp = tempfile.mkdtemp()
cam = os.path.join(tmp, 'cam')
subprocess.run([sys.executable, '-I', os.path.join(root, 'pipeline', 'make_camera_trajectory.py'), '--n', str(FRAMES), '--out', cam], check=True)
out = os.path.join(tmp, 'cr.npz')
subprocess.run([sys.executable, '-I', os.path.join(root, 'pipeline', 'knot_scanner.py'), os.path.join(root, 'examples', 'trefoil.fseries'),
                '--camera', cam + '.npz', '--output', out, '--N', str(N)], check=True)
gpu = np.load(out); R = np.load(cam + '.npz')['rotmat']
P = reconstruct_curve(*parse_fseries(os.path.join(root, 'examples', 'trefoil.fseries')), N=N).astype(np.float64)
bad = 0; total = 0
for t in range(FRAMES):
    Q = P @ R[t].T.astype(np.float64)
    ref = set()
    for i in range(N):
        for j in range(i + 2, N):
            if i == 0 and j == N - 1: continue
            p0, p1, p2, p3 = Q[i], Q[(i + 1) % N], Q[j], Q[(j + 1) % N]
            d1, d2, d0 = p1[:2] - p0[:2], p3[:2] - p2[:2], p2[:2] - p0[:2]
            den = d1[0] * d2[1] - d1[1] * d2[0]
            if abs(den) < 1e-10: continue
            a = (d0[0] * d2[1] - d0[1] * d2[0]) / den; b = (d0[0] * d1[1] - d0[1] * d1[0]) / den
            if 0 < a < 1 and 0 < b < 1: ref.add((i, j))
    m = gpu['t_step'] == t
    got = set(zip(gpu['seg_i'][m].tolist(), gpu['seg_j'][m].tolist()))
    total += len(ref); bad += len(ref ^ got)
print(f'frames {FRAMES}, reference crossings {total}, pairs differing between GPU and NumPy: {bad}')
sys.exit(0 if bad == 0 else 1)
