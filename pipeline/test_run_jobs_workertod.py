# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# Counter-case for run_jobs() from eich_vollkorpus_robust.py (Disc #704).
# Four cases. The predecessor (mp.Pool.imap_unordered) fails at B and C:
# it hangs silently. Exactly that must NOT happen any more here.
import os, signal, sys, time
import site_config
sys.path.insert(0, site_config.SIXD_DIR)
from eich_vollkorpus_robust import run_jobs

MARKE = "/tmp/gegenfall_runjobs_marke_"

def brav(job):
    return job[0], 0.1, False, f"CONVERGED {job[0]}"

def stirbt_einmal(job):
    if job[0] == "toedlich":
        m = MARKE + "einmal"
        if not os.path.exists(m):
            open(m, "w").close()
            time.sleep(0.3)
            os.kill(os.getpid(), signal.SIGTERM)   # earlyoom equivalent
            time.sleep(5)
        return job[0], 0.1, False, f"CONVERGED {job[0]} (zweiter Versuch)"
    return brav(job)

def stirbt_immer(job):
    if job[0] == "unrettbar":
        time.sleep(0.2)
        os.kill(os.getpid(), signal.SIGTERM)
        time.sleep(5)
    return brav(job)

def wirft(job):
    if job[0] == "wirft":
        raise ValueError("kaputte Kurve")
    return brav(job)

def jobs(namen):
    return [(n, "/kein/pfad", 0) for n in namen]

if __name__ == "__main__":
    fehler = 0
    def pruefe(name, ist, soll):
        global fehler
        ok = ist == soll
        fehler += (not ok)
        print(f"   {'✅' if ok else '🔴'} {name}: {ist}  (expected {soll})")

    for f in ("einmal",):
        try: os.remove(MARKE + f)
        except OSError: pass

    print("\n=== A  all well-behaved — no restart, no losses ===")
    r = run_jobs(jobs(["a", "b", "c", "d", "e"]), 1, fn=brav)
    pruefe("(ok, fail, aufgegeben, neustarts)", r, (5, 0, 0, 0))

    print("\n=== B  a worker dies ONCE — task must come back ===")
    t0 = time.time()
    r = run_jobs(jobs(["v", "toedlich", "w", "x"]), 1, fn=stirbt_einmal)
    print(f"   duration {time.time()-t0:.1f}s")
    pruefe("(ok, fail, aufgegeben, neustarts)", r, (4, 0, 0, 1))

    print("\n=== C  task kills the worker ALWAYS — must give up, not hang ===")
    t0 = time.time()
    r = run_jobs(jobs(["p", "unrettbar", "q"]), 1, max_versuche=3, fn=stirbt_immer)
    print(f"   duration {time.time()-t0:.1f}s")
    pruefe("ok",           r[0], 2)
    pruefe("aufgegeben",   r[2], 1)
    pruefe("Neustarts>=2", r[3] >= 2, True)

    print("\n=== D  exception IN the worker — count as fail, not as death ===")
    r = run_jobs(jobs(["y", "wirft", "z"]), 1, fn=wirft)
    pruefe("(ok, fail, aufgegeben, neustarts)", r, (2, 1, 0, 0))

    print(f"\n{'✅ ALL FOUR GREEN' if not fehler else f'🔴 {fehler} FAILED'}")
    sys.exit(1 if fehler else 0)
