# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""collide.py — which volume-equal knot pairs does θ separate?

Input:
  values  one or more files <prefix>_values.tsv.bz2 written by research/theta_census_gpu/theta_gpu_census.py --values
          (tab-separated columns: name, diagram crossings, θ_raw(K) mod p, θ_raw(mirror K) mod p); several files are given
          as one comma-separated argument;
  TOL     volume tolerance (we use 1e-9);
  volumes the tables of research/volume_census/volume_census.py (name, volume, flag, attempts); only rows with flag `pos`
          (a solution with all tetrahedra positively oriented) and a volume are used.
Mutation preserves the hyperbolic volume (Ruberman 1987), so every mutant pair is a volume-equal pair; not every volume-equal
pair is a mutant pair. Groups: knots sorted by volume, neighbours closer than TOL join one group (chained); every pair inside a
group is counted. θ_raw = θ / (Δ(T1)Δ(T2)Δ(T1T2)) is a knot invariant. θ "separates" A, B when the values prove both A ≠ B and
A ≠ mirror B: (θ_raw(A) ≢ θ_raw(B) or θ_raw(mA) ≢ θ_raw(mB)) and (θ_raw(A) ≢ θ_raw(mB) or θ_raw(mA) ≢ θ_raw(B)) mod p
(mirror-blind, as census names ignore mirror images). A difference mod p proves that θ_raw differs (reduction mod p), given a
correct kernel and diagram conversion; for pairs with equal Alexander polynomial (e.g. mutants) it proves that θ differs.
Agreement is only "equal at the point modulo p".
Output: <out>_stats.json (totals, per pair of crossing numbers) and <out>_not_separated.jsonl (one line per pair θ does not
separate, with theta_zero_both).
Counter-case built in: if 15nh_048139 and 15nh_048140 (a pair with θ = 0 and equal volume at 15 crossings) end up in one
volume group, they must be reported as a non-separated pair (assert); otherwise the counter-case is recorded as not applicable.
theta_zero_both looks at θ_raw(K) of both knots (under the sign rule θ(mK) = −θ(K) the mirror is then 0 too).
Separation is between ORIENTED knots: each census diagram fixes one orientation; A ≠ −B would need θ(−K) = θ(K), which
Bar-Natan and van der Veen state as a conjecture.
Usage: python -I collide.py values1.tsv.bz2[,values2.tsv.bz2,...] TOL out_prefix vol_11n.tsv vol_12n.tsv ...
       python -I collide.py --selftest
"""
import bz2, json, sys, collections, os, tempfile, math

def run(vals_paths, tol, pre, vol_paths):
    # pass 1: names with θ (a set of all names); pass 2: volumes and groups (all volumes in memory); pass 3: θ values of
    # group members only — this saves the θ strings, not the names and volumes (tens of millions of knots still need GBs)
    have = set()
    for vp in vals_paths:
        for k, line in enumerate(bz2.open(vp, 'rt'), 1):
            line = line.rstrip('\r\n')
            if line.count('\t') != 3: raise ValueError(f'{vp}:{k}: malformed values row: {line[:80]!r}')
            have.add(line.split('\t', 1)[0])
    vol = []; seen = {}               # a name in two volume tables is used once (and must have the same volume there)
    no_theta = 0
    for p in vol_paths:
        for line in open(p):
            if line.startswith('#') or not line.strip(): continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) < 3: raise ValueError(f'{p}: malformed volume row: {line[:80]!r}')
            if f[0] not in have and f[1] not in ('None', '') and f[2] == 'pos': no_theta += 1
            if f[0] in have and f[1] not in ('None', '') and f[2] == 'pos':
                x = float(f[1])
                if not math.isfinite(x): raise ValueError(f'{p}: non-finite volume for {f[0]}: {f[1]!r}')
                if f[0] in seen:
                    if seen[f[0]] != x: raise ValueError(f'{f[0]}: two different volumes in the volume tables ({seen[f[0]]!r}, {x!r})')
                    continue
                vol.append((x, f[0])); seen[f[0]] = x
    n_have = len(have); del have, seen
    vol.sort()
    groups, cur = [], [vol[0]] if vol else []
    for v in vol[1:]:
        if v[0] - cur[-1][0] < tol: cur.append(v)
        else:
            if len(cur) > 1: groups.append(cur)
            cur = [v]
    if len(cur) > 1: groups.append(cur)
    members = {n for g in groups for _, n in g}
    th = {}
    for vp in vals_paths:
        for line in bz2.open(vp, 'rt'):
            name, n, a, b = line.rstrip('\r\n').split('\t')
            if name in members:
                if name in th and th[name] != (a, b): raise ValueError(f'{name}: different θ values in two values files')
                th[name] = (a, b)
    cc = ('15nh_048139', '15nh_048140')
    cc_grouped = any(cc[0] in names and cc[1] in names for names in ({n for _, n in g} for g in groups))
    del members

    def sep(a, b):
        ta, tma = th[a]; tb, tmb = th[b]
        if 'None' in (ta, tb, tma, tmb): return None
        return (ta != tb or tma != tmb) and (ta != tmb or tma != tb)      # proves A ≠ B and A ≠ mirror B

    stats = collections.Counter(); blind = []; by_c = collections.defaultdict(collections.Counter)
    def crossings(name):                       # census names such as 15nh_048139, 15ah_…; anything else: '?'
        c = name.split('n')[0].split('a')[0]
        return c if c.isdigit() else '?'
    ckey = lambda a, b: '-'.join(sorted((crossings(a), crossings(b)), key=lambda c: (not c.isdigit(), int(c) if c.isdigit() else 0)))
    # key = the unordered pair of crossing numbers, smaller first
    for g in groups:
        names = [n for _, n in g]
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                r = sep(a, b)
                key = ('separated' if r else 'not_separated') if r is not None else 'none'
                stats[key] += 1; by_c[ckey(a, b)][key] += 1
                if r is False:
                    blind.append(dict(a=a, b=b, va=g[i][0], vb=g[j][0], theta_zero_both=th[a][0] == '0' and th[b][0] == '0'))
    stats['knots_with_volume'] = len(vol); stats['groups'] = len(groups); stats['knots_with_theta'] = n_have
    stats['volume_rows_without_theta'] = no_theta
    if cc_grouped:
        if 'None' in th[cc[0]] + th[cc[1]]:
            stats['countercase_15nh_048139_048140_not_applicable'] = 'a None value'
        else:
            hit = any({x['a'], x['b']} == set(cc) for x in blind)
            stats['countercase_15nh_048139_048140_not_separated'] = int(hit)
            if not hit:
                raise RuntimeError('counter-case failed: the θ = 0 pair 15nh_048139 / 15nh_048140 is not reported as non-separated')
    else:
        stats['countercase_15nh_048139_048140_not_applicable'] = 'not in one volume group'
    json.dump(dict(stats=stats, by_crossings={k: dict(v) for k, v in sorted(by_c.items())}, tol=tol), open(pre + '_stats.json', 'w'), indent=1)
    with open(pre + '_not_separated.jsonl', 'w') as f:
        for x in blind: f.write(json.dumps(x) + '\n')
    return dict(stats), blind

def selftest():
    d = tempfile.mkdtemp(prefix='collide_selftest_')
    vals = os.path.join(d, 'v_values.tsv.bz2'); vol = os.path.join(d, 'vol.tsv')
    with bz2.open(vals, 'wt') as f:      # name, crossings, θ(K), θ(mK)
        for row in [('11nh_1', 11, 5, 7), ('11nh_2', 11, 5, 7),          # equal θ          -> not separated
                    ('12nh_1', 12, 3, 9), ('12nh_2', 12, 9, 3),          # mirror-equal θ   -> not separated
                    ('13nh_1', 13, 1, 2), ('13nh_2', 13, 4, 8),          # different θ      -> separated
                    ('14nh_1', 14, 0, 0), ('14nh_2', 14, 0, 0),          # θ = 0 both       -> not separated, zero flag
                    ('15nh_1', 15, 6, 6), ('15nh_2', 15, 6, 6),          # volumes 2e-9 apart -> not one group at TOL 1e-9
                    ('16nh_1', 16, 'None', 3), ('16nh_2', 16, 3, 4),     # None -> counted as 'none'
                    ('17nh_1', 17, 5, 1), ('18nh_1', 18, 5, 2),          # θ(A) = θ(B) but θ(mA) ≠ θ(mB): separated, key 17-18
                    ('19nh_1', 19, 7, 3), ('19nh_2', 19, 7, 3), ('19nh_3', 19, 7, 3),   # 3-knot chain: steps 6e-10, ends 1.2e-9 apart
                    ('9nh_1', 9, 0, 5), ('10nh_1', 10, 5, 0),            # mirror-equal, θ = 0 on one side only; key 9-10 (numeric order)
                    ('11nh_9', 11, 2, 2), ('11nh_8', 11, 2, 2),          # 11nh_8 has flag 'neg' -> not used
                    ('15nh_048139', 15, 0, 0), ('15nh_048140', 15, 0, 0),    # the built-in counter-case, applicable
                    ('20nh_1', 20, 0, 5), ('20nh_2', 20, 0, 5),          # θ = 0 both, mirrors 5: not separated, flagged
                    ('21nh_1', 21, 1, 2), ('21nh_2', 21, 3, 1)]:         # θ(A) = θ(mB) but θ(mA) ≠ θ(B): separated
            f.write('\t'.join(map(str, row)) + '\n')
    with open(vol, 'w') as f:
        f.write('# synthetic\n')
        for n, v, flag in [('11nh_1', 2.0, 'pos'), ('11nh_2', 2.0 + 5e-10, 'pos'), ('12nh_1', 3.0, 'pos'), ('12nh_2', 3.0, 'pos'),
                           ('13nh_1', 4.0, 'pos'), ('13nh_2', 4.0, 'pos'), ('14nh_1', 5.0, 'pos'), ('14nh_2', 5.0, 'pos'),
                           ('15nh_1', 6.0, 'pos'), ('15nh_2', 6.0 + 2e-9, 'pos'), ('16nh_1', 7.0, 'pos'), ('16nh_2', 7.0, 'pos'),
                           ('18nh_1', 8.0, 'pos'), ('17nh_1', 8.0 + 5e-10, 'pos'),    # 18-crossing knot sorts first
                           ('19nh_1', 9.0, 'pos'), ('19nh_2', 9.0 + 6e-10, 'pos'), ('19nh_3', 9.0 + 1.2e-9, 'pos'),
                           ('10nh_1', 10.0, 'pos'), ('9nh_1', 10.0 + 1e-10, 'pos'), ('11nh_9', 11.0, 'pos'), ('11nh_8', 11.0, 'neg'),
                           ('15nh_048139', 12.0, 'pos'), ('15nh_048140', 12.0 + 1e-11, 'pos'),
                           ('20nh_1', 13.0, 'pos'), ('20nh_2', 13.0, 'pos'), ('21nh_1', 14.0, 'pos'), ('21nh_2', 14.0, 'pos'),
                           ('99nh_1', 15.0, 'pos')]:                              # a volume row without θ
            f.write(f'{n}\t{v!r}\t{flag}\t1\n')
    vals2 = os.path.join(d, 'v2_values.tsv.bz2'); vol2 = os.path.join(d, 'vol2.tsv')
    with bz2.open(vals2, 'wt') as f: f.write('11nh_1\t11\t5\t7\n')          # same name, same values: accepted
    with open(vol2, 'w') as f: f.write('11nh_1\t2.0\tpos\t1\n')                 # same name, same volume: used once
    stats, blind = run([vals, vals2], 1e-9, os.path.join(d, 'out'), [vol, vol2])
    pairs = {frozenset((x['a'], x['b'])): x for x in blind}
    st = json.load(open(os.path.join(d, 'out_stats.json')))
    def raises(fn):
        try: fn(); return False
        except (ValueError, RuntimeError): return True
    vbad = os.path.join(d, 'vbad_values.tsv.bz2')
    with bz2.open(vbad, 'wt') as f: f.write('11nh_1\t11\t6\t7\n')            # same name, DIFFERENT values
    conflict_raises = raises(lambda: run([vals, vbad], 1e-9, os.path.join(d, 'o2'), [vol]))
    vcc = os.path.join(d, 'vcc_values.tsv.bz2'); volcc = os.path.join(d, 'volcc.tsv')
    with bz2.open(vcc, 'wt') as f: f.write('15nh_048139\t15\t1\t2\n15nh_048140\t15\t3\t4\n')   # counter-case pair separated
    with open(volcc, 'w') as f: f.write('15nh_048139\t1.0\tpos\t1\n15nh_048140\t1.0\tpos\t1\n')
    cc_raises = raises(lambda: run([vcc], 1e-9, os.path.join(d, 'o3'), [volcc]))
    vt = os.path.join(d, 'vt_values.tsv.bz2'); volt = os.path.join(d, 'volt.tsv')
    with bz2.open(vt, 'wt') as f: f.write('11nh_1\t11\t1\t1\n11nh_2\t11\t1\t1\n')
    with open(volt, 'w') as f: f.write(f'11nh_1\t{0.5!r}\tpos\t1\n11nh_2\t{0.5 + 2 ** -30!r}\tpos\t1\n')
    strict_tol = run([vt], 2 ** -30, os.path.join(d, 'o4'), [volt])[0]['groups'] == 0     # gap == TOL: not grouped
    import shutil; shutil.rmtree(d, ignore_errors=True)
    checks = {
        'groups == 11': stats['groups'] == 11,
        'separated == 3 (13nh, 17nh/18nh, 21nh)': stats.get('separated') == 3,
        'None pair not among the non-separated': frozenset(('16nh_1', '16nh_2')) not in pairs,
        'θ = 0 flag reads θ_raw(K), not the mirror': pairs.get(frozenset(('20nh_1', '20nh_2')), {}).get('theta_zero_both') is True,
        'duplicate name used once': stats['knots_with_volume'] == 26,
        'volume row without θ counted': stats.get('volume_rows_without_theta') == 1,
        'knots with θ counted (incl. 11nh_8, θ but no pos volume)': stats['knots_with_theta'] == 27,
        'chain: 19nh_1 and 19nh_3 (1.2e-9 apart) in one group': frozenset(('19nh_1', '19nh_3')) in pairs,
        'cross-c key numeric order 9-10': '9-10' in st['by_crossings'] and '10-9' not in st['by_crossings'],
        'one-sided zero not flagged': pairs.get(frozenset(('9nh_1', '10nh_1')), {}).get('theta_zero_both') is False,
        'non-pos row ignored (11nh_8)': not any('11nh_9' in k for k in pairs),
        'counter-case applicable and passed': stats.get('countercase_15nh_048139_048140_not_separated') == 1,
        'conflicting θ values in two files raise': conflict_raises,
        'a failing counter-case raises': cc_raises,
        'a gap equal to TOL does not group (strict <)': strict_tol,
        'none == 1 (16nh)': stats.get('none') == 1,
        'cross-c key smaller crossing number first (17-18)': '17-18' in st['by_crossings'] and '18-17' not in st['by_crossings'],

        'not_separated == 9': stats.get('not_separated') == 9,
        'equal θ not separated': frozenset(('11nh_1', '11nh_2')) in pairs,
        'mirror-equal θ not separated': frozenset(('12nh_1', '12nh_2')) in pairs,
        'θ = 0 flagged': pairs.get(frozenset(('14nh_1', '14nh_2')), {}).get('theta_zero_both') is True,
        '2e-9 apart not grouped': frozenset(('15nh_1', '15nh_2')) not in pairs,
    }
    bad = [k for k, ok in checks.items() if not ok]
    print('selftest', 'OK' if not bad else f'FAILED {bad}', f'({len(checks)} checks)')
    return not bad

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        sys.exit(0 if selftest() else 1)
    st, _ = run(sys.argv[1].split(','), float(sys.argv[2]), sys.argv[3], sys.argv[4:])
    print(json.dumps(st), flush=True)
