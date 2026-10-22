# Replaces the authors' absolute paths in a released .py file by site_config names (22.10.2026). Token-based, no regex on code:
# - a group of adjacent string literals whose value EQUALS an AUTHORS value -> site_config.<NAME>
# - a group whose value CONTAINS AUTHORS values (embedded code, composed paths) -> placeholders @NAME@ inside the literals and
#   .replace('@NAME@', site_config.NAME) appended after the group (a trailer binds to the whole implicit concatenation)
# Skipped: docstrings, comments, f-strings, bytes. Adds `import site_config` after the first top-level import if anything changed.
# Usage: python3 -I apply_site_config.py FILE.py [--write]     (without --write: list what would change)
# The result must then pass path_equivalence.py against the original.
import ast, io, os, sys, tokenize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from path_equivalence import AUTH

NAMES = sorted(AUTH, key=lambda n: -len(AUTH[n]))   # longest value first: PIPELINE_DIR before FREEZE_DIR before ROOT

def docstring_lines(tree):
    s = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.body:
            b = n.body[0]
            if isinstance(b, ast.Expr) and isinstance(b.value, ast.Constant) and isinstance(b.value.value, str):
                s.update(range(b.lineno, b.end_lineno + 1))
    return s

def groups(tokens):
    """Runs of adjacent plain STRING tokens (only NL/COMMENT between them)."""
    run = []
    for i, t in enumerate(tokens):
        if t.type == tokenize.STRING: run.append(i)
        elif t.type in (tokenize.NL, tokenize.COMMENT) and run: continue
        else:
            if run: yield run
            run = []
    if run: yield run

def name_at(v, k):
    for n in NAMES:
        a = AUTH[n]
        if v.startswith(a, k):
            end = k + len(a)
            if end == len(v) or v[end] in "/'\"\n ": return n
    return None

def main():
    path = sys.argv[1]; write = '--write' in sys.argv
    src = open(path).read(); tree = ast.parse(src); doc = docstring_lines(tree)
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    edits = []   # (start(row,col), end(row,col), new_text)
    for g in groups(toks):
        ts = [toks[i] for i in g]
        if any(t.start[0] in doc for t in ts): continue
        if any(t.string.lstrip('rRuU')[:1] not in ('"', "'") for t in ts): continue   # f-strings, bytes
        val = ''.join(ast.literal_eval(t.string) for t in ts)
        exact = [n for n in NAMES if AUTH[n] == val]
        if exact:
            edits.append((ts[0].start, ts[-1].end, f'site_config.{exact[0]}')); continue
        used = []
        for t in ts:
            v = ast.literal_eval(t.string); new = t.string; k = 0; hit = False
            while k < len(v):
                n = name_at(v, k)
                if n:
                    if AUTH[n] not in new: raise SystemExit(f'{path}:{t.start[0]}: path not literal in token text (escapes?)')
                    new = new.replace(AUTH[n], f'@{n}@', 1); used.append(n); k += len(AUTH[n]); hit = True
                else: k += 1
            if hit: edits.append((t.start, t.end, new))
        if used:
            trailer = ''.join(f".replace('@{n}@', site_config.{n})" for n in dict.fromkeys(used))
            edits.append((ts[-1].end, ts[-1].end, trailer))
    if not edits: print('no paths', path); return
    lines = src.splitlines(keepends=True); offs = [0]
    for l in lines: offs.append(offs[-1] + len(l))
    pos = lambda rc: offs[rc[0] - 1] + rc[1]
    out = src
    for a, b, new in sorted(edits, key=lambda e: (pos(e[0]), pos(e[1])), reverse=True):
        out = out[:pos(a)] + new + out[pos(b):]
    imports = [s for s in tree.body if isinstance(s, (ast.Import, ast.ImportFrom)) and not (isinstance(s, ast.ImportFrom) and s.module == '__future__')]
    if 'import site_config' not in src:
        ol = out.splitlines(keepends=True)
        at = imports[0].end_lineno if imports else 0
        ol.insert(at, 'import site_config\n'); out = ''.join(ol)
    for a, b, new in sorted(edits, key=lambda e: pos(e[0])):
        print(f'{path}:{a[0]}: {new[:110]}')
    if write: open(path, 'w').write(out); print('WRITTEN', path)

if __name__ == '__main__': main()
