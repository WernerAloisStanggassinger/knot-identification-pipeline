# Applies a translation map to the release (23.10.2026). Token-based, so layout and comments survive and nothing is
# replaced inside a longer word. The proof that the result is the same program is NOT this script's job - that is
# tools/map_equivalence.py, run afterwards against the last proven revision.
# Map (JSON): {"identifiers": {old: new}, "strings": {old: new}, "keys": {old: new}, "files": {old_path: new_path}}
#   identifiers - Python NAME tokens (also inside f-string fields)
#   strings     - Python string tokens whose value is exactly `old` (one piece of an implicit concatenation counts alone)
#   keys        - like strings, and additionally replaced as whole words in .sh and .json files (data field names,
#                 flags, environment variables, file-name parts that shell scripts also spell out)
#   files       - git mv old -> new
# Usage: python3 -I tools/apply_translation.py MAP.json [--dry-run]
import ast, io, json, os, re, subprocess, sys, tokenize

def literal(tok, new):
    """New source literal in the old token's style: same prefix and quotes when the text needs no escaping, else repr()."""
    m = re.match(r'([rRbBuU]*)(\'\'\'|"""|\'|")', tok)
    prefix, q = m.group(1), m.group(2)
    if '\\' not in new and '\n' not in new and q[0] not in new: return prefix + q + new + q
    return re.sub('[rR]', '', prefix) + repr(new)

def rewrite_py(src, idents, strings):
    lines = src.splitlines(keepends=True); offs = [0]
    for l in lines: offs.append(offs[-1] + len(l))
    edits = []
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        s = None
        if t.type == tokenize.NAME and t.string in idents: s = idents[t.string]
        elif t.type == tokenize.STRING:
            try: v = ast.literal_eval(t.string)
            except Exception: v = None
            if isinstance(v, str) and v in strings: s = literal(t.string, strings[v])
        if s is not None:
            a = offs[t.start[0] - 1] + t.start[1]; b = offs[t.end[0] - 1] + t.end[1]
            assert src[a:b] == t.string, (t, src[a:b])
            edits.append((a, b, s))
    for a, b, s in reversed(edits): src = src[:a] + s + src[b:]
    return src, len(edits)

def main():
    m = json.load(open(sys.argv[1])); dry = '--dry-run' in sys.argv
    idents = m.get('identifiers', {}); keys = m.get('keys', {}); strings = {**m.get('strings', {}), **keys}
    for o, n in m.get('files', {}).items():
        print('git mv', o, n)
        if not dry: subprocess.run(['git', 'mv', o, n], check=True)
    files = subprocess.run(['git', 'ls-files'], capture_output=True, text=True).stdout.split()
    total = 0
    for f in files:
        if f.startswith('tools/') or not os.path.exists(f): continue
        src = open(f).read()
        if f.endswith('.py'):
            new, n = rewrite_py(src, idents, strings)
            ast.parse(new)
        elif f.endswith(('.sh', '.json')):
            new, n = src, 0
            for o, v in sorted(keys.items(), key=lambda x: -len(x[0])):
                new, k = re.subn(r'(?<![A-Za-z0-9_])' + re.escape(o) + r'(?![A-Za-z0-9_])', v, new); n += k
        else: continue
        if n:
            total += n; print(f'{n:5d}  {f}')
            if not dry: open(f, 'w').write(new)
    print('replacements:', total, '(dry run)' if dry else '')

if __name__ == '__main__': main()
