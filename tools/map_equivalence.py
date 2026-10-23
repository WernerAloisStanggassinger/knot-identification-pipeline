# Translation gate for the English release (23.10.2026).
# Proves that the working tree differs from a git revision that was already proven equivalent to the authors' originals
# (tools/strict_equivalence.py, tools/path_equivalence.py) ONLY by a consistent renaming, and writes that renaming out as
# the translation table. Comments and docstrings are free. Everything else is compared node by node; the differences
# allowed are:
#   - identifiers (names, attributes, parameters, functions, classes, keyword arguments, imported module names), if the
#     renaming is the same everywhere in the repository, one-to-one, does not hit a name that already exists, and only
#     touches attributes / keyword arguments that the repository itself defines (never those of a third-party library);
#   - string constants: a string that appears anywhere in a KEY position (subscript, dict key, comparison, membership,
#     .get/.pop/.setdefault/getenv, add_argument, endswith/startswith, path/file names, os.environ) must be renamed
#     everywhere or nowhere - no old occurrence may survive in any .py/.sh/.json file; a string only ever used as a VALUE
#     (free text written into results or printed) may change freely and is listed in the ledger;
#   - strings that are embedded Python code are compared recursively under the same rules.
# Renamed files are paired by --rename OLD_PATH=NEW_PATH. Exit 0 = green.
# Usage: python3 -I tools/map_equivalence.py REV [--rename old=new ...] [--out tools/translation_table.json]
import ast, builtins, json, os, re, subprocess, sys, textwrap

KEY_CALLS = {'get', 'pop', 'setdefault', 'getenv', 'add_argument', 'endswith', 'startswith', 'open', 'join', 'glob',
             'rglob', 'exists', 'isfile', 'isdir', 'Path', 'with_suffix', 'index', 'count', 'remove', 'getattr',
             'hasattr', 'setattr', 'run', 'Popen', 'check_output', 'check_call', 'call', 'replace', 'split', 'rsplit',
             'find', 'rfind', 'removesuffix', 'removeprefix', 'match', 'search', 'fullmatch', 'compile', 'sub'}
IDENT_FIELDS = {ast.Name: 'id', ast.Attribute: 'attr', ast.arg: 'arg', ast.FunctionDef: 'name',
                ast.AsyncFunctionDef: 'name', ast.ClassDef: 'name', ast.keyword: 'arg', ast.alias: 'name'}

def strip_docstrings(tree):
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.body:
            b = n.body[0]
            if isinstance(b, ast.Expr) and isinstance(b.value, ast.Constant) and isinstance(b.value.value, str):
                b.value.value = '<DOC>'
    return tree

def parse_code(s):
    if '\n' not in s: return None
    for src in (s, textwrap.dedent(s)):
        try: t = ast.parse(src)
        except SyntaxError: continue
        if any(isinstance(n, (ast.FunctionDef, ast.Import, ast.ImportFrom, ast.For, ast.While, ast.If, ast.Assign)) for n in ast.walk(t)):
            return strip_docstrings(t)
    return None

def key_strings(tree, acc):
    """Every string constant that sits in a key position."""
    def add(node):
        for c in ast.walk(node) if isinstance(node, (ast.Tuple, ast.List, ast.Set, ast.BinOp, ast.JoinedStr)) else [node]:
            if isinstance(c, ast.Constant) and isinstance(c.value, str): acc.add(c.value)
    for n in ast.walk(tree):
        if isinstance(n, ast.Subscript): add(n.slice)
        elif isinstance(n, ast.Dict): [add(k) for k in n.keys if k is not None]
        elif isinstance(n, ast.Compare): [add(x) for x in [n.left] + n.comparators]
        elif isinstance(n, ast.Call):
            f = n.func; name = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
            if name in KEY_CALLS: [add(a) for a in n.args]
            for k in n.keywords:
                if k.arg in ('dest', 'default', 'choices', 'prefix', 'suffix'): add(k.value)
        elif isinstance(n, ast.Constant) and isinstance(n.value, str):
            code = parse_code(n.value)
            if code is not None: key_strings(code, acc)
    return acc

def defined_names(tree, acc):
    """Names the repository itself defines (functions, classes, parameters, assigned attributes and names)."""
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)): acc.add(n.name)
        elif isinstance(n, ast.arg): acc.add(n.arg)
        elif isinstance(n, (ast.Assign, ast.AnnAssign, ast.AugAssign, ast.For, ast.comprehension, ast.With, ast.NamedExpr)):
            for t in ast.walk(n):
                if isinstance(t, ast.Name) and isinstance(t.ctx, ast.Store): acc.add(t.id)
                if isinstance(t, ast.Attribute) and isinstance(t.ctx, ast.Store): acc.add(t.attr)
        elif isinstance(n, ast.Constant) and isinstance(n.value, str):
            code = parse_code(n.value)
            if code is not None: defined_names(code, acc)
    return acc

def strip_comments_loose(s):
    """Remove '#' comments from code-like text that may not parse on its own (a piece of a template)."""
    out = []
    for line in s.split('\n'):
        q = None; i = 0; cut = len(line)
        while i < len(line):
            c = line[i]
            if q:
                if c == '\\': i += 2; continue
                if line.startswith(q, i): i += len(q); q = None; continue
            else:
                if line.startswith(('"""', "'''"), i): q = line[i:i + 3]; i += 3; continue
                if c in '"\'': q = c
                elif c == '#': cut = i; break
            i += 1
        out.append(line[:cut].rstrip())
    return '\n'.join(l for l in out if l.strip())

CODE_LINE = re.compile(r'^\s*(def |class |import |from \S+ import |for |while |if |elif |else:|try:|except|finally:|return\b|with |raise |'
                       r'del |assert |pass\b|yield\b|break\b|continue\b|global |nonlocal |lambda\b|[A-Za-z_][\w.]*\(|'
                       r'[A-Za-z_][\w.\[\]\'"]*\s*[-+*/]?=[^=])')

def looks_like_code(s):
    return '\n' in s and any(CODE_LINE.match(l) for l in s.split('\n'))

def template(js):
    return ''.join(v.value if isinstance(v, ast.Constant) else f'__FV{i}__' for i, v in enumerate(js.values))

def argparse_dests(tree, acc):
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'add_argument':
            dest = next((k.value.value for k in n.keywords if k.arg == 'dest' and isinstance(k.value, ast.Constant)), None)
            flags = [a.value for a in n.args if isinstance(a, ast.Constant) and isinstance(a.value, str)]
            if dest is None and flags:
                longf = [f for f in flags if f.startswith('--')]
                dest = (longf[0] if longf else flags[0]).lstrip('-').replace('-', '_')
            if dest: acc.add(dest)
    return acc

def repo_functions(tree, acc):
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            acc[n.name] = {a.arg for a in n.args.args + n.args.kwonlyargs + n.args.posonlyargs}
        elif isinstance(n, ast.ClassDef):
            acc.setdefault(n.name, set())
    return acc

def keyword_callees(tree, acc):
    for n in ast.walk(tree):
        if isinstance(n, ast.Call):
            f = n.func; name = f.id if isinstance(f, ast.Name) else (f.attr if isinstance(f, ast.Attribute) else None)
            for k in n.keywords:
                if k.arg: acc.setdefault(k.arg, set()).add(name)
    return acc

def name_strings(tree, acc):
    """Strings that spell an identifier: globals()/locals()/vars()[...] and getattr/setattr/hasattr/delattr(obj, 'name')."""
    for n in ast.walk(tree):
        if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Call) and isinstance(n.value.func, ast.Name) \
                and n.value.func.id in ('globals', 'locals', 'vars') and isinstance(n.slice, ast.Constant):
            acc.append(n.slice)
        elif isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in ('getattr', 'setattr', 'hasattr', 'delattr') \
                and len(n.args) > 1 and isinstance(n.args[1], ast.Constant):
            acc.append(n.args[1])
    return acc

class Walk:
    def __init__(self, fname): self.f = fname; self.idents = []; self.strings = []; self.errors = []
    def cmp(self, a, b, path='root'):
        if type(a) is not type(b):
            self.errors.append(f'{self.f}: {path}: node {type(a).__name__} vs {type(b).__name__}'); return
        if isinstance(a, list):
            if len(a) != len(b): self.errors.append(f'{self.f}: {path}: list length {len(a)} vs {len(b)}'); return
            for i, (x, y) in enumerate(zip(a, b)): self.cmp(x, y, f'{path}[{i}]')
            return
        if not isinstance(a, ast.AST):
            if a != b: self.errors.append(f'{self.f}: {path}: {a!r} vs {b!r}')
            return
        if isinstance(a, ast.Constant):
            if isinstance(a.value, str) and isinstance(b.value, str):
                if a.value != b.value:
                    ca, cb = parse_code(a.value), parse_code(b.value)
                    if ca is not None and cb is not None: self.cmp(ca, cb, path + '<code>')
                    elif looks_like_code(a.value) or looks_like_code(b.value):
                        if strip_comments_loose(a.value) != strip_comments_loose(b.value):
                            self.errors.append(f'{self.f}: {path}: multi-line string changed beyond comments')
                    else: self.strings.append((a.value, b.value, getattr(a, 'lineno', 0)))
            elif type(a.value) is not type(b.value) or a.value != b.value:
                self.errors.append(f'{self.f}: {path}: constant {a.value!r} vs {b.value!r}')
            return
        if isinstance(a, ast.JoinedStr):
            fa = [v for v in a.values if not isinstance(v, ast.Constant)]; fb = [v for v in b.values if not isinstance(v, ast.Constant)]
            self.cmp(fa, fb, path + '.fields')
            ta, tb = template(a), template(b)
            if ta != tb:
                ca, cb = parse_code(ta), parse_code(tb)
                if ca is not None and cb is not None: self.cmp(ca, cb, path + '<fcode>')
                elif looks_like_code(ta) or looks_like_code(tb):
                    if strip_comments_loose(ta) != strip_comments_loose(tb):
                        self.errors.append(f'{self.f}: {path}: multi-line f-string template changed beyond comments')
                else: self.strings.append((ta, tb, getattr(a, 'lineno', 0)))
            return
        idf = IDENT_FIELDS.get(type(a))
        if isinstance(a, ast.ImportFrom):
            self.idents.append((a.module, b.module, 'module', getattr(a, 'lineno', 0)))
        for field in a._fields:
            x, y = getattr(a, field, None), getattr(b, field, None)
            if field == idf or (isinstance(a, ast.alias) and field == 'asname') or (isinstance(a, ast.ImportFrom) and field == 'module') \
                    or (isinstance(a, (ast.Global, ast.Nonlocal)) and field == 'names'):
                if isinstance(a, ast.ImportFrom) and field == 'module': continue
                xs, ys = (x, y) if isinstance(x, list) else ([x], [y])
                kind = type(a).__name__ + ('.asname' if field == 'asname' else '')
                for u, v in zip(xs, ys):
                    if u is not None or v is not None: self.idents.append((u, v, kind, getattr(a, 'lineno', 0)))
                continue
            self.cmp(x, y, f'{path}.{field}')

def git_show(rev, path):
    r = subprocess.run(['git', 'show', f'{rev}:{path}'], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None

def main():
    rev = sys.argv[1]; renames = {}; out = 'tools/translation_table.json'
    args = sys.argv[2:]
    while args:
        a = args.pop(0)
        if a == '--rename': o, n = args.pop(0).split('=', 1); renames[o] = n
        elif a == '--out': out = args.pop(0)
    old_files = subprocess.run(['git', 'ls-tree', '-r', '--name-only', rev], capture_output=True, text=True).stdout.split()
    py_old = [p for p in old_files if p.endswith('.py') and p != 'tools/map_equivalence.py']
    errors, idents, strings = [], [], []
    old_trees, new_trees = {}, {}
    for p in py_old:
        q = renames.get(p, p)
        if not os.path.exists(q): errors.append(f'{p}: missing in working tree as {q}'); continue
        ta = strip_docstrings(ast.parse(git_show(rev, p))); tb = strip_docstrings(ast.parse(open(q).read()))
        old_trees[p], new_trees[q] = ta, tb
        w = Walk(q); w.cmp(ta, tb); errors += w.errors
        idents += [(u, v, k, q, l) for u, v, k, l in w.idents]; strings += [(u, v, q, l) for u, v, l in w.strings]
    # identifiers: one-to-one, uniform, fresh, repository-defined
    fwd, bwd = {}, {}
    for u, v, k, q, l in idents:
        if u is None or v is None: continue
        if fwd.setdefault(u, v) != v: errors.append(f'{q}:{l}: identifier {u!r} renamed to {v!r} here but to {fwd[u]!r} elsewhere')
        if bwd.setdefault((q, v), u) != u: errors.append(f'{q}:{l}: two identifiers ({bwd[(q, v)]!r}, {u!r}) renamed to the same {v!r}')
    all_old_idents = {u for u, v, k, q, l in idents if u}
    defined = set()
    for t in old_trees.values(): defined_names(t, defined)
    ren = {u: v for u, v in fwd.items() if u != v}
    undefined = []   # decided after the keys are known: an argparse dest or a data key renamed together is the repository's own
    per_file = {}
    for u, v, k, q, l in idents:
        if u: per_file.setdefault(q, set()).add(u)
    for q, names in per_file.items():
        for u in names:
            v = ren.get(u)
            if v and v in names and fwd.get(v, v) == v:
                errors.append(f'{q}: identifier {u!r} -> {v!r}: {v!r} already exists in this file and keeps its name (collision)')
    for u, v in ren.items():
        if hasattr(builtins, v): errors.append(f'identifier {u!r} -> {v!r}: shadows a builtin')
        repo_modules = {os.path.basename(x)[:-3] for x in py_old}
        if u not in defined and u not in repo_modules:
            undefined.append((u, v))
    # imports of repository modules must point at a module file that exists in the working tree
    mods_old = {os.path.basename(x)[:-3] for x in py_old}
    mods_new = {os.path.basename(renames.get(x, x))[:-3] for x in py_old if os.path.exists(renames.get(x, x))}
    for q, t in new_trees.items():
        for n in ast.walk(t):
            names = [a.name.split('.')[0] for a in n.names] if isinstance(n, ast.Import) else \
                    ([n.module.split('.')[0]] if isinstance(n, ast.ImportFrom) and n.module and n.level == 0 else [])
            for m in names:
                if m in mods_old and m not in mods_new:
                    errors.append(f'{q}:{n.lineno}: imports {m!r}, but that module file no longer exists (renamed?)')
    # argparse: an attribute that is an argparse dest may only be renamed together with its flag
    dests_old, dests_new, attrs_new = set(), set(), set()
    for t in old_trees.values(): argparse_dests(t, dests_old)
    for t in new_trees.values():
        argparse_dests(t, dests_new); attrs_new |= {n.attr for n in ast.walk(t) if isinstance(n, ast.Attribute)}
    for u, v in ren.items():
        if u in dests_old and v not in dests_new: errors.append(f'identifier {u!r} -> {v!r}: {u!r} is an argparse dest, but no option now yields {v!r}')
    for d in dests_old - dests_new:
        if d in attrs_new: errors.append(f'argparse dest {d!r} is gone (option renamed) but the attribute {d!r} is still read')
    # keyword arguments: renaming allowed only when every callee is a repository function with that parameter
    funcs, kwc = {}, {}
    for t in old_trees.values(): repo_functions(t, funcs); keyword_callees(t, kwc)
    # strings that spell identifiers must follow the identifier renaming
    for p in py_old:
        q = renames.get(p, p)
        if p not in old_trees: continue
        na, nb = name_strings(old_trees[p], []), name_strings(new_trees[q], [])
        for x, y in zip(na, nb):
            if isinstance(x.value, str) and x.value in ren and y.value != ren[x.value]:
                errors.append(f'{q}:{x.lineno}: string {x.value!r} names the identifier renamed to {ren[x.value]!r}, but reads {y.value!r}')
    # strings
    keys_old = set()
    for t in old_trees.values(): key_strings(t, keys_old)
    sfwd, ledger = {}, {}
    for u, v, q, l in strings:
        if u in keys_old:
            if sfwd.setdefault(u, v) != v: errors.append(f'{q}:{l}: key {u!r} renamed to {v!r} here but {sfwd[u]!r} elsewhere')
        else: ledger.setdefault(u, v)
    for u, v in sfwd.items():
        if v in keys_old and sfwd.get(v, v) == v: errors.append(f'key {u!r} -> {v!r}: {v!r} already used as a key (collision)')
    # (kwarg check needs the key renaming: dict(sorte=...) is a data key and may move only together with the key 'sorte')
    for u, v in ren.items():
        for callee in kwc.get(u, ()):
            if (callee not in funcs or (u not in funcs[callee] and funcs[callee])) and sfwd.get(u) != v:
                errors.append(f'keyword argument {u!r} -> {v!r} passed to {callee!r}, which is not a repository function with '
                              f'that parameter (a data key? then rename the key {u!r} -> {v!r} too)')
    for u, v in undefined:
        if u in dests_old or sfwd.get(u) == v: continue
        errors.append(f'identifier {u!r} -> {v!r}: not defined by the repository (third-party attribute or keyword?)')
    tracked = subprocess.run(['git', 'ls-files'], capture_output=True, text=True).stdout.split() + list(renames.values())
    for f in sorted(set(tracked)):
        if not f.endswith(('.py', '.sh', '.json')) or not os.path.exists(f) or f.startswith('tools/translation_table'): continue
        src = open(f, errors='replace').read()
        if f.endswith('.py'):
            consts = set()
            try:
                for n in ast.walk(ast.parse(src)):
                    if isinstance(n, ast.Constant) and isinstance(n.value, str):
                        consts.add(n.value); c = parse_code(n.value)
                        if c is not None: consts |= {m.value for m in ast.walk(c) if isinstance(m, ast.Constant) and isinstance(m.value, str)}
            except SyntaxError as e: errors.append(f'{f}: syntax error {e}'); continue
            for u in sfwd:
                if sfwd[u] != u and u in consts: errors.append(f'{f}: old key {u!r} still present (renamed to {sfwd[u]!r} elsewhere)')
        else:
            for u in sfwd:
                if sfwd[u] != u and re.search(r'(?<![A-Za-z0-9_])' + re.escape(u) + r'(?![A-Za-z0-9_])', src):
                    errors.append(f'{f}: old key {u!r} still present (renamed to {sfwd[u]!r})')
    table = {'revision': rev, 'files': {o: n for o, n in renames.items()},
             'identifiers': dict(sorted(ren.items())),
             'keys': {u: v for u, v in sorted(sfwd.items()) if u != v},
             'free_text': {u: v for u, v in ledger.items()}}
    if errors:
        for e in errors[:60]: print('RED', e)
        print(f'DIFFERENT: {len(errors)} problems'); sys.exit(1)
    json.dump(table, open(out, 'w'), indent=1, ensure_ascii=False)
    print(f'EQUIVALENT modulo renaming: {len(py_old)} files, {len(table["identifiers"])} identifiers, {len(table["keys"])} keys, '
          f'{len(table["free_text"])} free texts, {len(renames)} file renames -> {out}')

if __name__ == '__main__': main()
