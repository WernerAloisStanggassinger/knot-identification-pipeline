# Path gate for the release (22.10.2026): proves that a released file, with site_config in the AUTHORS profile, is the same program
# as the original pipeline file — up to the translation allowed by strict_equivalence (comments, docstrings, message texts).
# How: in BOTH files, every "path expression" is folded to its value — string constants, site_config.<NAME> (-> AUTHORS[NAME]),
# module-level names assigned exactly once to a foldable expression, Path(x) / PurePath(x), x / y on paths, os.path.join(...),
# str + str, and s.replace(a, b) on foldable strings. Then `import site_config` statements are removed from the released file
# and both trees are compared with strict_equivalence's normaliser (embedded code strings recursively).
# Usage: python3 -I path_equivalence.py ORIGINAL.py RELEASED.py [--drop-seal]   exit 0 = equivalent
import ast, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from strict_equivalence import Norm

def authors_table():
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [os.path.join(here, '..', 'pipeline', 'site_config.py'),                                   # inside the repository
             os.path.join(here, '..', 'knot-identification-pipeline', 'pipeline', 'site_config.py')]   # authors' staging layout
    src = open(next(c for c in cands if os.path.exists(c))).read()
    t = ast.parse(src); env = {}
    for st in t.body:  # evaluate only the literal part that defines _K, _CONDA and AUTHORS
        if isinstance(st, ast.Assign) and isinstance(st.targets[0], ast.Name) and st.targets[0].id in ('_K', '_CONDA', 'AUTHORS'):
            env[st.targets[0].id] = eval(compile(ast.Expression(st.value), 'site_config', 'eval'), {}, dict(env))
    return env['AUTHORS']

AUTH = authors_table()

class P(str):
    """A folded pathlib.Path value (kept distinct from a plain str)."""

def assigned_once(tree):
    count = {}
    for n in ast.walk(tree):
        tg = []
        if isinstance(n, ast.Assign): tg = n.targets
        elif isinstance(n, (ast.AugAssign, ast.AnnAssign)): tg = [n.target]
        elif isinstance(n, (ast.For, ast.comprehension)): tg = [n.target]
        elif isinstance(n, (ast.FunctionDef, ast.ClassDef)): count[n.name] = count.get(n.name, 0) + 2
        for t in tg:
            for x in ast.walk(t):
                if isinstance(x, ast.Name): count[x.id] = count.get(x.id, 0) + 1
    return {k for k, v in count.items() if v == 1}

def fold(node, env):
    """Value of a foldable path expression, or None."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str): return node.value
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == 'site_config':
        return AUTH[node.attr]
    if isinstance(node, ast.Name) and node.id in env: return env[node.id]
    if isinstance(node, ast.Call):
        f = node.func; name = ast.unparse(f)
        args = [fold(a, env) for a in node.args]
        if node.keywords or any(a is None for a in args): return None
        if name in ('Path', 'PurePath', 'pathlib.Path') and len(args) == 1: return P(args[0])
        if name == 'os.path.join' and args: return os.path.join(*args)
        if name == 'str' and len(args) == 1: return str(args[0])
        if isinstance(f, ast.Attribute) and f.attr == 'replace' and len(args) == 2:
            base = fold(f.value, env)
            if base is not None and not isinstance(base, P): return base.replace(args[0], args[1])
        return None
    if isinstance(node, ast.BinOp):
        a, b = fold(node.left, env), fold(node.right, env)
        if a is None or b is None: return None
        if isinstance(node.op, ast.Div) and isinstance(a, P): return P(os.path.join(a, b))
        if isinstance(node.op, ast.Add) and not isinstance(a, P) and not isinstance(b, P): return a + b
    return None

class Folder(ast.NodeTransformer):
    def __init__(self, env): self.env = env
    def generic_visit(self, node):
        if isinstance(node, ast.expr) and not isinstance(node, ast.Constant):
            v = fold(node, self.env)
            if v is not None:
                tag = 'PATH:' if isinstance(v, P) else ''
                return ast.copy_location(ast.Constant(tag + v if tag else v), node)
        return super().generic_visit(node)
    def visit_Assign(self, node):  # never fold assignment targets
        node.value = self.visit(node.value); return node

def fold_tree(tree):
    once = assigned_once(tree); env = {}
    for st in tree.body:  # module-level names assigned once to a foldable expression
        if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name) and st.targets[0].id in once:
            v = fold(st.value, env)
            if v is not None: env[st.targets[0].id] = v
    return Folder(env).visit(tree)

def drop_site_config_import(tree):
    def is_imp(s): return isinstance(s, ast.Import) and [a.name for a in s.names] == ['site_config'] and s.names[0].asname is None
    for n in ast.walk(tree):
        for field in ('body', 'orelse', 'finalbody'):
            b = getattr(n, field, None)
            if isinstance(b, list): setattr(n, field, [s for s in b if not is_imp(s)] or b)
    return tree

def main():
    a, b = sys.argv[1], sys.argv[2]; drop = '--drop-seal' in sys.argv
    ta = fold_tree(ast.parse(open(a).read()))
    tb = fold_tree(drop_site_config_import(ast.parse(open(b).read())))
    A = Norm(drop).run(ta); B = Norm(False).run(tb)
    if A == B: print('EQUIVALENT', b); sys.exit(0)
    import difflib
    for l in list(difflib.unified_diff(A.split(', '), B.split(', '), lineterm='', n=2))[:40]: print(l)
    print('DIFFERENT', b); sys.exit(1)

if __name__ == '__main__': main()
