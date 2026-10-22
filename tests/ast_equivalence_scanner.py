# Proves that the published scanner files differ from the pipeline originals ONLY in string constants
# (docstrings, messages, help texts) and comments, plus the removed blind-test seal (3 statements).
# Usage: python3 -I tests/ast_equivalence.py ORIGINAL_DIR PUBLISHED_DIR
import ast, sys, os
class Blank(ast.NodeTransformer):
    def visit_Constant(self, n):
        return ast.copy_location(ast.Constant('' if isinstance(n.value, str) else n.value), n)
    def visit_JoinedStr(self, n):                     # f-strings: keep the formatted expressions, drop the literal text
        vals = [self.visit(v) for v in n.values if isinstance(v, ast.FormattedValue)]
        return ast.copy_location(ast.JoinedStr(values=vals), n)
def norm(src, drop_seal):
    t = ast.parse(src)
    if drop_seal:
        def is_seal(s):
            d = ast.dump(s)
            return 'seal_all' in d or 'pipeline_seal' in d or ("sys.path" in ast.unparse(s) and 'insert' in ast.unparse(s))
        t.body = [s for s in t.body if not is_seal(s)]
    return ast.dump(Blank().visit(t))
ok = True
for f in ('knot_scanner.py', 'make_camera_trajectory.py'):
    a = norm(open(os.path.join(sys.argv[1], f)).read(), True)
    b = norm(open(os.path.join(sys.argv[2], f)).read(), False)
    print(f, 'code-identical (strings blanked, seal removed):', a == b); ok &= a == b
# The CUDA kernel is a string constant, so the check above does not see it: compare it byte for byte.
ns = {}
for d, tag in ((sys.argv[1], 'orig'), (sys.argv[2], 'pub')):
    src = open(os.path.join(d, 'knot_scanner.py')).read()
    t = ast.parse(src)
    ns[tag] = [n.value.value for n in t.body if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', '') == 'CROSSING_KERNEL'][0]
print('CUDA kernel byte-identical:', ns['orig'] == ns['pub']); ok &= ns['orig'] == ns['pub']
sys.exit(0 if ok else 1)
