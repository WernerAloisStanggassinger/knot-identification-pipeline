# Strict equivalence gate for the English release of the pipeline scripts (21.10.2026).
# A translated file may differ from its original ONLY in:
#   - comments (invisible to the AST),
#   - docstrings (first statement of module/class/function if it is a string),
#   - the literal text of "message" strings: arguments of print(), log*/warn*/error*/info*/debug*(), argparse help=/description=/epilog=,
#     and of exception constructors in `raise` statements, including the literal parts of f-strings there.
# Every OTHER string constant (dict keys, file names, column names, JSON fields, SQL, regexes, embedded code ...) must be
# byte-identical. A string that parses as Python and contains a def/import/for/while/if is treated as EMBEDDED CODE and
# compared recursively with the same rules. Exit 0 = equivalent.
# Usage: python3 -I strict_equivalence.py ORIGINAL.py TRANSLATED.py [--drop-seal]
import ast, sys

MSG_FUNCS = {'print', 'exit'}   # sys.exit("text") prints the text and exits with status 1; integer exit codes stay compared
MSG_PREFIXES = ('log', 'warn', 'error', 'info', 'debug', 'err', 'fail', 'die')
MSG_KW = {'help', 'description', 'epilog'}

def is_code_string(s):
    if '\n' not in s: return False
    try: t = ast.parse(s)
    except SyntaxError:
        try: t = ast.parse(__import__('textwrap').dedent(s))
        except SyntaxError: return False
    return any(isinstance(n, (ast.FunctionDef, ast.Import, ast.ImportFrom, ast.For, ast.While, ast.If, ast.Assign)) for n in ast.walk(t))

def norm_code_string(s, drop_seal=False):
    import textwrap
    try: t = ast.parse(s)
    except SyntaxError: t = ast.parse(textwrap.dedent(s))
    return 'CODE[' + Norm(drop_seal).run(t) + ']'

class Norm(ast.NodeTransformer):
    def __init__(self, drop_seal=False): self.drop_seal = drop_seal; self.msg_depth = 0
    def run(self, tree):
        if self.drop_seal:
            def seal(st):
                u = ast.unparse(st)
                return 'seal_all' in u or 'pipeline_seal' in u or ('sys.path.insert' in u and 'freeze_manifest' in u)
            tree.body = [s for s in tree.body if not seal(s)]
        self.strip_docstrings(tree)
        return ast.dump(self.visit(tree))
    def strip_docstrings(self, tree):
        for n in ast.walk(tree):
            if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.body:
                b = n.body[0]
                if isinstance(b, ast.Expr) and isinstance(b.value, ast.Constant) and isinstance(b.value.value, str):
                    b.value.value = '<DOC>'
    def is_msg_call(self, n):
        f = n.func
        name = f.id if isinstance(f, ast.Name) else (f.attr if isinstance(f, ast.Attribute) else '')
        return name in MSG_FUNCS or name.lower().startswith(MSG_PREFIXES)
    def visit_Call(self, n):
        msg = self.is_msg_call(n)
        n.func = self.visit(n.func)
        if msg: self.msg_depth += 1
        n.args = [self.visit(a) for a in n.args]
        if msg: self.msg_depth -= 1
        new_kw = []
        for k in n.keywords:
            if k.arg in MSG_KW: self.msg_depth += 1; k.value = self.visit(k.value); self.msg_depth -= 1
            else: k.value = self.visit(k.value)
            new_kw.append(k)
        n.keywords = new_kw
        return n
    def visit_Raise(self, n):
        self.msg_depth += 1; n = self.generic_visit(n); self.msg_depth -= 1; return n
    def visit_JoinedStr(self, n):
        if self.msg_depth:
            return ast.copy_location(ast.JoinedStr(values=[self.visit(v) for v in n.values if isinstance(v, ast.FormattedValue)]), n)
        return self.generic_visit(n)
    def visit_Constant(self, n):
        if isinstance(n.value, str):
            if self.msg_depth: return ast.copy_location(ast.Constant('<MSG>'), n)
            if is_code_string(n.value): return ast.copy_location(ast.Constant(norm_code_string(n.value)), n)
        return n

def main():
    a, b = sys.argv[1], sys.argv[2]; drop = '--drop-seal' in sys.argv
    A = Norm(drop).run(ast.parse(open(a).read())); B = Norm(False).run(ast.parse(open(b).read()))
    if A == B: print('EQUIVALENT', b); sys.exit(0)
    import difflib
    for l in list(difflib.unified_diff(A.split(', '), B.split(', '), lineterm='', n=2))[:40]: print(l)
    print('DIFFERENT', b); sys.exit(1)
if __name__ == '__main__': main()
