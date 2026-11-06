#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Recursion-free evaluator for the huge theta_full expressions.

Reason: `eval()` compiles `a + b + c + …` into a LEFT-DEEP tree; with 20,000 summands
the CPython *compiler* dies on its C stack. `sys.setrecursionlimit` does NOT help there
(that is the interpreter stack, not the compiler's) — exactly why the fix of
26.07. did not close the gap, but only shifted it.

Here: tokenizer + shunting-yard, iterative, exact in Fraction. No recursion, no compile().
"""
import re
from fractions import Fraction as F

_TOK = re.compile(r"\s*(\d+|[A-Za-z]\w*|\*\*|\^|[-+*/()])")


def tokens(s):
    i, n = 0, len(s)
    while i < n:
        m = _TOK.match(s, i)
        if not m:
            if s[i].isspace():
                i += 1
                continue
            raise ValueError("unknown character %r at %d" % (s[i], i))
        i = m.end()
        yield m.group(1)


_PREC = {"+": 1, "-": 1, "*": 2, "/": 2, "u-": 3, "^": 4}
_RECHTS = {"^", "u-"}


def evaluiere(ausdruck, umgebung):
    """Exact value of the expression. umgebung: {'t1': Fraction, 't2': Fraction}."""
    werte, ops = [], []
    potenz = {}   # (name, exponent) -> value, saves a lot with 20k monomials

    def anwenden():
        op = ops.pop()
        if op == "u-":
            werte.append(-werte.pop())
            return
        b = werte.pop()
        a = werte.pop()
        if op == "+":
            werte.append(a + b)
        elif op == "-":
            werte.append(a - b)
        elif op == "*":
            werte.append(a * b)
        elif op == "/":
            werte.append(F(a) / b)
        elif op == "^":
            werte.append(a ** int(b))

    vorher_wert = False   # was there a value directly before? (distinguishes binary from unary -)
    for t in tokens(ausdruck):
        if t.isdigit():
            werte.append(F(int(t)))
            vorher_wert = True
        elif t[0].isalpha():
            if t not in umgebung:
                raise ValueError("unknown name %r (known: %s)"
                                 % (t, ", ".join(sorted(umgebung))))
            werte.append(umgebung[t])
            vorher_wert = True
        elif t == "(":
            ops.append(t)
            vorher_wert = False
        elif t == ")":
            while ops and ops[-1] != "(":
                anwenden()
            if not ops:
                raise ValueError("bracket imbalance")
            ops.pop()
            vorher_wert = True
        else:
            if t == "**":
                t = "^"
            if t == "-" and not vorher_wert:
                t = "u-"
            while (ops and ops[-1] != "(" and
                   (_PREC[ops[-1]] > _PREC[t] or
                    (_PREC[ops[-1]] == _PREC[t] and t not in _RECHTS))):
                anwenden()
            ops.append(t)
            vorher_wert = False
    while ops:
        if ops[-1] == "(":
            raise ValueError("bracket imbalance")
        anwenden()
    if len(werte) != 1:
        raise ValueError("expression not resolved (%d values left)" % len(werte))
    return werte[0]


def _selbsttest():
    """Two duties: (1) on small expressions BIT-IDENTICAL with `eval`,
    (2) on an expression on which `eval` DIES, deliver a result."""
    import sys
    umg = {"t1": F(22, 7), "t2": F(21, 13), "x": F(3)}
    faelle = ["2+3*4", "(2+3)*4", "-2^2", "2^3^2", "(t1+t2)/(t1*t2)", "-t1^2*t2 + 3",
              "1 - 2 - 3", "(1-2)/(3-5)", "x^4 - 3*x^2 + 1", "-(t1 - t2)^3/(t1*t2)^2"]
    ok = True
    for a in faelle:
        soll = F(eval(a.replace("^", "**"), {"__builtins__": {}}, dict(umg)))
        ist = evaluiere(a, umg)
        gut = ist == soll
        ok &= gut
        print(("  ✔ " if gut else "  🔴 ") + "%-26s %s" % (a, ist if gut else "%s != %s" % (ist, soll)))
    # (2) the case this is about: 60,000 summands — `eval` dies on the compiler stack
    gross = " + ".join("%d*t1^%d*t2^%d" % (i % 7 + 1, i % 5, i % 3) for i in range(60000))
    try:
        eval(gross.replace("^", "**"), {"__builtins__": {}}, dict(umg))
        print("  ⚠️  `eval` MANAGED the giant expression — the counter-case no longer applies on this "
              "Python version, test without meaning there.")
    except RecursionError:
        print("  ✔ `eval` dies on the giant expression (RecursionError) — exactly the case")
    soll = sum((F(i % 7 + 1) * umg["t1"] ** (i % 5) * umg["t2"] ** (i % 3) for i in range(60000)), F(0))
    ist = evaluiere(gross, umg)
    gut = ist == soll
    ok &= gut
    print(("  ✔ " if gut else "  🔴 ") + "Giant expression (60,000 summands) agrees with "
          "the reference sum")
    print("=== gross_eval SELFTEST " + ("GREEN" if ok else "RED") + " ===")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    import sys
    if "--selbsttest" in sys.argv:
        _selbsttest()
    print(__doc__)
