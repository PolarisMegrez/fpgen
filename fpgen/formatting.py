from __future__ import annotations
"""
Readable (latex-like) formatting helpers for presenting FP terms.

These helpers are UI-only and do not affect the algebra. They are used by the
API to render human-friendly terms in addition to raw tuples.
"""
import re
from typing import List, Tuple
from fpgen.parsing import tokenize, parse, extract_ops, _normalize_expr, _split_top_level_sum

_GREEK_SEQ = [
    'α','β','γ','δ','ε','ζ','η','θ','ι','κ','λ','μ','ν','ξ','ο','π','ρ','σ','τ','υ','φ','χ','ψ','ω'
]

__all__ = [
    "_GREEK_SEQ",
    "_get_modes_in_expr",
    "_fmt_power",
    "_fmt_deriv",
    "_format_term_blocks_readable",
]

def _get_modes_in_expr(expr_latex: str) -> List[str]:
    norm = _normalize_expr(expr_latex)
    parts = _split_top_level_sum(norm)
    modes: set[str] = set()
    for _sgn, term in parts:
        # Strip optional leading numeric coefficient before tokenizing
        term_str = term.lstrip()
        mcoeff = re.match(r"^(\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?\s*(?:\*)?\s*", term_str)
        if mcoeff:
            term_str = term_str[mcoeff.end():]
        toks = tokenize(term_str)
        ast = parse(toks)
        ops = extract_ops(ast)
        for op in ops:
            modes.add(op.mode)
    return sorted(modes)

def _fmt_power(symbol: str, p: int) -> str:
    if p <= 0:
        return ''
    if p == 1:
        return symbol
    return f"{symbol}^{{{p}}}"

def _fmt_deriv(var: str, p: int) -> str:
    if p <= 0:
        return ''
    if p == 1:
        return f"∂_{{{var}}}"
    return f"∂_{{{var}}}^{{{p}}}"

def _format_term_blocks_readable(coeff: float, blocks: List[Tuple[int,int,int,int]], greek_vars: List[str]) -> str:
    pieces: List[str] = []
    # Build per-mode factor string: ∂_{g}^m ∂_{g*}^n g^i (g*)^j
    for (ex, g) in zip(blocks, greek_vars):
        m, n, i, j = ex
        segs: List[str] = []
        s1 = _fmt_deriv(g, m)
        if s1:
            segs.append(s1)
        s2 = _fmt_deriv(f"{g}*", n)
        if s2:
            segs.append(s2)
        s3 = _fmt_power(g, i)
        if s3:
            segs.append(s3)
        s4 = _fmt_power(f"{g}*", j)
        if s4:
            segs.append(s4)
        if segs:
            pieces.append(' '.join(segs))
    coeff_str = f"{coeff:g}"
    if not pieces:
        return coeff_str
    return coeff_str + " · " + " · ".join(pieces)