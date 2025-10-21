from __future__ import annotations
"""
Builders: glue from parsing/substitutions to unimode/multimode algebra.

This module exposes helpers to translate parsed operator expressions into
unimode sequences and finally into multi-mode polynomials. Substitutions are
+narrowly scoped and default to the Wigner (s=0) map via SUBS; higher-level API
+should pass explicit s/representation when needed.
"""
import re
from typing import List, Tuple, Dict
from fpgen.parsing import Op, tokenize, parse, extract_ops, _normalize_expr, _split_top_level_sum
from fpgen.substitutions import SUBS
from fpgen.errors import ParseError
from fpgen.algebra import Elem, UniPoly

__all__ = [
    "expr_to_elem_sequence",
    "expr_to_elem_polynomial",
    "elems_to_unimode_sequences",
    "term_ops_to_mode_sequences",
]

def expr_to_elem_sequence(expr_latex: str) -> List[List[Elem]]:
    """Convert a restricted LaTeX operator expression into a sequence of Elem-choice lists.

    Each operator around \rho expands to a small list of Elem (sum) at that position.
    The output is a list matching the product order: [ [Elem,...], [Elem,...], ... ].
    """
    norm = _normalize_expr(expr_latex)
    parts = _split_top_level_sum(norm)
    if len(parts) != 1:
        # Only single term supported by this helper; sums handled at a higher level
        raise ParseError('expr_to_elem_sequence expects a single term (no top-level +/-)')
    _, term = parts[0]
    toks = tokenize(term)
    ast = parse(toks)
    ops = extract_ops(ast)
    seq: List[List[Elem]] = []
    for op in ops:
        fn = SUBS.get((op.side, op.kind))
        if not fn:
            raise ParseError('Internal: unknown substitution')
        # Op.power already expanded by tokenizer for grouped powers; keep power for safety
        for _ in range(op.power):
            seq.append(fn(op.mode))
    return seq

def expr_to_elem_polynomial(expr_latex: str) -> List[Tuple[float, List[List[Elem]]]]:
    """Parse expression with top-level +/- into a list of (scalar, Elem-sequence).

    Each entry is a scalar coefficient and a sequence of Elem-choice lists for that term.
    """
    norm = _normalize_expr(expr_latex)
    parts = _split_top_level_sum(norm)
    out: List[Tuple[float, List[List[Elem]]]] = []
    for sgn, term in parts:
        coeff_scalar = 1.0
        term_str = term.lstrip()
        mcoeff = re.match(r"^(\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?\s*(?:\*)?\s*", term_str)
        if mcoeff:
            coeff_scalar = float(mcoeff.group(0).replace('*', ' ').strip())
            term_str = term_str[mcoeff.end():]
        toks = tokenize(term_str)
        ast = parse(toks)
        ops = extract_ops(ast)
        seq: List[List[Elem]] = []
        for op in ops:
            fn = SUBS.get((op.side, op.kind))
            if not fn:
                raise ParseError('Internal: unknown substitution')
            for _ in range(op.power):
                seq.append(fn(op.mode))
        out.append((float(sgn)*coeff_scalar, seq))
    return out

def elems_to_unimode_sequences(seq: List[List[Elem]]) -> Dict[str, List['UniPoly']]:
    """Split a sequence of Elem-choice lists into per-mode UniPoly sequences.

    Returns a dict: mode -> [UniPoly_at_pos1, UniPoly_at_pos2, ...].
    """
    by_mode: Dict[str, List[UniPoly]] = {}
    for choices in seq:
        # Each position contributes a small UniPoly (sum over its Elem choices)
        poly: UniPoly = UniPoly.from_elem_choices(choices)
        # All Elems in this choices list share the same mode by construction
        if not choices:
            continue
        mode = choices[0].mode
        by_mode.setdefault(mode, []).append(poly)
    return by_mode

def term_ops_to_mode_sequences(ops: List['Op']) -> Dict[str, List['UniPoly']]:
    """For a given term (sequence of ops around rho), build per-mode sequence of UniPoly.
    Each op contributes a small UniPoly (sum of up to 2 UniMonos) appended in multiplication order.
    """
    by_mode: Dict[str, List[UniPoly]] = {}
    for op in ops:
        fn = SUBS.get((op.side, op.kind))
        if not fn:
            raise ParseError('Internal: unknown substitution')
        elems = fn(op.mode)
        # Convert Elem list into a UniPoly (sum) for this position
        poly: UniPoly = UniPoly.from_elem_choices(elems)
        lst = by_mode.setdefault(op.mode, [])
        # For powers, Op already repeats; so just append poly once per Op instance
        lst.append(poly)
    return by_mode