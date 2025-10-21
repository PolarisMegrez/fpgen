from __future__ import annotations
"""Public API for the s-ordered (Wigner/P/Q) Fokker–Planck toolkit.

Stable entrypoints:
- opr2fp: Readable output (latex-like terms and/or tuples)
- to_fp_normal_form: Programmatic tuple output (normal or flux order)

All other modules under fpgen.* are internal and may change.
"""
from typing import Dict, List, Tuple, Optional
import re
from fpgen.parsing import tokenize, parse, extract_ops, _normalize_expr, _split_top_level_sum
from fpgen.substitutions import _resolve_s_and_subs
from fpgen.algebra import UniPoly, MultiMono, MultiPoly
from fpgen.conversion import convert_normal_to_flux, convert_flux_to_normal
from fpgen.errors import WignerFPError, ParseError
from fpgen.formatting import _GREEK_SEQ, _get_modes_in_expr, _format_term_blocks_readable

__all__ = [
    "opr2fp",
    "to_fp_normal_form",
    "convert_normal_to_flux",
    "convert_flux_to_normal",
]

def opr2fp(
    expr_latex: str,
    *,
    output_format: str = 'flux',
    max_order: int = -1,
    output_style: str = 'both',  # 'both' | 'tuple' | 'latex'
    representation: Optional[str] = None,
    s: Optional[float] = None,
) -> str:
    """Top-level user API: readable Fokker–Planck operator output for s-ordered representations.

    - output_format: 'flux' (default) or 'normal'.
    - Representation selection: by 's' in [-1,1] or by 'representation' ('Wigner'|'P'|'Q' with aliases).
      If both provided, 's' wins (warning issued). Header shows either 'Wigner/P/Q' or 'Cahill-Glauber s=...'.
    - output_style: 'both' | 'tuple' | 'latex'. Bullets use '>' to avoid confusion with minus sign.
    """
    # Determine modes (sorted) and map to Greek variables
    modes = _get_modes_in_expr(expr_latex)
    greek_vars = [_GREEK_SEQ[i] if i < len(_GREEK_SEQ) else f"ζ{i+1}" for i in range(len(modes))]
    mode_map = {m: g for m, g in zip(modes, greek_vars)}
    # Resolve representation
    s_val, _subs = _resolve_s_and_subs(s, representation)
    # Determine header name per user rule
    if s is not None:
        rep_header = f"Cahill-Glauber s={s:g}"
    else:
        rep_header = 'Wigner' if abs(s_val) < 1e-15 else ('P' if s_val > 0 else 'Q')
    # Get tuples using existing pipeline (already lexicographically sorted)
    tuples = _to_fp_with_representation(
        expr_latex,
        output_format=output_format,
        max_order=max_order,
        representation=representation,
        s=s,
    )
    # Header
    lines: List[str] = []
    lines.append(f"Representation: {rep_header}")
    if modes:
        mapping_str = ', '.join(f"{m}->{mode_map[m]}" for m in modes)
        lines.append(f"Mode mapping: {mapping_str}")
        headers = ' | '.join(
            f"[∂_{{{mode_map[m]}}}, ∂_{{{mode_map[m]}*}}, {mode_map[m]}, {mode_map[m]}*]" for m in modes
        )
        lines.append(f"Exponent headers (per-mode order): {headers}")
    else:
        lines.append("Exponent headers: [∂_{α}, ∂_{α*}, α, α*]")

    if not tuples:
        lines.append("(empty)")
        return "\n".join(lines)

    # Sections based on output_style
    if output_style not in ('both','tuple','latex'):
        output_style = 'both'

    if output_style in ('both','tuple'):
        lines.append("Tuples:")
        for t in tuples:
            lines.append("  > " + str(t))

    if output_style in ('both','latex'):
        lines.append("Terms:")
        for t in tuples:
            c = float(t[0])
            blocks = list(t[1:])
            if len(blocks) < len(greek_vars):
                blocks = blocks + [(0,0,0,0)] * (len(greek_vars) - len(blocks))
            term_str = _format_term_blocks_readable(c, blocks, greek_vars)
            lines.append("  > " + term_str)

    return "\n".join(lines)

def to_fp_normal_form(
    expr_latex: str,
    *,
    output_format: str = 'flux',
    max_order: int = -1,
    representation: Optional[str] = None,
    s: Optional[float] = None,
) -> List[Tuple]:
    """Primary programmatic API: compute FP tuples in s-ordered (Wigner/P/Q) representations.

    - Choose output_format='flux' or 'normal'.
    - Select representation either by 's' in [-1,1] or by 'representation' name (aliases supported).
      If both are provided, 's' takes precedence and a warning is issued.
    """
    return _to_fp_with_representation(
        expr_latex,
        output_format=output_format,
        max_order=max_order,
        representation=representation,
        s=s,
    )

def _to_fp_with_representation(
    expr_latex: str,
    *,
    output_format: str = 'flux',
    max_order: int = -1,
    representation: Optional[str] = None,
    s: Optional[float] = None,
) -> List[Tuple]:
    try:
        s_val, subs = _resolve_s_and_subs(s, representation)
        norm = _normalize_expr(expr_latex)
        # Support top-level sums with + and -
        parts = _split_top_level_sum(norm)
        # Build MultiPoly per term using the class-based unimode pipeline
        term_polys: List[MultiPoly] = []
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
            # per-mode UniPoly sequences with correct order:
            # - Left of rho: multiply left-to-right
            # - Right of rho: multiply right-to-left
            by_mode: Dict[str, List[UniPoly]] = {}
            left_ops = [op for op in ops if op.side == 'L']
            right_ops = [op for op in ops if op.side == 'R']
            # process left in given order
            for op in left_ops:
                fn = subs.get((op.side, op.kind))
                if not fn:
                    raise ParseError('Internal: unknown substitution')
                elems = fn(op.mode)
                poly = UniPoly.from_elem_choices(elems)
                by_mode.setdefault(op.mode, []).append(poly)
            # process right in reverse order
            for op in reversed(right_ops):
                fn = subs.get((op.side, op.kind))
                if not fn:
                    raise ParseError('Internal: unknown substitution')
                elems = fn(op.mode)
                poly = UniPoly.from_elem_choices(elems)
                by_mode.setdefault(op.mode, []).append(poly)
            # multiply sequences per mode
            mode_polys: Dict[str, UniPoly] = {m: UniPoly.multiply_sequence(seq) for m, seq in by_mode.items()}
            multi_poly = MultiPoly.from_per_mode_unipolys(mode_polys)
            # apply overall scalar and sign
            term_polys.append(multi_poly.scale(float(sgn) * coeff_scalar))
        # Sum and normalize
        if term_polys:
            all_terms: List[MultiMono] = [t for p in term_polys for t in p.terms]
            normal_poly = MultiPoly(all_terms).normalize().truncate(max_order)
        else:
            normal_poly = MultiPoly([])
        if output_format == 'normal':
            return normal_poly.to_legacy()
        # Convert normal to flux order using existing converter
        return convert_normal_to_flux(normal_poly.to_legacy(), max_order=max_order)
    except ParseError as e:
        raise WignerFPError(str(e))