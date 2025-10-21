from __future__ import annotations
"""
Cahill–Glauber s-ordered substitution builders and representation helpers.

Exposes:
- _make_subs_for_s(s): build left/right maps for (a, a†) around ρ.
- SUBS: default Wigner (s=0) substitutions for internal helpers.
- _normalize_representation_name/_rep_to_s/_resolve_s_and_subs: resolve s vs. named (Wigner/P/Q).

Notes:
- Higher-level APIs should prefer passing s explicitly; SUBS is a convenience default.
"""
import warnings
from typing import Dict, List, Tuple, Optional
from fpgen.errors import WignerFPError
from fpgen.algebra import Elem

# Substitution builders for general s-representation (Cahill–Glauber)
def _make_subs_for_s(s: float) -> Dict[Tuple[str,str], 'callable']:
    # Left: a ρ -> α + (1-s)/2 ∂_{α*}; a† ρ -> α* - (1+s)/2 ∂_α
    # Right: ρ a -> α - (1+s)/2 ∂_{α*}; ρ a† -> α* + (1-s)/2 ∂_α
    c_L_a = 0.5 * (1.0 - s)
    c_L_adag = -0.5 * (1.0 + s)
    c_R_a = -0.5 * (1.0 + s)
    c_R_adag = 0.5 * (1.0 - s)
    def left_a(mode: str) -> List[Elem]:
        return [Elem(mode, alpha_pow=1), Elem(mode, dalpha_star=1, coeff=c_L_a)]
    def left_adag(mode: str) -> List[Elem]:
        return [Elem(mode, alpha_star_pow=1), Elem(mode, dalpha=1, coeff=c_L_adag)]
    def right_a(mode: str) -> List[Elem]:
        return [Elem(mode, alpha_pow=1), Elem(mode, dalpha_star=1, coeff=c_R_a)]
    def right_adag(mode: str) -> List[Elem]:
        return [Elem(mode, alpha_star_pow=1), Elem(mode, dalpha=1, coeff=c_R_adag)]
    return {
        ('L','a'): left_a,
        ('L','adag'): left_adag,
        ('R','a'): right_a,
        ('R','adag'): right_adag,
    }

# Default substitutions (Wigner, s=0) for legacy helpers
SUBS = _make_subs_for_s(0.0)

def _normalize_representation_name(rep: Optional[str]) -> Optional[str]:
    if rep is None:
        return None
    r = rep.strip().lower()
    aliases_p = {'p','glauber','glauber p','glauber-sudarshan','glauber sudarshan'}
    aliases_q = {'q','husimi','husimi q'}
    aliases_w = {'w','wigner','wigner-weyl','weyl','weyl-wigner'}
    if r in aliases_p:
        return 'p'
    if r in aliases_q:
        return 'q'
    if r in aliases_w:
        return 'wigner'
    return r

def _rep_to_s(rep_norm: Optional[str]) -> float:
    if rep_norm is None or rep_norm == 'wigner':
        return 0.0
    if rep_norm == 'p':
        return 1.0
    if rep_norm == 'q':
        return -1.0
    # Unknown: default to Wigner
    return 0.0

def _resolve_s_and_subs(s: Optional[float], representation: Optional[str]) -> Tuple[float, Dict[Tuple[str,str], 'callable']]:
    rep_norm = _normalize_representation_name(representation)
    if s is not None and rep_norm is not None:
        warnings.warn("Both 's' and 'representation' provided; using 's' and ignoring 'representation'.")
    if s is None:
        s_val = _rep_to_s(rep_norm)
    else:
        s_val = float(s)
    if not (-1.0 <= s_val <= 1.0):
        raise WignerFPError("Parameter 's' must be within [-1, 1].")
    subs = _make_subs_for_s(s_val)
    return s_val, subs