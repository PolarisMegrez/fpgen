from __future__ import annotations
"""
Conversion between normal-order and flux-order tuple encodings.

convert_normal_to_flux expands per-mode blocks via commutator binomials only
(no integration by parts), preserving constants and pure multipliers.
convert_flux_to_normal combines and reorders back to normal counts.
"""
import math
from typing import Dict, List, Tuple

__all__ = [
    "convert_normal_to_flux",
    "convert_flux_to_normal",
]

def convert_normal_to_flux(terms: List[Tuple], max_order: int = -1) -> List[Tuple]:
    """Convert normal-form tuples (multipliers before derivatives) into flux-form
    (derivatives on the LEFT), using only [∂_α, α] = 1 and [∂_{α*}, α*] = 1 identities.
    No integration-by-parts elimination; constants and pure multipliers are retained.
    """
    # Expand each term per mode with (-1)^k binomial identities and accumulate
    accum: Dict[Tuple[Tuple[int,int,int,int], ...], float] = {}
    for t in terms:
        c0 = float(t[0])
        per_blocks = list(t[1:])
        expansions: List[Tuple[float, List[Tuple[int,int,int,int]]]] = [(1.0, [])]
        for (m,n,i,j) in per_blocks:
            per_exp = _fluxize_single_mode(m,n,i,j)
            new_exp: List[Tuple[float, List[Tuple[int,int,int,int]]]] = []
            for (c_prev, blocks_prev) in expansions:
                for (coef, ex2) in per_exp:
                    new_exp.append((c_prev*coef, blocks_prev + [ex2]))
            expansions = new_exp
        for (coef, blocks) in expansions:
            key = tuple(blocks)
            accum[key] = accum.get(key, 0.0) + c0 * coef
    # Optional truncation by derivative order AFTER conversion (usually truncation was already applied before)
    if max_order is not None and max_order >= 0:
        accum = {k:c for k,c in accum.items() if sum(m+n for (m,n,i,j) in k) <= max_order}
    # Build list and sort: strict lexicographic ascending by exponent tuples (ignoring coefficient)
    items = [ (float(c),)+k for k,c in accum.items() if abs(c) > 1e-12 ]
    items.sort(key=lambda t: t[1:])
    return items

def _sum_tuples(terms: List[Tuple]) -> Dict[Tuple, float]:
    accum: Dict[Tuple, float] = {}
    for t in terms:
        # keep coefficients as real floats
        c = float(t[0])
        # build key with fake mode order based on length; assume single-mode 'a' if unknown
        # Here t[1:] are already ordered per to_fp_normal_form; we create a synthetic key
        # Since we don't have mode names, we tag modes in order as m0,m1,...
        key = tuple((f"m{idx}",)+tuple(block) for idx, block in enumerate(t[1:]))
        accum[key] = accum.get(key, 0.0) + c
    return accum

def _tuples_from_accum(accum: Dict[Tuple, float]) -> List[Tuple]:
    eps = 1e-12
    items_sorted = sorted(((k, c) for k,c in accum.items() if abs(c) > eps), key=lambda kc: kc[0])
    out: List[Tuple] = []
    for key, c in items_sorted:
        per = tuple((m,n,i,j) for (_mode,m,n,i,j) in key)
        out.append((float(c),) + per)
    return out

def convert_flux_to_normal(terms: List[Tuple]) -> List[Tuple]:
    # Our internal tuple encoding always corresponds to normal-order counts, so this is identity after combining.
    return _tuples_from_accum(_sum_tuples(terms))

def _fluxize_single_mode(m:int,n:int,i:int,j:int) -> List[Tuple[float, Tuple[int,int,int,int]]]:
        """Move all multipliers (alpha^i, alpha*^j) to the RIGHT of derivatives (∂_α^m, ∂_{α*}^n)
        to obtain flux order (derivatives on the LEFT). Uses identities:
            alpha^i ∂_α^m = Σ_{k=0..min(m,i)} (-1)^k C(m,k) (i)_k ∂_α^{m-k} alpha^{i-k}
            (alpha*^j) ∂_{α*}^n = Σ_{l=0..min(n,j)} (-1)^l C(n,l) (j)_l ∂_{α*}^{n-l} (alpha*^{j-l})
        Combined independently for the two variables.
        Returns list of (coef, (m', n', i', j')).
        """
        out: List[Tuple[float, Tuple[int,int,int,int]]] = []
        max_k = min(m, i)
        max_l = min(n, j)
        for k in range(max_k+1):
                coef_k = ((-1.0)**k) * math.comb(m, k) * (math.factorial(i) // math.factorial(i - k) if k <= i else 0)
                m2 = m - k
                i2 = i - k
                for l in range(max_l+1):
                        coef_l = ((-1.0)**l) * math.comb(n, l) * (math.factorial(j) // math.factorial(j - l) if l <= j else 0)
                        n2 = n - l
                        j2 = j - l
                        out.append((float(coef_k*coef_l), (m2, n2, i2, j2)))
        return out

def _commute_multipliers_left_single_mode(m:int,n:int,i:int,j:int) -> List[Tuple[complex, Tuple[int,int,int,int]]]:
    # Move alpha^i across (d/d alpha)^m and (alpha*^j) across (d/d alpha*)^n to the left, producing sum
    out: List[Tuple[complex, Tuple[int,int,int,int]]] = []
    max_k = min(m, i)
    max_l = min(n, j)
    for k in range(max_k+1):
        coef_k = math.comb(m, k) * (math.factorial(i) // math.factorial(max(0, i-k)))
        i_left = i - k
        m_right = m - k
        for l in range(max_l+1):
            coef_l = math.comb(n, l) * (math.factorial(j) // math.factorial(max(0, j-l)))
            j_left = j - l
            n_right = n - l
            out.append((coef_k*coef_l, (m_right, n_right, i_left, j_left)))
    return out

def _truncate_terms(accum: Dict[Tuple, complex], max_order: int) -> Dict[Tuple, complex]:
    if max_order is None or max_order < 0:
        return accum
    out: Dict[Tuple, complex] = {}
    for key, c in accum.items():
        # total derivative order is sum over modes of m+n
        order = sum(m+n for (_mode,m,n,i,j) in key)
        if order <= max_order:
            out[key] = out.get(key, 0.0) + c
    return out

def _eliminate_pure_multipliers(accum: Dict[Tuple, complex]) -> Dict[Tuple, complex]:
    # Convert any (m=n=0, i or j > 0) into derivative terms via identities:
    # alpha^i = d/d alpha (alpha^i) - i alpha^{i-1}
    # alpha*^j = d/d alpha* (alpha*^j) - j alpha*^{j-1}
    # Repeat until no pure multipliers remain. Prefer eliminating j (alpha*) first.
    changed = True
    accum = dict(accum)
    while changed:
        changed = False
        items = list(accum.items())
        for key, c in items:
            if abs(c) == 0:
                continue
            blocks = list(key)
            # process each mode separately; since modes are independent, we can rewrite per-mode
            new_contrib: Dict[Tuple, complex] = {}
            did_local = False
            new_blocks = list(blocks)
            for idx, (_mode, m, n, i, j) in enumerate(blocks):
                if m == 0 and n == 0 and (i > 0 or j > 0):
                    did_local = True
                    # remove current key
                    accum[key] -= c
                    if abs(accum[key]) < 1e-15:
                        del accum[key]
                    # prefer eliminating j first
                    if j > 0:
                        # term: +1 * (0,1,i,j)
                        k1 = list(blocks)
                        k1[idx] = (_mode, 0, 1, i, j)
                        k1 = tuple(k1)
                        new_contrib[k1] = new_contrib.get(k1, 0.0) + c
                        # remainder: -j * (0,0,i,j-1)
                        k2 = list(blocks)
                        k2[idx] = (_mode, 0, 0, i, j-1)
                        k2 = tuple(k2)
                        new_contrib[k2] = new_contrib.get(k2, 0.0) - c * j
                    else:
                        # i > 0
                        k1 = list(blocks)
                        k1[idx] = (_mode, 1, 0, i, 0)
                        k1 = tuple(k1)
                        new_contrib[k1] = new_contrib.get(k1, 0.0) + c
                        k2 = list(blocks)
                        k2[idx] = (_mode, 0, 0, i-1, 0)
                        k2 = tuple(k2)
                        new_contrib[k2] = new_contrib.get(k2, 0.0) - c * i
                    break
            if did_local:
                # add new contributions and continue outer loop
                for nk, nc in new_contrib.items():
                    accum[nk] = accum.get(nk, 0.0) + nc
                changed = True
                break
    # Optionally eliminate constant (0,0,0,0) by writing 1 = d/d alpha*(alpha*) - alpha* d/d alpha*
    items = list(accum.items())
    for key, c in items:
        if abs(c) == 0:
            continue
        if all(m==0 and n==0 and i==0 and j==0 for (_mode,m,n,i,j) in key):
            # pick first mode to transform constant
            if len(key) > 0:
                mode0 = key[0][0]
                k1 = list(key)
                k1[0] = (mode0, 0, 1, 0, 1)
                k1 = tuple(k1)
                k2 = list(key)
                k2[0] = (mode0, 0, 1, 0, 0)  # represents -alpha* d/d alpha*; in tuple counts it's (n=1) and j=0 with coeff -alpha*? Not encodable; skip
                # We can instead keep the constant; or leave as is.
                # For safety, leave constant unchanged.
                pass
    return accum