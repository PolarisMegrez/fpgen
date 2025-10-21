from __future__ import annotations
"""
Algebraic data structures and operations for unimode/multimode polynomials.

This module defines the core data types and algebra needed by the s-ordered
Fokker–Planck pipeline:
- Elem: elementary operator choice at a given position (α, α*, ∂α, ∂α* or identity)
- UniMono, UniPoly: unimode monomial/polynomial with normal-order exponent tuples (m,n,i,j)
- MultiMono, MultiPoly: cartesian product across modes

Notes:
- All coefficients are real floats in the public API.
"""
from dataclasses import dataclass
from typing import Dict, List, Tuple
import math

__all__ = [
    "Elem",
    "UniMono",
    "UniPoly",
    "MultiMono",
    "MultiPoly",
]

# Elementary operators for a single mode
@dataclass(frozen=True)
class Elem:
    # Each element corresponds to one of: alpha, alpha*, d/d alpha, d/d alpha*, or identity.
    mode: str
    dalpha: int = 0
    dalpha_star: int = 0
    alpha_pow: int = 0
    alpha_star_pow: int = 0
    coeff: complex = 1.0

    def to_unimono(self) -> Tuple[float, Tuple[int, int, int, int]]:
        """Convert this Elem into a UniMono tuple (coeff, (m,n,i,j))."""
        return (float(self.coeff), (self.dalpha, self.dalpha_star, self.alpha_pow, self.alpha_star_pow))

@dataclass(frozen=True)
class UniMono:
    coeff: float
    ex: Tuple[int, int, int, int]  # (m, n, i, j)

    def multiply(self, other: 'UniMono') -> 'UniPoly':
        parts = UniMono._multiply_exponents(self.ex, other.ex)
        terms = [UniMono(self.coeff * other.coeff * float(coef), ex) for (coef, ex) in parts]
        return UniPoly(terms).normalize()

    def to_tuple(self) -> Tuple[float, Tuple[int,int,int,int]]:
        return (float(self.coeff), self.ex)

    @staticmethod
    def _multiply_exponents(exA: Tuple[int, int, int, int], exB: Tuple[int, int, int, int]) -> List[Tuple[float, Tuple[int, int, int, int]]]:
        """Multiply two unimode normal-order blocks: exA then exB, and return the normal-ordered expansion (derivatives on the RIGHT).

        We need to move the multipliers of the right factor across the derivatives of the left factor using:
          d^m a^i = sum_{k=0..min(m,i)} C(m,k) (i)_k a^{i-k} d^{m-k}
        and similarly for the starred variables. No (-1)^k signs here (normal-order direction).
        """
        mA, nA, iA, jA = exA
        mB, nB, iB, jB = exB
        out: List[Tuple[float, Tuple[int, int, int, int]]] = []
        max_k = min(mA, iB)
        max_l = min(nA, jB)
        for k in range(max_k + 1):
            coef_k = math.comb(mA, k) * (math.factorial(iB) // math.factorial(iB - k) if k <= iB else 0)
            iB_left = iB - k
            mA_right = mA - k
            for l in range(max_l + 1):
                coef_l = math.comb(nA, l) * (math.factorial(jB) // math.factorial(jB - l) if l <= jB else 0)
                jB_left = jB - l
                nA_right = nA - l
                m_tot = mA_right + mB
                n_tot = nA_right + nB
                i_tot = iA + iB_left
                j_tot = jA + jB_left
                out.append((float(coef_k * coef_l), (m_tot, n_tot, i_tot, j_tot)))
        return out

@dataclass
class UniPoly:
    terms: List[UniMono]

    def multiply(self, other: 'UniPoly') -> 'UniPoly':
        out: List[UniMono] = []
        for a in self.terms:
            for b in other.terms:
                out.extend(a.multiply(b).terms)
        return UniPoly(out).normalize()

    def normalize(self) -> 'UniPoly':
        acc: Dict[Tuple[int,int,int,int], float] = {}
        for t in self.terms:
            acc[t.ex] = acc.get(t.ex, 0.0) + float(t.coeff)
        items = [(c, ex) for ex, c in acc.items() if abs(c) > 1e-12]
        items.sort(key=lambda ce: ce[1])
        return UniPoly([UniMono(float(c), ex) for (c, ex) in items])

    def sort_key(self) -> Tuple:
        return tuple(t.ex for t in self.terms)

    @staticmethod
    def from_elem_choices(choices: List[Elem]) -> 'UniPoly':
        return UniPoly([UniMono(*e.to_unimono()) for e in choices])

    def to_legacy(self) -> List[Tuple[float, Tuple[int,int,int,int]]]:
        return [t.to_tuple() for t in self.terms]

    @staticmethod
    def multiply_sequence(seq: List['UniPoly']) -> 'UniPoly':
        """Multiply a sequence of UniPoly in order.

        This provides a class-centralized implementation so higher layers don't
        need a separate helper function.
        """
        if not seq:
            return UniPoly([UniMono(1.0, (0, 0, 0, 0))])
        acc = seq[0]
        for p in seq[1:]:
            acc = acc.multiply(p)
        return acc

@dataclass(frozen=True)
class MultiMono:
    coeff: float
    blocks: Tuple[Tuple[int,int,int,int], ...]

    def to_tuple(self) -> Tuple:
        return (float(self.coeff),) + self.blocks

@dataclass
class MultiPoly:
    terms: List[MultiMono]

    def normalize(self) -> 'MultiPoly':
        acc: Dict[Tuple[Tuple[int,int,int,int], ...], float] = {}
        for t in self.terms:
            key = t.blocks
            acc[key] = acc.get(key, 0.0) + float(t.coeff)
        items = [ (c, key) for key, c in acc.items() if abs(c) > 1e-12 ]
        items.sort(key=lambda ck: ck[1])
        return MultiPoly([MultiMono(float(c), key) for (c, key) in items])

    def sort(self) -> 'MultiPoly':
        # Strict lexicographic ascending order by exponent tuples (m,n,i,j) per mode, ignoring coefficients
        return MultiPoly(sorted(self.terms, key=lambda mm: mm.blocks))

    def scale(self, s: float) -> 'MultiPoly':
        return MultiPoly([MultiMono(t.coeff * s, t.blocks) for t in self.terms])

    def truncate(self, max_order: int) -> 'MultiPoly':
        if max_order is None or max_order < 0:
            return self
        kept = []
        for t in self.terms:
            per = t.blocks
            tot = sum(m+n for (m,n,i,j) in per)
            if tot <= max_order:
                kept.append(t)
        return MultiPoly(kept).normalize()

    @staticmethod
    def from_per_mode_unipolys(mode_polys: Dict[str, UniPoly]) -> 'MultiPoly':
        modes_sorted = sorted(mode_polys.keys())
        if not modes_sorted:
            return MultiPoly([MultiMono(1.0, ((0,0,0,0),))])
        acc: List[MultiMono] = [MultiMono(1.0, tuple())]
        for m in modes_sorted:
            poly = mode_polys[m]
            new_acc: List[MultiMono] = []
            for mm in acc:
                for um in poly.terms:
                    new_acc.append(MultiMono(mm.coeff * um.coeff, mm.blocks + (um.ex,)))
            acc = new_acc
        return MultiPoly(acc).normalize()

    def to_legacy(self) -> List[Tuple]:
        return [t.to_tuple() for t in self.terms]