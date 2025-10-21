"""
fpgen: s-ordered (Wigner/P/Q) Fokker–Planck toolkit.

Stable public entrypoints:
- opr2fp: Readable, latex-like presentation of the FP operator.
- to_fp_normal_form: Programmatic tuple output (normal or flux order).

All other modules in this package (parsing, substitutions, algebra, builders,
conversion, formatting, errors) are internal implementation details and may
change. Prefer importing only the two entrypoints above from this package.
"""

__version__ = "0.1.0"

from .api import (
    opr2fp,
    to_fp_normal_form,
)

__all__ = [
    "opr2fp",
    "to_fp_normal_form",
]
