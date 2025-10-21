# fpgen — Generate Fokker–Planck terms (Cahill–Glauber s‑ordered)

fpgen generates the corresponding term(s) of the phase‑space Fokker–Planck equation from a quantum Liouvillian operator expression, in the Cahill–Glauber s‑parameterized quasi‑probability representations (Wigner: s=0, P: s=+1, Q: s=−1). It accepts restricted LaTeX‑like operator expressions and produces either tuple form or a readable LaTeX‑like form.

Authors: Yu Xue-Hao and Qiao Cong-Feng (University of Chinese Academy of Sciences, UCAS)

- Supported representations:
  - Wigner (s = 0)
  - P / Glauber–Sudarshan (s = +1)
  - Q / Husimi (s = −1)
- Multi-mode is supported; each latin letter (a, b, c, ...) denotes an independent mode.
- Input must contain exactly one `\rho`.

## Substitution rules (single mode)

For a general s-parameter (−1 ≤ s ≤ 1), the left/right operator products map to:

- Left multiplication
  - a ρ → α + (1 − s)/2 ∂/∂α*
  - a† ρ → α* − (1 + s)/2 ∂/∂α
- Right multiplication
  - ρ a → α − (1 + s)/2 ∂/∂α*
  - ρ a† → α* + (1 − s)/2 ∂/∂α

Special cases:
- Wigner (s = 0)
  - a ρ → α + 1/2 ∂/∂α*,  a† ρ → α* − 1/2 ∂/∂α
  - ρ a → α − 1/2 ∂/∂α*,  ρ a† → α* + 1/2 ∂/∂α
- P / Glauber–Sudarshan (s = +1)
  - a ρ → α,               a† ρ → α* − ∂/∂α
  - ρ a → α − ∂/∂α*,       ρ a† → α*
- Q / Husimi (s = −1)
  - a ρ → α + ∂/∂α*,       a† ρ → α*
  - ρ a → α,               ρ a† → α* + ∂/∂α

## API

Two main entry points:

### Main (readable + tuples):
opr2fp(expr, output_format='flux'|'normal', max_order=-1, output_style='both'|'tuple'|'latex', representation=None, s=None)

- Returns a formatted string that includes a header (representation, mode mapping, exponent headers) and either/both:
  - Tuples (when output_style includes 'tuple')
  - LaTeX-like terms (when output_style includes 'latex')
- Representation selection:
  - By s in [−1, 1], or by representation name ('Wigner', 'P' (aliases: 'Glauber', 'Glauber P', 'Glauber-Sudarshan'), 'Q' (aliases: 'Husimi', 'Husimi Q')).
  - If both are provided, s takes precedence and representation is ignored with a warning.

### Programmatic tuples only:
- _to_fp_with_representation(expr, output_format='flux'|'normal', max_order=-1, representation=None, s=None) -> List[tuple]
- to_wigner_fp_normal_form(expr, output_format='flux'|'normal', max_order=-1) -> List[tuple] (legacy, Wigner by default)

## Readable printing

- Use opr2fp(...) as the top-level convenience function.
- The header shows the selected representation and the per-mode exponent headers in the order [∂_α, ∂_{α*}, α, α*].
- Modes are mapped to Greek letters α, β, γ, ... for readability.

## Examples

```python
from fpgen import opr2fp, to_fp_normal_form

expr = "a^{\\dagger} \\rho a"

# Wigner (default)
print(opr2fp(expr))

# Explicit P representation
print(opr2fp(expr, representation='P'))

# Explicit Q representation
print(opr2fp(expr, representation='Husimi Q'))

# Using s parameter (takes precedence over representation)
print(opr2fp(expr, s=1.0, representation='Q'))

# Programmatic tuples only, P representation
tuples_p = to_fp_normal_form(expr, output_format='flux', representation='p')
print(tuples_p)
```

## Notes

- Sorting is lexicographic by exponent tuples (ignoring coefficients).
- No integration by parts is performed in flux conversion; constants and pure multipliers are preserved.
- The parser accepts LaTeX-like input with a single `\\rho` and supports powers (e.g., `a^{2}`, `(a^{\\dagger})^{2}`).

## Acknowledgements and License

This project was developed with the assistance of Copilot. Licensed under the MIT License. See LICENSE for details.
