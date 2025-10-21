from __future__ import annotations
"""
Parsing utilities for restricted LaTeX-like operator expressions.

Responsibilities:
- Tokenization (Token)
- Minimal AST with parentheses (Node)
- Operator extraction around a single \rho (Op)
- Expression normalization and top-level +/- splitting

Errors:
- ParseError on structural or token issues.
"""
import re
from typing import List, Tuple, Optional
from dataclasses import dataclass
from fpgen.errors import ParseError

TOKEN_RE = re.compile(r"\\rho|[a-zA-Z](?:\^{\\dagger|\+|\*})?|\^\{\d+}|[()]")
_GROUP_CREATOR_PWR = re.compile(r"\(([a-zA-Z])\^\{(\\dagger|\+|\*)\}\)\^\{(\d+)\}")
_GROUP_ANN_PWR = re.compile(r"\(([a-zA-Z])\)\^\{(\d+)\}")

__all__ = [
    "TOKEN_RE",
    "ParseError",
    "Token",
    "Node",
    "tokenize",
    "parse",
    "extract_ops",
    "_normalize_expr",
    "_split_top_level_sum",
]

# Tokenize LaTeX-like input

@dataclass
class Token:
    kind: str
    text: str
    mode: Optional[str] = None
    power: int = 1

def tokenize(expr: str) -> List[Token]:
    # Note: top-level '+'/'-' are split before tokenization; '*' multiplication is normalized out.
    # We avoid pre-validating characters too strictly; the parser below will raise on unsupported tokens.
    tokens: List[Token] = []
    i = 0
    while i < len(expr):
        ch = expr[i]
        if ch.isspace():
            i += 1
            continue
        if expr.startswith('\\rho', i):
            tokens.append(Token('rho', '\\rho'))
            i += 4
            continue
        # mode letters
        if ch.isalpha():
            mode = ch
            i += 1
            # possible superscript dagger/+/* using ^{...}
            if i < len(expr) and expr[i] == '^':
                # expect ^{...}
                if i+1 < len(expr) and expr[i+1] == '{':
                    j = expr.find('}', i+2)
                    if j == -1:
                        raise ParseError('Missing closing } for superscript')
                    sup = expr[i+2:j]
                    # Support a^{\dagger n} by splitting into creator with implicit power n
                    m = re.match(r"^(\\dagger|\+|\*)\s*(\d+)$", sup)
                    if sup in ('\\dagger', '+', '*'):
                        tokens.append(Token('creator', f"{mode}^{{{sup}}}", mode=mode))
                        i = j+1
                    elif m:
                        # creator with power, we will record as repeated later by emitting 'creator' followed by (n-1) more creators in parsing stage
                        # Here we store as a special ann_pwr-like token carrying power for creator
                        power = int(m.group(2))
                        if power <= 0:
                            raise ParseError('Non-positive power not allowed')
                        # We'll encode this as multiple creator atoms in downstream extraction by expanding now to repeated creators
                        # Expand inline for simplicity
                        # Replace the current token stream position with creator repeated 'power' times by pushing tokens
                        for _ in range(power):
                            tokens.append(Token('creator', f"{mode}^{{\\dagger}}", mode=mode))
                        i = j+1
                    else:
                        # power on annihilation, e.g., a^{2} -> expand to repeated 'ann'
                        if sup.isdigit():
                            power = int(sup)
                            if power <= 0:
                                raise ParseError('Non-positive power not allowed')
                            for _ in range(power):
                                tokens.append(Token('ann', mode, mode=mode))
                            i = j+1
                        else:
                            raise ParseError('Invalid superscript after ^{...}')
                else:
                    raise ParseError('Expected ^{...} after mode letter')
            else:
                tokens.append(Token('ann', mode, mode=mode))
            continue
        if ch in '()':
            tokens.append(Token('paren', ch))
            i += 1
            continue
        # stray symbols
        raise ParseError('Unexpected character')
    return tokens

# AST nodes (very small because we only need concatenation and parentheses)

@dataclass
class Node:
    kind: str  # 'seq' or 'atom'
    value: Optional[Token] = None
    children: Optional[List['Node']] = None

def parse(tokens: List[Token]) -> Node:
    # Recursive descent for parentheses and sequences
    def parse_seq(i: int) -> Tuple[Node, int]:
        items: List[Node] = []
        while i < len(tokens):
            t = tokens[i]
            if t.kind == 'paren' and t.text == ')':
                break
            if t.kind == 'paren' and t.text == '(':
                sub, j = parse_seq(i+1)
                if j >= len(tokens) or tokens[j].kind != 'paren' or tokens[j].text != ')':
                    raise ParseError('Missing closing )')
                items.append(sub)
                i = j+1
                continue
            items.append(Node('atom', value=t))
            i += 1
        return Node('seq', children=items), i
    root, idx = parse_seq(0)
    if idx != len(tokens):
        if tokens[idx].kind == 'paren' and tokens[idx].text == ')':
            raise ParseError('Unmatched )')
        raise ParseError('Unexpected tokens at end')
    return root

# Build operator choices list, marking positions w.r.t. rho

@dataclass
class Op:
    side: str  # 'L' or 'R' relative to rho
    kind: str  # 'a' or 'adag'
    mode: str
    power: int = 1

def extract_ops(root: Node) -> List[Op]:
    # Flatten sequence, find rho, everything before is left-ops (acts on rho from left), after is right-ops
    flat: List[Token] = []
    def flatten(n: Node):
        if n.kind == 'atom':
            flat.append(n.value)  # type: ignore[arg-type]
        else:
            for c in (n.children or []):
                flatten(c)
    flatten(root)
    if not any(t.kind == 'rho' for t in flat):
        raise ParseError('Expression must contain \\rho exactly once')
    # Ensure only one rho
    if sum(1 for t in flat if t.kind == 'rho') != 1:
        raise ParseError('Expression must contain exactly one \\rho')
    rho_idx = next(i for i,t in enumerate(flat) if t.kind == 'rho')
    left = flat[:rho_idx]
    right = flat[rho_idx+1:]
    ops: List[Op] = []
    # Parse left side (in order): a, adag, and powers
    i = 0
    while i < len(left):
        t = left[i]
        if t.kind == 'ann':
            ops.append(Op('L', 'a', t.mode or '?', 1))
        elif t.kind == 'ann_pwr':
            ops.append(Op('L', 'a', t.mode or '?', t.power))
        elif t.kind == 'creator':
            ops.append(Op('L', 'adag', t.mode or '?', 1))
        else:
            raise ParseError('Unsupported token on left of \\rho')
        i += 1
    # Parse right side (in order)
    i = 0
    while i < len(right):
        t = right[i]
        if t.kind == 'ann':
            ops.append(Op('R', 'a', t.mode or '?', 1))
        elif t.kind == 'ann_pwr':
            ops.append(Op('R', 'a', t.mode or '?', t.power))
        elif t.kind == 'creator':
            ops.append(Op('R', 'adag', t.mode or '?', 1))
        else:
            raise ParseError('Unsupported token on right of \\rho')
        i += 1
    return ops

# Expression normalization: handle explicit '*' multiplication, sums, and parenthesized powers like (a^{\dagger})^{2}

def _expand_group_powers(s: str) -> str:
    # Expand (a^{dagger})^{n} -> a^{dagger} a^{dagger} ... (n times)
    while True:
        m = _GROUP_CREATOR_PWR.search(s)
        if not m:
            break
        letter, sup, pwr = m.group(1), m.group(2), int(m.group(3))
        if pwr <= 0:
            raise ParseError('Non-positive power not allowed')
        repl = ' '.join([f"{letter}^{{{sup}}}" for _ in range(pwr)])
        s = s[:m.start()] + repl + s[m.end():]
    while True:
        m = _GROUP_ANN_PWR.search(s)
        if not m:
            break
        letter, pwr = m.group(1), int(m.group(2))
        if pwr <= 0:
            raise ParseError('Non-positive power not allowed')
        repl = ' '.join([letter for _ in range(pwr)])
        s = s[:m.start()] + repl + s[m.end():]
    return s

def _normalize_expr(expr: str) -> str:
    # Remove explicit multiplication signs and normalize grouped powers
    s = expr.replace('*', ' ')
    s = _expand_group_powers(s)
    return s

def _split_top_level_sum(expr: str) -> List[Tuple[int, str]]:
    # Split expr into terms by top-level +/-. Return list of (sign, term)
    terms: List[Tuple[int,str]] = []
    depth_paren = 0
    depth_brace = 0
    start = 0
    sign = 1
    i = 0
    # Trim leading spaces
    expr = expr.strip()
    # Handle leading sign
    if expr.startswith('+'):
        expr = expr[1:].lstrip()
    elif expr.startswith('-'):
        sign = -1
        expr = expr[1:].lstrip()
    i = 0
    while i < len(expr):
        ch = expr[i]
        if ch == '{':
            depth_brace += 1
        elif ch == '}':
            depth_brace = max(0, depth_brace-1)
        elif ch == '(':
            depth_paren += 1
        elif ch == ')':
            depth_paren = max(0, depth_paren-1)
        elif ch in ['+','-'] and depth_paren == 0 and depth_brace == 0:
            term = expr[:i].rstrip()
            if term:
                terms.append((sign, term))
            # Prepare for next term
            sign = 1 if ch == '+' else -1
            expr = expr[i+1:].lstrip()
            i = 0
            continue
        i += 1
    # Remainder
    term = expr.strip()
    if term:
        terms.append((sign, term))
    return terms