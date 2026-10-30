"""Codes beyond Steane, as ftlib.Code objects. Data qubits are numbered from 0.

COLOR_17: the [[17,1,5]] colour code of Liou & Lai 2024 (arXiv:2407.00607, Eq. (1), p.5);
their qubit k is our qubit k-1. Self-dual: the X and Z generators have the same supports.
Their Eq. (2) gives X_17, Z_17 as logicals, but X_17 anticommutes with g7 = Z8 Z12 Z16 Z17, so
we use the weight-5 logical on their qubits {1, 2, 5, 9, 13} found by css_distance_check
(test_color17.py, 2026-10-29).
"""
import itertools

import numpy as np

from .ftlib import Code

_LIOU_17 = [(1, 2, 3, 4), (1, 3, 5, 6), (5, 6, 9, 10), (7, 8, 11, 12), (9, 10, 13, 14),
            (11, 12, 15, 16), (8, 12, 16, 17), (3, 4, 6, 7, 10, 11, 14, 15)]
COLOR_17_SUPPORTS = [tuple(q - 1 for q in s) for s in _LIOU_17]
COLOR_17 = Code("[[17,1,5]] colour code (Liou 2024)", 17, COLOR_17_SUPPORTS, COLOR_17_SUPPORTS,
                logical_x=[(0, 1, 4, 8, 12)], logical_z=[(0, 1, 4, 8, 12)])


def _matrix(supports, n):
    m = np.zeros((len(supports), n), dtype=np.uint8)
    for i, s in enumerate(supports):
        m[i, list(s)] = 1
    return m


def gf2_rank(m):
    m = m.copy() % 2
    r = 0
    for c in range(m.shape[1]):
        rows = np.nonzero(m[r:, c])[0]
        if len(rows) == 0:
            continue
        p = r + rows[0]
        m[[r, p]] = m[[p, r]]
        for k in range(m.shape[0]):
            if k != r and m[k, c]:
                m[k] ^= m[r]
        r += 1
        if r == m.shape[0]:
            break
    return r


def css_distance_check(code, max_weight):
    """Lowest-weight nontrivial logical operator of each type, searched up to max_weight.

    For a CSS code a Pauli with X part x and Z part z is an undetectable logical error iff
    Hz x = 0, Hx z = 0 and (x not in rowspace Hx or z not in rowspace Hz); its weight is at
    least max(|x|, |z|), so the minimum over Paulis is the minimum over X-type and Z-type
    vectors separately. Returns {"X": (weight, support) or None, "Z": ...}.
    """
    n = code.n_data
    hx, hz = _matrix(code.x_supports, n), _matrix(code.z_supports, n)
    rx, rz = gf2_rank(hx), gf2_rank(hz)
    out = {}
    for kind, check, stab, rk in (("X", hz, hx, rx), ("Z", hx, hz, rz)):
        out[kind] = None
        for w in range(1, max_weight + 1):
            for sup in itertools.combinations(range(n), w):
                v = np.zeros(n, dtype=np.uint8)
                v[list(sup)] = 1
                if (check @ v % 2).any():
                    continue
                if gf2_rank(np.vstack([stab, v])) > rk:  # not a stabilizer: a logical
                    out[kind] = (w, sup)
                    break
            if out[kind]:
                break
    return out
