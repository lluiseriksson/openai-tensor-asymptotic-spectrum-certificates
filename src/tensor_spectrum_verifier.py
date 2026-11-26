"""TENSOR-SPEC-CERT-1.0: Exact rational/integer certificate verifier for
OpenAI's 'Euclidean geometry in the asymptotic spectrum of tensors' (Nov 24, 2026).

Strictly zero floating-point arithmetic, zero ast.Div / /=, zero math imports.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any

SCHEMA_VERSION = "TENSOR-SPEC-CERT-1.0"
MAX_JSON_BYTES = 1_000_000
MAX_RATIONAL_ABS = 10**12
MAX_MAT_DIM = 32
MAX_CUTS = 16
MAX_SITES = 16
MAX_EXT_Q = 6
MAX_LIST_LEN = 64


class VerificationError(ValueError):
    """Raised when a TENSOR-SPEC-CERT-1.0 certificate fails any check."""


def _no_dup_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise VerificationError(f"Duplicate JSON key: {k!r}")
        out[k] = v
    return out


def _gcd(a: int, b: int) -> int:
    a = -a if a < 0 else a
    b = -b if b < 0 else b
    while b != 0:
        a, b = b, a % b
    return a


def _req_dict(obj: Any, name: str) -> dict[str, Any]:
    if not isinstance(obj, dict):
        raise VerificationError(f"{name} must be a JSON object")
    return obj


def _req_list(obj: Any, name: str, max_len: int = MAX_LIST_LEN) -> list[Any]:
    if not isinstance(obj, list):
        raise VerificationError(f"{name} must be a JSON array")
    if len(obj) > max_len:
        raise VerificationError(f"{name} exceeds max length {max_len}")
    return obj


def _req_str(obj: Any, name: str) -> str:
    if not isinstance(obj, str) or not obj.strip():
        raise VerificationError(f"{name} must be a non-empty string")
    return obj


def _req_int(obj: Any, name: str, lo: int | None = None, hi: int | None = None) -> int:
    if isinstance(obj, bool) or not isinstance(obj, int):
        raise VerificationError(f"{name} must be an integer")
    if lo is not None and obj < lo:
        raise VerificationError(f"{name}={obj} < {lo}")
    if hi is not None and obj > hi:
        raise VerificationError(f"{name}={obj} > {hi}")
    return obj


def _req_bool(obj: Any, name: str) -> bool:
    if not isinstance(obj, bool):
        raise VerificationError(f"{name} must be a boolean")
    return obj


def _req_frac(obj: Any, name: str) -> Fraction:
    """Parse exact lowest-terms rational [p, q] with q > 0, gcd(|p|, q) = 1."""
    if not isinstance(obj, list) or len(obj) != 2:
        raise VerificationError(f"{name} must be a [p, q] rational pair")
    p = _req_int(obj[0], f"{name}[0]", -MAX_RATIONAL_ABS, MAX_RATIONAL_ABS)
    q = _req_int(obj[1], f"{name}[1]", 1, MAX_RATIONAL_ABS)
    if _gcd(p, q) != 1:
        raise VerificationError(f"{name}=[{p},{q}] is not in lowest terms")
    return Fraction(p, q)


def _req_frac_vec(obj: Any, name: str, expected_len: int | None = None) -> list[Fraction]:
    arr = _req_list(obj, name)
    if expected_len is not None and len(arr) != expected_len:
        raise VerificationError(f"{name} length {len(arr)} != {expected_len}")
    return [_req_frac(x, f"{name}[{i}]") for i, x in enumerate(arr)]


def _req_int_mat(obj: Any, name: str, rows: int | None = None, cols: int | None = None) -> list[list[int]]:
    arr = _req_list(obj, name, max_len=MAX_MAT_DIM)
    if len(arr) == 0:
        raise VerificationError(f"{name} must be non-empty")
    if rows is not None and len(arr) != rows:
        raise VerificationError(f"{name} row count {len(arr)} != {rows}")
    c0: int | None = cols
    mat: list[list[int]] = []
    for i, r in enumerate(arr):
        rl = _req_list(r, f"{name}[{i}]", max_len=MAX_MAT_DIM)
        if len(rl) == 0:
            raise VerificationError(f"{name}[{i}] must be non-empty")
        if c0 is None:
            c0 = len(rl)
        elif len(rl) != c0:
            raise VerificationError(f"{name}[{i}] col count {len(rl)} != {c0}")
        mat.append([_req_int(v, f"{name}[{i}][{j}]", -MAX_RATIONAL_ABS, MAX_RATIONAL_ABS) for j, v in enumerate(rl)])
    return mat


def _req_frac_mat(obj: Any, name: str, rows: int | None = None, cols: int | None = None) -> list[list[Fraction]]:
    arr = _req_list(obj, name, max_len=MAX_MAT_DIM)
    if len(arr) == 0:
        raise VerificationError(f"{name} must be non-empty")
    if rows is not None and len(arr) != rows:
        raise VerificationError(f"{name} row count {len(arr)} != {rows}")
    c0: int | None = cols
    mat: list[list[Fraction]] = []
    for i, r in enumerate(arr):
        rl = _req_list(r, f"{name}[{i}]", max_len=MAX_MAT_DIM)
        if len(rl) == 0:
            raise VerificationError(f"{name}[{i}] must be non-empty")
        if c0 is None:
            c0 = len(rl)
        elif len(rl) != c0:
            raise VerificationError(f"{name}[{i}] col count {len(rl)} != {c0}")
        mat.append([_req_frac(v, f"{name}[{i}][{j}]") for j, v in enumerate(rl)])
    return mat


def _is_prime(p: int) -> bool:
    if p < 2:
        return False
    d = 2
    while d * d <= p:
        if p % d == 0:
            return False
        d += 1
    return True


def _mat_rank_q(mat: list[list[int]]) -> int:
    """Exact rank over Q using Fraction Gaussian elimination without float or /."""
    m = len(mat)
    n = len(mat[0])
    a = [[Fraction(v, 1) for v in row] for row in mat]
    rank = 0
    for col in range(n):
        piv = -1
        for r in range(rank, m):
            if a[r][col] != 0:
                piv = r
                break
        if piv == -1:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        pv = a[rank][col]
        inv_pv = Fraction(pv.denominator, pv.numerator)
        for r in range(rank + 1, m):
            if a[r][col] != 0:
                factor = a[r][col] * inv_pv
                for c in range(col, n):
                    a[r][c] = a[r][c] - factor * a[rank][c]
        rank += 1
        if rank == m:
            break
    return rank


def _mat_rank_fp(mat: list[list[int]], p: int) -> int:
    """Exact rank over F_p using modular arithmetic."""
    m = len(mat)
    n = len(mat[0])
    a = [[v % p for v in row] for row in mat]
    rank = 0
    for col in range(n):
        piv = -1
        for r in range(rank, m):
            if a[r][col] % p != 0:
                piv = r
                break
        if piv == -1:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        inv_pv = pow(a[rank][col], p - 2, p)
        for r in range(rank + 1, m):
            if a[r][col] != 0:
                factor = (a[r][col] * inv_pv) % p
                for c in range(col, n):
                    a[r][c] = (a[r][c] - factor * a[rank][c]) % p
        rank += 1
        if rank == m:
            break
    return rank


def _mat_det_int(mat: list[list[int]]) -> int:
    """Exact integer determinant via Fraction Gaussian elimination (no float, no /)."""
    n = len(mat)
    a = [[Fraction(v, 1) for v in row] for row in mat]
    sign = 1
    det = Fraction(1, 1)
    for col in range(n):
        piv = -1
        for r in range(col, n):
            if a[r][col] != 0:
                piv = r
                break
        if piv == -1:
            return 0
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            sign = -sign
        pv = a[col][col]
        det = det * pv
        inv_pv = Fraction(pv.denominator, pv.numerator)
        for r in range(col + 1, n):
            if a[r][col] != 0:
                factor = a[r][col] * inv_pv
                for c in range(col, n):
                    a[r][c] = a[r][c] - factor * a[col][c]
    det = det * sign
    if det.denominator != 1:
        raise VerificationError("Internal determinant denominator != 1 for integer matrix")
    return det.numerator


def _block_diag(u: list[list[int]], v: list[list[int]]) -> list[list[int]]:
    ru, cu = len(u), len(u[0])
    rv, cv = len(v), len(v[0])
    out: list[list[int]] = []
    for r in range(ru):
        out.append(list(u[r]) + [0] * cv)
    for r in range(rv):
        out.append([0] * cu + list(v[r]))
    return out


def _kron(u: list[list[int]], v: list[list[int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for ur in u:
        for vr in v:
            row: list[int] = []
            for x in ur:
                for y in vr:
                    row.append(x * y)
            out.append(row)
    return out


def _binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    num = 1
    den = 1
    kk = k if k <= n - k else n - k
    for i in range(1, kk + 1):
        num *= n - kk + i
        den *= i
    return num // den


def _subsets(q: int, b: int) -> list[tuple[int, ...]]:
    """Lexicographically ordered b-element subsets of {0, ..., q-1}."""
    if b == 0:
        return [()]
    if b < 0 or b > q:
        return []
    res: list[tuple[int, ...]] = []

    def rec(start: int, cur: list[int]) -> None:
        if len(cur) == b:
            res.append(tuple(cur))
            return
        for x in range(start, q - (b - len(cur)) + 1):
            cur.append(x)
            rec(x + 1, cur)
            cur.pop()

    rec(0, [])
    return res


def _exterior_ops(q: int, b: int) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], list[tuple[int, ...]], dict[tuple[int, int, int, int], int], dict[tuple[int, int, int, int], int]]:
    """Construct exact integer matrices (D_j L_i)_{beta, alpha} and (L_i D_j)_{beta, alpha}
    on Lambda^b Z^q.
    Returns bases (B_b, B_plus, B_minus) and dicts mapping (i, j, alpha, beta) -> int entry.
    """
    b_cur = _subsets(q, b)
    b_plus = _subsets(q, b + 1)
    b_minus = _subsets(q, b - 1)
    idx_plus = {s: idx for idx, s in enumerate(b_plus)}
    idx_minus = {s: idx for idx, s in enumerate(b_minus)}

    # L_i: B_b -> B_plus
    # L_i(e_S) = (-1)^(pos of i in S union {i}) e_{S union {i}} if i not in S else 0
    l_up: list[list[tuple[int, int]]] = [[(-1, 0) for _ in b_cur] for _ in range(q)]
    for i in range(q):
        for a_idx, s in enumerate(b_cur):
            if i in s:
                continue
            pos = 0
            while pos < len(s) and s[pos] < i:
                pos += 1
            ns = s[:pos] + (i,) + s[pos:]
            sign = 1 if (pos % 2 == 0) else -1
            l_up[i][a_idx] = (idx_plus[ns], sign)

    # D_j: B_plus -> B_cur
    # D_j(e_T) = (-1)^(pos of j in T) e_{T \ {j}} if j in T else 0
    idx_cur = {s: idx for idx, s in enumerate(b_cur)}
    d_down_from_plus: list[list[tuple[int, int]]] = [[(-1, 0) for _ in b_plus] for _ in range(q)]
    for j in range(q):
        for t_idx, t in enumerate(b_plus):
            if j not in t:
                continue
            pos = t.index(j)
            ns = t[:pos] + t[pos + 1 :]
            sign = 1 if (pos % 2 == 0) else -1
            d_down_from_plus[j][t_idx] = (idx_cur[ns], sign)

    # D_j: B_cur -> B_minus
    d_down_from_cur: list[list[tuple[int, int]]] = [[(-1, 0) for _ in b_cur] for _ in range(q)]
    for j in range(q):
        for a_idx, s in enumerate(b_cur):
            if j not in s:
                continue
            pos = s.index(j)
            ns = s[:pos] + s[pos + 1 :]
            sign = 1 if (pos % 2 == 0) else -1
            d_down_from_cur[j][a_idx] = (idx_minus[ns], sign)

    # L_i: B_minus -> B_cur
    l_up_from_minus: list[list[tuple[int, int]]] = [[(-1, 0) for _ in b_minus] for _ in range(q)]
    for i in range(q):
        for m_idx, ms in enumerate(b_minus):
            if i in ms:
                continue
            pos = 0
            while pos < len(ms) and ms[pos] < i:
                pos += 1
            ns = ms[:pos] + (i,) + ms[pos:]
            sign = 1 if (pos % 2 == 0) else -1
            l_up_from_minus[i][m_idx] = (idx_cur[ns], sign)

    dl_map: dict[tuple[int, int, int, int], int] = {}
    ld_map: dict[tuple[int, int, int, int], int] = {}
    for i in range(q):
        for j in range(q):
            for a_idx in range(len(b_cur)):
                t_idx, s1 = l_up[i][a_idx]
                if t_idx >= 0:
                    b_idx, s2 = d_down_from_plus[j][t_idx]
                    if b_idx >= 0:
                        dl_map[(i, j, a_idx, b_idx)] = dl_map.get((i, j, a_idx, b_idx), 0) + s1 * s2
                m_idx, s3 = d_down_from_cur[j][a_idx]
                if m_idx >= 0:
                    b_idx, s4 = l_up_from_minus[i][m_idx]
                    if b_idx >= 0:
                        ld_map[(i, j, a_idx, b_idx)] = ld_map.get((i, j, a_idx, b_idx), 0) + s3 * s4
    return b_cur, b_plus, b_minus, dl_map, ld_map


# ---------------------------------------------------------------------------
# Mode 1: flattening_cut_nonadditivity (Sec. 1.1 & Sec. 2.1)
# ---------------------------------------------------------------------------
def _verify_mode1(cert: dict[str, Any]) -> dict[str, Any]:
    field_p = _req_int(cert.get("field_char"), "field_char", 0, 10007)
    if field_p != 0 and not _is_prime(field_p):
        raise VerificationError(f"field_char={field_p} must be 0 (Q) or prime")

    cuts = _req_list(cert.get("cuts"), "cuts", max_len=MAX_CUTS)
    if len(cuts) < 2:
        raise VerificationError("cuts must contain at least 2 flattening cuts")

    rank_fn = _mat_rank_q if field_p == 0 else (lambda m: _mat_rank_fp(m, field_p))

    prod_u = 1
    prod_v = 1
    prod_sum = 1
    cut_ranks: list[dict[str, int]] = []
    for idx, c_raw in enumerate(cuts):
        c = _req_dict(c_raw, f"cuts[{idx}]")
        _req_str(c.get("cut_label"), f"cuts[{idx}].cut_label")
        mat_u = _req_int_mat(c.get("matrix_U"), f"cuts[{idx}].matrix_U")
        mat_v = _req_int_mat(c.get("matrix_V"), f"cuts[{idx}].matrix_V")
        exp_ru = _req_int(c.get("expected_rank_U"), f"cuts[{idx}].expected_rank_U", 1, MAX_MAT_DIM)
        exp_rv = _req_int(c.get("expected_rank_V"), f"cuts[{idx}].expected_rank_V", 1, MAX_MAT_DIM)

        ru = rank_fn(mat_u)
        rv = rank_fn(mat_v)
        if ru != exp_ru:
            raise VerificationError(f"cuts[{idx}] rank(U)={ru} != expected {exp_ru}")
        if rv != exp_rv:
            raise VerificationError(f"cuts[{idx}] rank(V)={rv} != expected {exp_rv}")

        # Verify direct-sum additivity and tensor-product multiplicativity on this cut
        mat_dsum = _block_diag(mat_u, mat_v)
        r_dsum = rank_fn(mat_dsum)
        if r_dsum != ru + rv:
            raise VerificationError(f"cuts[{idx}] direct-sum rank {r_dsum} != {ru + rv}")

        if len(mat_u) * len(mat_v) <= MAX_MAT_DIM and len(mat_u[0]) * len(mat_v[0]) <= MAX_MAT_DIM:
            mat_tprod = _kron(mat_u, mat_v)
            r_tprod = rank_fn(mat_tprod)
            if r_tprod != ru * rv:
                raise VerificationError(f"cuts[{idx}] tensor-product rank {r_tprod} != {ru * rv}")

        prod_u *= ru
        prod_v *= rv
        prod_sum *= r_dsum
        cut_ranks.append({"rank_U": ru, "rank_V": rv, "rank_U_oplus_V": r_dsum})

    # Verify geometric-mean non-additivity for 2-cut product or general k-cut Holder/Minkowski strict gap:
    # For k cuts: G(U+V)^k = prod(r_U + r_V), and if k = 2:
    #   G(U+V) > G(U) + G(V) <=> prod_sum > prod_u + prod_v and (prod_sum - prod_u - prod_v)^2 > 4*prod_u*prod_v.
    # For general k >= 2, check exact k-th root lower witness or exact k=2 identity.
    expect_nonadd = _req_bool(cert.get("certify_nonadditivity"), "certify_nonadditivity")
    nonadd_certified = False
    if len(cuts) == 2:
        diff = prod_sum - prod_u - prod_v
        is_strict = (diff > 0) and (diff * diff > 4 * prod_u * prod_v)
        if expect_nonadd and not is_strict:
            raise VerificationError(
                f"2-cut geometric mean non-additivity failed: prod_sum={prod_sum}, prod_u={prod_u}, prod_v={prod_v}"
            )
        if not expect_nonadd and is_strict:
            raise VerificationError("Expected additive 2-cut example but got strict non-additivity")
        nonadd_certified = is_strict
    else:
        # For k >= 3 cuts, check via exact integer rational brackets u_hi^k >= prod_u, v_hi^k >= prod_v,
        # (u_hi + v_hi)^k < prod_sum.
        if expect_nonadd:
            u_hi = _req_frac(cert.get("geom_upper_U"), "geom_upper_U")
            v_hi = _req_frac(cert.get("geom_upper_V"), "geom_upper_V")
            if u_hi <= 0 or v_hi <= 0:
                raise VerificationError("geom_upper_U and geom_upper_V must be positive")
            k = len(cuts)
            if u_hi**k < prod_u or v_hi**k < prod_v:
                raise VerificationError("Upper brackets do not dominate G(U)^k or G(V)^k")
            if not ((u_hi + v_hi) ** k < prod_sum):
                raise VerificationError("(u_hi + v_hi)^k is not strictly less than prod_sum")
            nonadd_certified = True
        else:
            # All ratios proportional: ru_i * rv_0 == ru_0 * rv_i
            r0u = cut_ranks[0]["rank_U"]
            r0v = cut_ranks[0]["rank_V"]
            for cr in cut_ranks[1:]:
                if cr["rank_U"] * r0v != r0u * cr["rank_V"]:
                    raise VerificationError("Expected proportional additive cuts, found non-proportional ranks")

    # Optional Hamming-square / 4-point CAT(0) violation check (Sec. 1.1 & Eq. 4.5)
    hamming_violation = False
    if "four_point_test" in cert and cert["four_point_test"] is not None:
        fp = _req_dict(cert["four_point_test"], "four_point_test")
        d_ab = _req_frac(fp.get("d_AB"), "four_point_test.d_AB")
        d_bc = _req_frac(fp.get("d_BC"), "four_point_test.d_BC")
        d_cd = _req_frac(fp.get("d_CD"), "four_point_test.d_CD")
        d_da = _req_frac(fp.get("d_DA"), "four_point_test.d_DA")
        d_ac = _req_frac(fp.get("d_AC"), "four_point_test.d_AC")
        d_bd = _req_frac(fp.get("d_BD"), "four_point_test.d_BD")
        for nm, val in [("d_AB", d_ab), ("d_BC", d_bc), ("d_CD", d_cd), ("d_DA", d_da), ("d_AC", d_ac), ("d_BD", d_bd)]:
            if val < 0:
                raise VerificationError(f"four_point_test.{nm} must be nonnegative")
        # Verify all 12 triangle inequalities across the 4-point symmetric distance matrix
        d4 = [
            [Fraction(0, 1), d_ab, d_ac, d_da],
            [d_ab, Fraction(0, 1), d_bc, d_bd],
            [d_ac, d_bc, Fraction(0, 1), d_cd],
            [d_da, d_bd, d_cd, Fraction(0, 1)],
        ]
        for i4 in range(4):
            for j4 in range(4):
                for k4 in range(4):
                    if d4[i4][k4] > d4[i4][j4] + d4[j4][k4]:
                        raise VerificationError(
                            f"four_point_test violates triangle inequality on ({i4},{j4},{k4})"
                        )
        lhs = d_ac * d_ac + d_bd * d_bd
        rhs = d_ab * d_ab + d_bc * d_bc + d_cd * d_cd + d_da * d_da
        exp_viol = _req_bool(fp.get("expect_cat0_violation"), "four_point_test.expect_cat0_violation")
        if exp_viol:
            if not (lhs > rhs):
                raise VerificationError(f"Expected Eq. (4.5) violation lhs={lhs} > rhs={rhs}")
            hamming_violation = True
        else:
            if lhs > rhs:
                raise VerificationError(f"Unexpected Eq. (4.5) violation lhs={lhs} > rhs={rhs}")

    return {
        "mode": "flattening_cut_nonadditivity",
        "field_char": field_p,
        "num_cuts": len(cuts),
        "prod_U": prod_u,
        "prod_V": prod_v,
        "prod_U_oplus_V": prod_sum,
        "nonadditivity_certified": nonadd_certified,
        "cat0_violation_certified": hamming_violation,
    }


# ---------------------------------------------------------------------------
# Mode 2: exterior_crossing_fermion (Eq. 3.2, Eq. 3.3, Theorem 4.1)
# ---------------------------------------------------------------------------
def _verify_mode2(cert: dict[str, Any]) -> dict[str, Any]:
    submode = _req_str(cert.get("submode"), "submode")
    if submode not in ("exterior_crossing", "vandermonde_unimodular", "explicit_crossing_factors"):
        raise VerificationError(f"Unknown Mode 2 submode: {submode!r}")

    if submode == "exterior_crossing":
        q = _req_int(cert.get("Q"), "Q", 2, MAX_EXT_Q)
        b = _req_int(cert.get("b"), "b", 1, q - 1)
        exp_r = _req_int(cert.get("expected_R"), "expected_R", 1)
        exp_rp = _req_int(cert.get("expected_R_plus"), "expected_R_plus", 1)
        exp_rm = _req_int(cert.get("expected_R_minus"), "expected_R_minus", 1)

        b_cur, b_plus, b_minus, dl_map, ld_map = _exterior_ops(q, b)
        r_val = len(b_cur)
        rp_val = len(b_plus)
        rm_val = len(b_minus)
        if r_val != exp_r or rp_val != exp_rp or rm_val != exp_rm:
            raise VerificationError(
                f"Exterior dimensions ({r_val},{rp_val},{rm_val}) != expected ({exp_r},{exp_rp},{exp_rm})"
            )
        # Verify Eq. (3.3) exact binomial identities
        if (b + 1) * rp_val != (q - b) * r_val:
            raise VerificationError("Binomial ratio R+/R identity failed")
        if (q - b + 1) * rm_val != b * r_val:
            raise VerificationError("Binomial ratio R-/R identity failed")

        # Verify Eq. (3.2) D_j L_i + L_i D_j = delta_{ij} I on ALL (i, j, alpha, beta) over Z!
        nonzeros_dl = 0
        nonzeros_ld = 0
        for i in range(q):
            for j in range(q):
                for a_idx in range(r_val):
                    for b_idx in range(r_val):
                        v_dl = dl_map.get((i, j, a_idx, b_idx), 0)
                        v_ld = ld_map.get((i, j, a_idx, b_idx), 0)
                        if v_dl != 0:
                            nonzeros_dl += 1
                        if v_ld != 0:
                            nonzeros_ld += 1
                        target = 1 if (i == j and a_idx == b_idx) else 0
                        if v_dl + v_ld != target:
                            raise VerificationError(
                                f"Eq. (3.2) failed at (i={i},j={j},a={a_idx},b={b_idx}): {v_dl}+{v_ld}!={target}"
                            )
        if nonzeros_dl == 0 or nonzeros_ld == 0:
            raise VerificationError("Vacuous exterior decomposition: zero summand")

        # Always compute exact flattening ranks of T_+ on cut (i, a) | (j, b) and T_- on cut (i, b) | (j, a)
        # For MAX_EXT_Q = 6, max Q * R = 6 * 20 = 120, fast in exact Fraction elimination.
        dim_flat = q * r_val
        t_plus = [[0] * dim_flat for _ in range(dim_flat)]
        t_minus = [[0] * dim_flat for _ in range(dim_flat)]
        for i in range(q):
            for j in range(q):
                for a_idx in range(r_val):
                    for b_idx in range(r_val):
                        t_plus[i * r_val + a_idx][j * r_val + b_idx] = dl_map.get((i, j, a_idx, b_idx), 0)
                        t_minus[i * r_val + b_idx][j * r_val + a_idx] = ld_map.get((i, j, a_idx, b_idx), 0)
        rk_p = _mat_rank_q(t_plus)
        rk_m = _mat_rank_q(t_minus)
        if rk_p != rp_val or rk_m != rm_val:
            raise VerificationError(
                f"Crossing flattening ranks ({rk_p},{rk_m}) != ({rp_val},{rm_val})"
            )

        return {
            "mode": "exterior_crossing_fermion",
            "submode": submode,
            "Q": q,
            "b": b,
            "R": r_val,
            "R_plus": rp_val,
            "R_minus": rm_val,
            "verified_crossing_ranks": [rk_p, rk_m],
        }

    if submode == "vandermonde_unimodular":
        # Proof of Theorem 4.1 (p. 19): P_{i,l} = C(i, l), (P P^T)_{i,j} = C(i+j, i),
        # and m-rowed minor matrix P^[m] restricted to energy downset h(I) = sum(i + 1/2) <= H
        # is unit lower triangular in lexicographic order, hence det(C^[m]) = +-1 over Z.
        two_h = _req_int(cert.get("two_H"), "two_H", 1, 20)
        m = _req_int(cert.get("block_size_m"), "block_size_m", 1, 4)
        exp_num_states = _req_int(cert.get("expected_num_states"), "expected_num_states", 1, MAX_MAT_DIM)

        # Enumerate m-element subsets I of {0, 1, ...} with sum_{i in I} (2i + 1) <= two_h
        max_mode = (two_h - 1) // 2
        states: list[tuple[int, ...]] = []
        for s in _subsets(max_mode + 1, m):
            energy2 = sum(2 * x + 1 for x in s)
            if energy2 <= two_h:
                states.append(s)
        if len(states) != exp_num_states:
            raise VerificationError(f"Energy downset size {len(states)} != expected {exp_num_states}")
        if len(states) > MAX_MAT_DIM:
            raise VerificationError("Energy downset exceeds MAX_MAT_DIM")

        # Verify P P^T = [C(i+j, i)] on 0..max_mode
        n_modes = max_mode + 1
        p_mat = [[_binom(i, l) for l in range(n_modes)] for i in range(n_modes)]
        for i in range(n_modes):
            for j in range(n_modes):
                dot = sum(p_mat[i][l] * p_mat[j][l] for l in range(n_modes))
                if dot != _binom(i + j, i):
                    raise VerificationError(f"Vandermonde failed at ({i},{j})")

        # Build compound matrix P^[m] and C^[m] on states (lexicographic order)
        k_st = len(states)
        p_comp = [[0] * k_st for _ in range(k_st)]
        c_block = [[0] * k_st for _ in range(k_st)]
        for r_idx, row_s in enumerate(states):
            for c_idx, col_s in enumerate(states):
                sub_p = [[p_mat[i][l] for l in col_s] for i in row_s]
                p_comp[r_idx][c_idx] = _mat_det_int(sub_p)
                sub_k = [[_binom(i + j, i) for j in col_s] for i in row_s]
                sign_exp = (m * (m - 1)) // 2 + sum(row_s)
                sgn = 1 if (sign_exp % 2 == 0) else -1
                c_block[r_idx][c_idx] = sgn * _mat_det_int(sub_k)

        # Check P^[m] is unit lower triangular in lexicographic order
        for r_idx in range(k_st):
            if p_comp[r_idx][r_idx] != 1:
                raise VerificationError(f"P^[m] diagonal at {r_idx} is {p_comp[r_idx][r_idx]} != 1")
            for c_idx in range(r_idx + 1, k_st):
                if p_comp[r_idx][c_idx] != 0:
                    raise VerificationError(f"P^[m] not lower triangular at ({r_idx},{c_idx})")

        # Verify Cauchy-Binet identity C_IJ = sgn(I) * (P^[m] (P^[m])^T)_{I,J}
        for r_idx, row_s in enumerate(states):
            sign_exp = (m * (m - 1)) // 2 + sum(row_s)
            sgn = 1 if (sign_exp % 2 == 0) else -1
            for c_idx in range(k_st):
                cb = sum(p_comp[r_idx][t] * p_comp[c_idx][t] for t in range(k_st))
                if c_block[r_idx][c_idx] != sgn * cb:
                    raise VerificationError(f"Cauchy-Binet mismatch at ({r_idx},{c_idx})")

        det_c = _mat_det_int(c_block)
        if det_c not in (1, -1):
            raise VerificationError(f"det(C^[m]) = {det_c}, not unimodular")

        return {
            "mode": "exterior_crossing_fermion",
            "submode": submode,
            "two_H": two_h,
            "block_size_m": m,
            "num_states": k_st,
            "det_C_block": det_c,
        }

    # submode == "explicit_crossing_factors"
    dx = _req_int(cert.get("d_X"), "d_X", 2, 8)
    dy = _req_int(cert.get("d_Y"), "d_Y", 2, 8)
    if dx * dy > MAX_MAT_DIM:
        raise VerificationError("d_X * d_Y exceeds MAX_MAT_DIM")
    dim = dx * dy
    rp = _req_int(cert.get("r_plus"), "r_plus", 1, dim - 1)
    rm = _req_int(cert.get("r_minus"), "r_minus", 1, dim - 1)
    if rp + rm >= dim:
        raise VerificationError(
            f"Trivial crossing split: r_plus + r_minus = {rp + rm} >= d_X*d_Y = {dim}"
        )

    a_plus = _req_int_mat(cert.get("A_plus"), "A_plus", rows=dim, cols=rp)
    b_plus = _req_int_mat(cert.get("B_plus"), "B_plus", rows=rp, cols=dim)
    a_minus = _req_int_mat(cert.get("A_minus"), "A_minus", rows=dim, cols=rm)
    b_minus = _req_int_mat(cert.get("B_minus"), "B_minus", rows=rm, cols=dim)

    # T_+ indexed by row = (x+, y+), col = (x-, y-)
    # T_- indexed by row = (x+, y-), col = (x-, y+)
    nz_p = 0
    nz_m = 0
    for xp in range(dx):
        for xm in range(dx):
            for yp in range(dy):
                for ym in range(dy):
                    row_p = xp * dy + yp
                    col_p = xm * dy + ym
                    val_p = sum(a_plus[row_p][k] * b_plus[k][col_p] for k in range(rp))
                    row_m = xp * dy + ym
                    col_m = xm * dy + yp
                    val_m = sum(a_minus[row_m][k] * b_minus[k][col_m] for k in range(rm))
                    if val_p != 0:
                        nz_p += 1
                    if val_m != 0:
                        nz_m += 1
                    target = 1 if (xp == xm and yp == ym) else 0
                    if val_p + val_m != target:
                        raise VerificationError(
                            f"Crossing sum mismatch at ({xp},{xm},{yp},{ym}): {val_p}+{val_m}!={target}"
                        )
    if nz_p == 0 or nz_m == 0:
        raise VerificationError("Vacuous crossing factorization: one summand is identically zero")

    return {
        "mode": "exterior_crossing_fermion",
        "submode": submode,
        "d_X": dx,
        "d_Y": dy,
        "r_plus": rp,
        "r_minus": rm,
    }


# ---------------------------------------------------------------------------
# Mode 3: cat0_flat_mean_width (Secs. 2.2, 4.2, 4.3, 4.4, 5.1, 5.2)
# ---------------------------------------------------------------------------
def _verify_mode3(cert: dict[str, Any]) -> dict[str, Any]:
    submode = _req_str(cert.get("submode"), "submode")
    if submode not in ("cat0_quadrilateral_midpoint", "norm_rigidity_violation", "polygon_mean_width_copying", "copying_quadratic_degeneration"):
        raise VerificationError(f"Unknown Mode 3 submode: {submode!r}")

    if submode == "cat0_quadrilateral_midpoint":
        n = _req_int(cert.get("num_sites"), "num_sites", 4, MAX_SITES)
        d_mat = _req_frac_mat(cert.get("distance_matrix"), "distance_matrix", rows=n, cols=n)
        # Verify pseudometric axioms
        for i in range(n):
            if d_mat[i][i] != 0:
                raise VerificationError(f"d[{i}][{i}] != 0")
            for j in range(n):
                if d_mat[i][j] < 0 or d_mat[i][j] != d_mat[j][i]:
                    raise VerificationError(f"Invalid distance at ({i},{j})")
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if d_mat[i][k] > d_mat[i][j] + d_mat[j][k]:
                        raise VerificationError(f"Triangle inequality failed at ({i},{j},{k})")

        # Check Eq. (4.5) on all ordered 4-tuples of distinct sites
        quads_checked = 0
        exact_rect_count = 0
        for a in range(n):
            for b in range(n):
                if b == a:
                    continue
                for c in range(n):
                    if c in (a, b):
                        continue
                    for d in range(n):
                        if d in (a, b, c):
                            continue
                        lhs = d_mat[a][c] ** 2 + d_mat[b][d] ** 2
                        rhs = d_mat[a][b] ** 2 + d_mat[b][c] ** 2 + d_mat[c][d] ** 2 + d_mat[d][a] ** 2
                        if lhs > rhs:
                            raise VerificationError(f"CAT(0) Eq. (4.5) failed on ({a},{b},{c},{d}): {lhs} > {rhs}")
                        if lhs == rhs:
                            exact_rect_count += 1
                        quads_checked += 1

        # Check Lemma 4.4 / Eq. (4.4) on specified geodesic interpolation triples
        interp_checks = _req_list(cert.get("geodesic_checks"), "geodesic_checks")
        if len(interp_checks) == 0:
            raise VerificationError("geodesic_checks must be non-empty")
        for idx, g_raw in enumerate(interp_checks):
            g = _req_dict(g_raw, f"geodesic_checks[{idx}]")
            x = _req_int(g.get("x"), f"geodesic_checks[{idx}].x", 0, n - 1)
            y = _req_int(g.get("y"), f"geodesic_checks[{idx}].y", 0, n - 1)
            m = _req_int(g.get("m"), f"geodesic_checks[{idx}].m", 0, n - 1)
            z = _req_int(g.get("z"), f"geodesic_checks[{idx}].z", 0, n - 1)
            s = _req_frac(g.get("s"), f"geodesic_checks[{idx}].s")
            if s <= 0 or s >= 1:
                raise VerificationError(f"geodesic_checks[{idx}].s must lie in (0, 1)")
            if x == y or d_mat[x][y] <= 0:
                raise VerificationError(f"geodesic_checks[{idx}] requires distinct endpoints with d(x,y) > 0")
            if z in (x, y, m):
                raise VerificationError(f"geodesic_checks[{idx}] degenerate observer z in {{x,y,m}}")
            # Verify m is the s-interval point [x, y]_s
            if d_mat[x][m] != s * d_mat[x][y] or d_mat[m][y] != (1 - s) * d_mat[x][y]:
                raise VerificationError(f"geodesic_checks[{idx}] point m={m} is not [x,y]_s")
            # Eq. (4.4): d(z, m)^2 <= (1-s) d(z,x)^2 + s d(z,y)^2 - s(1-s) d(x,y)^2
            lhs = d_mat[z][m] ** 2
            rhs = (1 - s) * (d_mat[z][x] ** 2) + s * (d_mat[z][y] ** 2) - s * (1 - s) * (d_mat[x][y] ** 2)
            if lhs > rhs:
                raise VerificationError(f"CAT(0) Eq. (4.4) failed at geodesic_checks[{idx}]: {lhs} > {rhs}")
            expect_eq = _req_bool(g.get("expect_equality"), f"geodesic_checks[{idx}].expect_equality")
            if expect_eq and lhs != rhs:
                raise VerificationError(f"Expected CAT(0) equality at geodesic_checks[{idx}], got {lhs} < {rhs}")
            if not expect_eq and lhs == rhs:
                raise VerificationError(f"Expected strict negative curvature at geodesic_checks[{idx}], got equality")

        return {
            "mode": "cat0_flat_mean_width",
            "submode": submode,
            "num_sites": n,
            "quads_checked": quads_checked,
            "exact_rect_quads": exact_rect_count,
            "geodesic_checks": len(interp_checks),
        }

    if submode == "norm_rigidity_violation":
        # Lemma 4.2 Eq. (4.2): ||t f + u g||_*^2 + ||t f - u g||_*^2 > 2(t^2 + u^2)
        # Dual norm on R^2: l_inf (dual of l_1) or l_1 (dual of l_inf)
        dual_norm = _req_str(cert.get("dual_norm"), "dual_norm")
        if dual_norm not in ("l_inf", "l_1"):
            raise VerificationError("dual_norm must be 'l_inf' or 'l_1'")
        t = _req_frac(cert.get("t"), "t")
        u = _req_frac(cert.get("u"), "u")
        if t <= 0 or u <= 0:
            raise VerificationError("t and u must be strictly positive")
        f = _req_frac_vec(cert.get("f"), "f", expected_len=2)
        g = _req_frac_vec(cert.get("g"), "g", expected_len=2)
        delta = _req_frac(cert.get("delta"), "delta")

        def dnorm(v: list[Fraction]) -> Fraction:
            a0 = -v[0] if v[0] < 0 else v[0]
            a1 = -v[1] if v[1] < 0 else v[1]
            if dual_norm == "l_inf":
                return a0 if a0 >= a1 else a1
            return a0 + a1

        if dnorm(f) != 1 or dnorm(g) != 1:
            raise VerificationError("f and g must lie on the dual unit sphere")
        vp = [t * f[0] + u * g[0], t * f[1] + u * g[1]]
        vm = [t * f[0] - u * g[0], t * f[1] - u * g[1]]
        a_plus = dnorm(vp)
        a_minus = dnorm(vm)
        s_val = t * t + u * u
        if not (a_plus * a_plus + a_minus * a_minus > 2 * s_val):
            raise VerificationError("Dual parallelogram inequality Eq. (4.2) was not violated")
        if not (-s_val < delta < s_val):
            raise VerificationError("delta must satisfy -S < delta < S")
        if not (s_val - a_minus * a_minus < delta < a_plus * a_plus - s_val):
            raise VerificationError("delta does not lie in strict saving interval (S - A_-^2, A_+^2 - S)")
        return {
            "mode": "cat0_flat_mean_width",
            "submode": submode,
            "dual_norm": dual_norm,
            "S": [s_val.numerator, s_val.denominator],
            "A_plus_sq": [(a_plus * a_plus).numerator, (a_plus * a_plus).denominator],
            "A_minus_sq": [(a_minus * a_minus).numerator, (a_minus * a_minus).denominator],
            "delta": [delta.numerator, delta.denominator],
        }

    if submode == "polygon_mean_width_copying":
        # Theorem 5.3 planar polygon W(K) = per(K)/2 with exact rational edge lengths,
        # Proposition 5.2 copying increment, submodularity, and convex-hull zero-cost insertion.
        verts_raw = _req_list(cert.get("polygon_vertices"), "polygon_vertices", max_len=MAX_SITES)
        if len(verts_raw) < 3:
            raise VerificationError("polygon_vertices must have >= 3 vertices")
        verts = [_req_frac_vec(v, f"polygon_vertices[{i}]", expected_len=2) for i, v in enumerate(verts_raw)]
        k = len(verts)
        edge_lens = _req_frac_vec(cert.get("edge_lengths"), "edge_lengths", expected_len=k)

        # Verify global strict CCW convexity (every non-edge vertex strictly left of every directed edge)
        # This rules out double-wound polygons, star polygons, repeated vertices, and collinear triples.
        per = Fraction(0, 1)
        for i in range(k):
            p0 = verts[i]
            p1 = verts[(i + 1) % k]
            dx1, dy1 = p1[0] - p0[0], p1[1] - p0[1]
            for j in range(k):
                if j == i or j == (i + 1) % k:
                    continue
                pj = verts[j]
                cross = dx1 * (pj[1] - p0[1]) - dy1 * (pj[0] - p0[0])
                if cross <= 0:
                    raise VerificationError(f"Polygon not globally strictly CCW convex at edge {i} vs vertex {j}")
            if edge_lens[i] <= 0 or edge_lens[i] * edge_lens[i] != dx1 * dx1 + dy1 * dy1:
                raise VerificationError(f"Exact edge length mismatch at edge {i}")
            per = per + edge_lens[i]

        exp_w = _req_frac(cert.get("expected_mean_width"), "expected_mean_width")
        if exp_w != per * Fraction(1, 2):
            raise VerificationError(f"Mean width {exp_w} != per(K)/2 = {per * Fraction(1, 2)}")

        # Verify Eq. (2.1): diam(K) <= W(K) <= spanning tree (per - max_edge)
        max_diag_sq = Fraction(0, 1)
        for i in range(k):
            for j in range(i + 1, k):
                dsq = (verts[i][0] - verts[j][0]) ** 2 + (verts[i][1] - verts[j][1]) ** 2
                if dsq > max_diag_sq:
                    max_diag_sq = dsq
        if exp_w * exp_w < max_diag_sq:
            raise VerificationError("Mean width violates diam(K) <= W(K)")
        max_edge = max(edge_lens)
        tree_bound = per - max_edge
        if exp_w > tree_bound:
            raise VerificationError("Mean width exceeds spanning-tree bound")

        # Verify convex-hull zero-cost insertion z in conv(verts)
        ins_checked = 0
        if "interior_point" in cert and cert["interior_point"] is not None:
            ins = _req_dict(cert["interior_point"], "interior_point")
            z = _req_frac_vec(ins.get("z"), "interior_point.z", expected_len=2)
            lambdas = _req_frac_vec(ins.get("convex_weights"), "interior_point.convex_weights", expected_len=k)
            if sum(lambdas, Fraction(0, 1)) != 1:
                raise VerificationError("convex_weights must sum to 1")
            for idx_l, lam in enumerate(lambdas):
                if lam < 0:
                    raise VerificationError(f"convex_weights[{idx_l}] < 0")
            zx = sum((lambdas[i] * verts[i][0] for i in range(k)), Fraction(0, 1))
            zy = sum((lambdas[i] * verts[i][1] for i in range(k)), Fraction(0, 1))
            if zx != z[0] or zy != z[1]:
                raise VerificationError("Convex combination does not equal z")
            delta_with_z = _req_frac(ins.get("delta_with_z"), "interior_point.delta_with_z")
            if delta_with_z != exp_w:
                raise VerificationError("Adding z in conv(K) must preserve delta(K) = W(K)")
            ins_checked = 1

        # Optional submodularity check on two overlapping subsets A, B of polygon vertices (Prop. 5.2)
        submod_checked = 0
        if "submodularity_check" in cert and cert["submodularity_check"] is not None:
            sm = _req_dict(cert["submodularity_check"], "submodularity_check")
            d_a = _req_frac(sm.get("delta_A"), "submodularity_check.delta_A")
            d_b = _req_frac(sm.get("delta_B"), "submodularity_check.delta_B")
            d_inter = _req_frac(sm.get("delta_A_inter_B"), "submodularity_check.delta_A_inter_B")
            d_union = _req_frac(sm.get("delta_A_union_B"), "submodularity_check.delta_A_union_B")
            idx_a = [_req_int(x, "idx_A", 0, k - 1) for x in _req_list(sm.get("indices_A"), "indices_A")]
            idx_b = [_req_int(x, "idx_B", 0, k - 1) for x in _req_list(sm.get("indices_B"), "indices_B")]
            set_a, set_b = set(idx_a), set(idx_b)
            if len(set_a) != len(idx_a) or len(set_b) != len(idx_b):
                raise VerificationError("Duplicate vertex indices in indices_A or indices_B")
            if len(set_a & set_b) == 0:
                raise VerificationError("Prop. 5.2 submodularity requires A intersect B nonempty")
            if len(set_a | set_b) != k:
                raise VerificationError("submodularity_check A union B must cover the polygon vertices")
            if d_union != exp_w:
                raise VerificationError("delta_A_union_B != W(K)")
            # Exact perimeter/2 recomputation for subsets of vertices on a convex polygon
            # (subset vertices in cyclic order form a convex sub-polygon!)
            def _sub_mean_width(idx_set: set[int]) -> Fraction:
                ordered = sorted(idx_set)
                m_sub = len(ordered)
                if m_sub <= 1:
                    return Fraction(0, 1)
                if m_sub == 2:
                    p0, p1 = verts[ordered[0]], verts[ordered[1]]
                    dsq = (p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2
                    # Must match an exact rational distance; check via integer sqrt on cleared fraction
                    return _exact_frac_sqrt(dsq)
                tot = Fraction(0, 1)
                for ii in range(m_sub):
                    p0 = verts[ordered[ii]]
                    p1 = verts[ordered[(ii + 1) % m_sub]]
                    dsq = (p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2
                    tot = tot + _exact_frac_sqrt(dsq)
                return tot * Fraction(1, 2)

            if _sub_mean_width(set_a) != d_a:
                raise VerificationError("delta_A mismatch against exact sub-polygon mean width")
            if _sub_mean_width(set_b) != d_b:
                raise VerificationError("delta_B mismatch against exact sub-polygon mean width")
            if _sub_mean_width(set_a & set_b) != d_inter:
                raise VerificationError("delta_A_inter_B mismatch against exact intersection mean width")
            if d_union + d_inter > d_a + d_b:
                raise VerificationError("Prop. 5.2 submodularity violated")
            submod_checked = 1

        return {
            "mode": "cat0_flat_mean_width",
            "submode": submode,
            "num_vertices": k,
            "mean_width": [exp_w.numerator, exp_w.denominator],
            "interior_insertion_checked": ins_checked,
            "submodularity_checked": submod_checked,
        }

    # submode == "copying_quadratic_degeneration" (Prop. 5.2, p. 24)
    r = _req_int(cert.get("r"), "r", 2, MAX_SITES)
    q = _req_int(cert.get("q"), "q", 2, 1000)
    h_val = _req_int(cert.get("h"), "h", 0, q - 1)
    v_list = [_req_int(x, f"v[{i}]", 0, q - 1) for i, x in enumerate(_req_list(cert.get("v"), "v"))]
    if len(v_list) != r:
        raise VerificationError(f"v length {len(v_list)} != r={r}")
    g_val = _req_int(cert.get("g"), "g", 0, q - 1)
    # Check Prop. 5.2 linear equation g = h + sum_{j=2}^r v_j - (r-1)*v_1
    rhs_lin = h_val + sum(v_list[1:]) - (r - 1) * v_list[0]
    if g_val != rhs_lin:
        raise VerificationError(f"Prop. 5.2 linear constraint failed: g={g_val} != {rhs_lin}")
    # Check monomial degree cap L = 2*r*(q-1) + 1
    exp_l = _req_int(cert.get("expected_L"), "expected_L", 1)
    if exp_l != 2 * r * (q - 1) + 1:
        raise VerificationError(f"expected_L={exp_l} != 2r(q-1)+1={2 * r * (q - 1) + 1}")
    # First degeneration weight: (g - h)^2 = g^2 - 2*h*(h + sum v_j - (r-1)v_1) + h^2
    w1 = (g_val - h_val) ** 2
    w1_local = g_val * g_val - 2 * h_val * (h_val + sum(v_list[1:]) - (r - 1) * v_list[0]) + h_val * h_val
    if w1 != w1_local or w1 < 0:
        raise VerificationError("First quadratic degeneration weight identity failed")
    exp_surv1 = _req_bool(cert.get("survives_first_degeneration"), "survives_first_degeneration")
    if exp_surv1 != (w1 == 0):
        raise VerificationError("survives_first_degeneration flag mismatch")

    surv2 = False
    w2 = 0
    if w1 != 0:
        if "survives_second_degeneration" in cert and _req_bool(
            cert.get("survives_second_degeneration"), "survives_second_degeneration"
        ):
            raise VerificationError("survives_second_degeneration cannot be True when w1 != 0")
    else:
        # Remaining equation: sum_{j=2}^r v_j = (r-1) v_1
        if sum(v_list[1:]) != (r - 1) * v_list[0]:
            raise VerificationError("Hyperplane invariant violated after g=h")
        w2_local = sum(v * v for v in v_list[1:]) - (r - 1) * v_list[0] * v_list[0]
        w2 = sum((v - v_list[0]) ** 2 for v in v_list[1:])
        if w2_local != w2 or w2 < 0:
            raise VerificationError("Second quadratic degeneration identity failed")
        exp_surv2 = _req_bool(cert.get("survives_second_degeneration"), "survives_second_degeneration")
        surv2 = w2 == 0
        if exp_surv2 != surv2:
            raise VerificationError("survives_second_degeneration flag mismatch")
        if surv2 and any(v != v_list[0] for v in v_list[1:]):
            raise VerificationError("Zero second weight must imply all v_j equal")

    return {
        "mode": "cat0_flat_mean_width",
        "submode": submode,
        "r": r,
        "q": q,
        "L": exp_l,
        "weight1": w1,
        "weight2": w2,
        "survives_both": (w1 == 0 and surv2),
    }


def _exact_isqrt(n: int) -> int:
    if n < 0:
        raise VerificationError("Negative value inside exact square root")
    if n == 0:
        return 0
    lo, hi = 1, n
    while lo <= hi:
        mid = (lo + hi) // 2
        sq = mid * mid
        if sq == n:
            return mid
        if sq < n:
            lo = mid + 1
        else:
            hi = mid - 1
    raise VerificationError(f"Non-square integer {n} in exact rational Euclidean length")


def _exact_frac_sqrt(f: Fraction) -> Fraction:
    return Fraction(_exact_isqrt(f.numerator), _exact_isqrt(f.denominator))


# ---------------------------------------------------------------------------
# Mode 4: catalyst_gap_routing_budget (Secs. 8, 9, 16)
# ---------------------------------------------------------------------------
def _verify_mode4(cert: dict[str, Any]) -> dict[str, Any]:
    submode = _req_str(cert.get("submode"), "submode")
    if submode not in ("tag_perm_and_catalyst_gap", "routing_budget_and_bregman", "diagonal_contradiction_gap"):
        raise VerificationError(f"Unknown Mode 4 submode: {submode!r}")

    if submode == "tag_perm_and_catalyst_gap":
        # Definition 8.1 Tag(h, w) + Section 8.1 / Lemma 8.2 catalyst gap & substitution
        h = _req_int(cert.get("h"), "h", 1, 4)
        w = _req_int(cert.get("w"), "w", 1, 16)
        n = 1 << h
        payloads = _req_list(cert.get("payloads"), "payloads", max_len=n)
        records = _req_list(cert.get("records"), "records", max_len=n)
        if len(payloads) != n or len(records) != n:
            raise VerificationError(f"Tag(h={h},w={w}) requires n={n} payloads and records")
        p_vals = [_req_int(x, f"payloads[{i}]", 0, (1 << w) - 1) for i, x in enumerate(payloads)]
        seen_addr: set[int] = set()
        valid_perm = True
        for i, rec_raw in enumerate(records):
            rec = _req_dict(rec_raw, f"records[{i}]")
            ai = _req_int(rec.get("address"), f"records[{i}].address", 0, n - 1)
            qi = _req_int(rec.get("payload"), f"records[{i}].payload", 0, (1 << w) - 1)
            if ai in seen_addr or qi != p_vals[ai]:
                valid_perm = False
            seen_addr.add(ai)
        if len(seen_addr) != n:
            valid_perm = False
        exp_entry = _req_int(cert.get("expected_tag_entry"), "expected_tag_entry", 0, 1)
        actual_entry = 1 if valid_perm else 0
        if actual_entry != exp_entry:
            raise VerificationError(f"Tag(h,w) entry {actual_entry} != expected {exp_entry}")

        # Lemma 8.2 substitution gap
        k_dim = _req_int(cert.get("flat_dim_k"), "flat_dim_k", 1, 16)
        b_pairs = [_req_int(x, f"catalyst_pairs[{i}]", 1, 10**6) for i, x in enumerate(_req_list(cert.get("catalyst_pairs"), "catalyst_pairs"))]
        if len(b_pairs) != k_dim:
            raise VerificationError("catalyst_pairs length must equal flat_dim_k")
        h_x = 2 * sum(b_pairs)
        exp_hx = _req_int(cert.get("expected_frame_price_hX"), "expected_frame_price_hX", 1)
        if h_x != exp_hx:
            raise VerificationError(f"h_X={h_x} != expected {exp_hx}")
        base_gap = _req_frac(cert.get("base_gap_D"), "base_gap_D")
        sub_gaps = _req_frac_vec(cert.get("occurrence_gaps"), "occurrence_gaps")
        if len(sub_gaps) == 0:
            raise VerificationError("occurrence_gaps must be non-empty")
        for idx, dg in enumerate(sub_gaps):
            if dg < 0:
                raise VerificationError(f"occurrence_gaps[{idx}] < 0")
        m_occ = len(sub_gaps)
        min_gap = base_gap
        for dg in sub_gaps:
            if dg < min_gap:
                min_gap = dg
        result_gap = min_gap - Fraction(2 * k_dim * m_occ, 1)
        exp_res_gap = _req_frac(cert.get("expected_substituted_gap"), "expected_substituted_gap")
        if result_gap != exp_res_gap:
            raise VerificationError(f"Lemma 8.2 substituted gap {result_gap} != {exp_res_gap}")
        if result_gap <= 0:
            raise VerificationError("Substituted gap must remain strictly positive for non-vacuity")

        return {
            "mode": "catalyst_gap_routing_budget",
            "submode": submode,
            "n_records": n,
            "tag_entry": actual_entry,
            "frame_price_hX": h_x,
            "substituted_gap": [result_gap.numerator, result_gap.denominator],
        }

    if submode == "routing_budget_and_bregman":
        # Eq. (9.1) multidirectional routing budget + Prop. 9.3 Eq. (9.2) Bregman potential decrease
        n = _req_int(cert.get("n"), "n", 1, 10**6)
        w = _req_int(cert.get("w"), "w", 1, 10**6)
        m_len = _req_int(cert.get("M"), "M", 2, 10**4)
        n_c = _req_int(cert.get("N_c"), "N_c", 1, 10**9)
        links = _req_list(cert.get("links"), "links", max_len=MAX_LIST_LEN)
        if len(links) == 0:
            raise VerificationError("links must be non-empty")
        total_charge = Fraction(0, 1)
        for idx, lk_raw in enumerate(links):
            lk = _req_dict(lk_raw, f"links[{idx}]")
            lam = _req_frac(lk.get("lambda"), f"links[{idx}].lambda")
            y = _req_int(lk.get("y"), f"links[{idx}].y", 0, 10**6)
            tp = _req_frac(lk.get("t_pow_p"), f"links[{idx}].t_pow_p")
            if lam < 0 or tp < 0:
                raise VerificationError(f"links[{idx}] negative lambda or t_pow_p")
            # Def. 9.1 item 3: at positive scale (t_pow_p > 0), length charge must equal 1
            if tp > 0 and lam != 1:
                raise VerificationError(f"links[{idx}] violates Def. 9.1: positive scale requires lambda = 1")
            total_charge = total_charge + lam * y + tp
        budget = Fraction(10 * n, 1) + Fraction(2 * m_len * n_c, w)
        if total_charge > budget:
            raise VerificationError(f"Eq. (9.1) budget exceeded: {total_charge} > {budget}")

        # Prop. 9.3 quadratic strongly convex score f(x) = (mu / 2) * (1/n_comm) * sum x_i^2
        # with box [0, a_i], inventory z, time s in (0, 1], x = x_s(z), and consistent x'
        s_time = _req_frac(cert.get("s"), "s")
        mu = _req_frac(cert.get("mu"), "mu")
        p_fee = _req_frac(cert.get("probe_fee_P"), "probe_fee_P")
        if s_time <= 0 or s_time > 1 or mu <= 0 or p_fee <= 0:
            raise VerificationError("Require s in (0,1], mu > 0, probe_fee_P > 0")
        a_box = _req_frac_vec(cert.get("box_a"), "box_a")
        m_dim = len(a_box)
        if m_dim == 0:
            raise VerificationError("box_a must be non-empty")
        z_vec = _req_frac_vec(cert.get("inventory_z"), "inventory_z", expected_len=m_dim)
        x_vec = _req_frac_vec(cert.get("optimizer_x"), "optimizer_x", expected_len=m_dim)
        xp_vec = _req_frac_vec(cert.get("probed_x_prime"), "probed_x_prime", expected_len=m_dim)

        for i in range(m_dim):
            if a_box[i] <= 0 or z_vec[i] < 0:
                raise VerificationError(f"Invalid box_a[{i}] or inventory_z[{i}]")
            if not (0 <= x_vec[i] <= a_box[i]) or not (0 <= xp_vec[i] <= a_box[i]):
                raise VerificationError(f"x[{i}] or x'[{i}] outside [0, a_{i}]")
            # Exact KKT condition for maximizing <z, x> - s * (mu/2) <x, x> on [0, a]:
            # x_i = min(a_i, z_i / (s * mu))
            unclamped = z_vec[i] * Fraction(s_time.denominator, s_time.numerator) * Fraction(mu.denominator, mu.numerator)
            opt_i = a_box[i] if unclamped >= a_box[i] else unclamped
            if x_vec[i] != opt_i:
                raise VerificationError(f"optimizer_x[{i}]={x_vec[i]} != exact KKT argmax {opt_i}")

        inv_m = Fraction(1, m_dim)
        # Bregman divergence for f(x) = (mu/2) * (1/m) * sum x_i^2 is exactly (mu/2) * (1/m) * sum (x'_i - x_i)^2
        breg = Fraction(1, 2) * mu * inv_m * sum(((xp_vec[i] - x_vec[i]) ** 2 for i in range(m_dim)), Fraction(0, 1))
        margin = Fraction(1, 2) * s_time * breg - p_fee
        if margin <= 0:
            raise VerificationError(f"Prop. 9.3 Eq. (9.2) Bregman margin not strictly positive: {margin}")

        return {
            "mode": "catalyst_gap_routing_budget",
            "submode": submode,
            "routing_charge": [total_charge.numerator, total_charge.denominator],
            "routing_budget": [budget.numerator, budget.denominator],
            "bregman_margin": [margin.numerator, margin.denominator],
        }

    # submode == "diagonal_contradiction_gap" (Sec. 16.1 Eq. 16.1)
    k_dim = _req_int(cert.get("flat_dim_k"), "flat_dim_k", 1, 16)
    g = _req_int(cert.get("grid_cells_g"), "grid_cells_g", 2, 10**4)
    r_scale = _req_int(cert.get("scale_R"), "scale_R", 8 * g, 10**9)
    e_pairs = _req_int(cert.get("external_pairs_E"), "external_pairs_E", 1, 10**9)
    d_gap = _req_frac(cert.get("cell_gap_D"), "cell_gap_D")
    c_k = _req_frac(cert.get("grid_const_Ck"), "grid_const_Ck")
    ell = _req_int(cert.get("selector_len_ell"), "selector_len_ell", 1, 10**9)
    if d_gap <= 0 or c_k <= 0:
        raise VerificationError("cell_gap_D and grid_const_Ck must be positive")
    h_id = Fraction(2 * r_scale * e_pairs, 1)
    grid_err = c_k * Fraction(e_pairs * g * g, r_scale)
    sel_err = Fraction(2 * k_dim * r_scale * (ell + g**k_dim), 1)
    upper = h_id - d_gap + grid_err + sel_err
    if upper < 0:
        raise VerificationError(f"Eq. (16.1) upper bound {upper} < 0 (out of spectral price regime; require D <= H(I))")
    if not (upper < h_id):
        raise VerificationError(f"Eq. (16.1) strict contradiction bound failed: {upper} >= {h_id}")
    net_saving = h_id - upper
    return {
        "mode": "catalyst_gap_routing_budget",
        "submode": submode,
        "H_identity": [h_id.numerator, h_id.denominator],
        "upper_bound": [upper.numerator, upper.denominator],
        "net_saving": [net_saving.numerator, net_saving.denominator],
    }


def verify_certificate(cert: Any) -> dict[str, Any]:
    d = _req_dict(cert, "certificate")
    schema = _req_str(d.get("schema"), "schema")
    if schema != SCHEMA_VERSION:
        raise VerificationError(f"Unsupported schema {schema!r}; expected {SCHEMA_VERSION!r}")
    cert_id = _req_str(d.get("certificate_id"), "certificate_id")
    mode = _req_str(d.get("mode"), "mode")
    _req_str(d.get("claim_ref"), "claim_ref")

    if mode == "flattening_cut_nonadditivity":
        res = _verify_mode1(d)
    elif mode == "exterior_crossing_fermion":
        res = _verify_mode2(d)
    elif mode == "cat0_flat_mean_width":
        res = _verify_mode3(d)
    elif mode == "catalyst_gap_routing_budget":
        res = _verify_mode4(d)
    else:
        raise VerificationError(f"Unknown mode: {mode!r}")

    res["certificate_id"] = cert_id
    res["schema"] = schema
    res["status"] = "VERIFIED"
    return res


def verify_file(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    if len(raw) > MAX_JSON_BYTES:
        raise VerificationError(f"Certificate file exceeds MAX_JSON_BYTES ({MAX_JSON_BYTES})")
    try:
        text = raw.decode("utf-8")
        obj = json.loads(text, object_pairs_hook=_no_dup_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"Malformed JSON: {exc}") from exc
    return verify_certificate(obj)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify TENSOR-SPEC-CERT-1.0 certificates.")
    parser.add_argument("paths", nargs="+", type=Path, help="Certificate JSON file(s) or directory")
    args = parser.parse_args(argv)

    files: list[Path] = []
    for p in args.paths:
        if p.is_dir():
            files.extend(sorted(p.glob("*.json")))
        else:
            files.append(p)
    if not files:
        print("No certificate files found.", file=sys.stderr)
        return 2

    ok = 0
    for f in files:
        try:
            summary = verify_file(f)
            print(f"[VERIFIED] {f.name}: {summary['mode']} ({summary.get('submode', 'core')})")
            ok += 1
        except VerificationError as exc:
            print(f"[REJECTED] {f.name}: {exc}", file=sys.stderr)
            return 1
    print(f"Verified {ok}/{len(files)} certificates under {SCHEMA_VERSION}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
