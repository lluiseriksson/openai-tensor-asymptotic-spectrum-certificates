"""Generate the 16 canonical TENSOR-SPEC-CERT-1.0 certificates for Paper 68."""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.tensor_spectrum_verifier import SCHEMA_VERSION, _exterior_ops, verify_certificate


def _build_explicit_crossing(q: int, b: int) -> tuple[int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[list[int]]]:
    """Build explicit integer crossing matrices (A+, B+, A-, B-) from exterior algebra Lambda^b Z^q."""
    b_cur, b_plus, b_minus, _, _ = _exterior_ops(q, b)
    r_val = len(b_cur)
    rp_val = len(b_plus)
    rm_val = len(b_minus)
    idx_plus = {s: idx for idx, s in enumerate(b_plus)}
    idx_minus = {s: idx for idx, s in enumerate(b_minus)}
    dim = q * r_val

    # A_plus[(i, a), t] = coeff of e_t in L_i(e_a)
    # B_plus[t, (j, b)] = coeff of e_b in D_j(e_t)
    a_plus = [[0] * rp_val for _ in range(dim)]
    b_plus_m = [[0] * dim for _ in range(rp_val)]
    for i in range(q):
        for a_idx, s in enumerate(b_cur):
            if i not in s:
                pos = 0
                while pos < len(s) and s[pos] < i:
                    pos += 1
                ns = s[:pos] + (i,) + s[pos:]
                a_plus[i * r_val + a_idx][idx_plus[ns]] = 1 if (pos % 2 == 0) else -1
    for t_idx, t in enumerate(b_plus):
        for j in range(q):
            if j in t:
                pos = t.index(j)
                ns = t[:pos] + t[pos + 1 :]
                b_idx = b_cur.index(ns)
                b_plus_m[t_idx][j * r_val + b_idx] = 1 if (pos % 2 == 0) else -1

    # A_minus[(i, b), m] = coeff of e_b in L_i(e_m)
    # B_minus[m, (j, a)] = coeff of e_m in D_j(e_a)
    a_minus = [[0] * rm_val for _ in range(dim)]
    b_minus_m = [[0] * dim for _ in range(rm_val)]
    for i in range(q):
        for m_idx, ms in enumerate(b_minus):
            if i not in ms:
                pos = 0
                while pos < len(ms) and ms[pos] < i:
                    pos += 1
                ns = ms[:pos] + (i,) + ms[pos:]
                b_idx = b_cur.index(ns)
                a_minus[i * r_val + b_idx][m_idx] = 1 if (pos % 2 == 0) else -1
    for j in range(q):
        for a_idx, s in enumerate(b_cur):
            if j in s:
                pos = s.index(j)
                ns = s[:pos] + s[pos + 1 :]
                b_minus_m[idx_minus[ns]][j * r_val + a_idx] = 1 if (pos % 2 == 0) else -1

    return q, r_val, rp_val, rm_val, a_plus, b_plus_m, a_minus, b_minus_m


def build_certificates() -> list[tuple[str, dict]]:
    certs: list[tuple[str, dict]] = []

    # ---- MODE 1: Flattening cut semiring & Sec. 1.1 geometric-mean non-additivity (C01-C04) ----
    # C01: Canonical Sec. 1.1 square example (Q): U=(2,1), V=(1,2), U+V=(3,3), G^2=9 > 8, Hamming square violates Eq. (4.5)
    certs.append((
        "C01_sec11_square_nonadditivity_Q.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C01",
            "mode": "flattening_cut_nonadditivity",
            "claim_ref": "Sec. 1.1 & Eq. (4.5): square midline geometric mean 3 > 2*sqrt(2) and Hamming square 4-point violation over Q",
            "field_char": 0,
            "certify_nonadditivity": True,
            "cuts": [
                {
                    "cut_label": "vertical_midline",
                    "matrix_U": [[1, 0], [0, 1]],
                    "matrix_V": [[1, 1], [1, 1]],
                    "expected_rank_U": 2,
                    "expected_rank_V": 1,
                },
                {
                    "cut_label": "horizontal_midline",
                    "matrix_U": [[1, 1], [1, 1]],
                    "matrix_V": [[1, 0], [0, 1]],
                    "expected_rank_U": 1,
                    "expected_rank_V": 2,
                },
            ],
            "four_point_test": {
                "d_AB": [1, 1],
                "d_BC": [1, 1],
                "d_CD": [1, 1],
                "d_DA": [1, 1],
                "d_AC": [2, 1],
                "d_BD": [2, 1],
                "expect_cat0_violation": True,
            },
        },
    ))

    # C02: Characteristic-2 (F_2) verification of Sec. 1.1 square non-additivity (field-independent obstruction)
    certs.append((
        "C02_sec11_square_nonadditivity_F2.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C02",
            "mode": "flattening_cut_nonadditivity",
            "claim_ref": "Sec. 1.1 over F_2: characteristic-2 invariance of midline cut non-additivity",
            "field_char": 2,
            "certify_nonadditivity": True,
            "cuts": [
                {
                    "cut_label": "vertical_midline_F2",
                    "matrix_U": [[1, 0], [0, 1]],
                    "matrix_V": [[1, 0], [0, 0]],
                    "expected_rank_U": 2,
                    "expected_rank_V": 1,
                },
                {
                    "cut_label": "horizontal_midline_F2",
                    "matrix_U": [[1, 0], [0, 0]],
                    "matrix_V": [[1, 0], [0, 1]],
                    "expected_rank_U": 1,
                    "expected_rank_V": 2,
                },
            ],
            "four_point_test": {
                "d_AB": [1, 2],
                "d_BC": [1, 2],
                "d_CD": [1, 2],
                "d_DA": [1, 2],
                "d_AC": [1, 1],
                "d_BD": [1, 1],
                "expect_cat0_violation": True,
            },
        },
    ))

    # C03: 3-cut cubic axis-plane obstruction (U=(4,1,1), V=(1,4,1), W-sum=(5,5,2), prod_sum=50 > (4^(1/3)+4^(1/3))^3 = 32)
    certs.append((
        "C03_3cut_cube_nonadditivity_Q.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C03",
            "mode": "flattening_cut_nonadditivity",
            "claim_ref": "Sec. 1.1 3-cut extension: geometric mean over 3 orthogonal cuts gives 50 > (8/5 + 8/5)^3 = 4096/125",
            "field_char": 0,
            "certify_nonadditivity": True,
            "geom_upper_U": [8, 5],
            "geom_upper_V": [8, 5],
            "cuts": [
                {
                    "cut_label": "x_cut",
                    "matrix_U": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
                    "matrix_V": [[1, 0], [0, 0]],
                    "expected_rank_U": 4,
                    "expected_rank_V": 1,
                },
                {
                    "cut_label": "y_cut",
                    "matrix_U": [[1, 0], [0, 0]],
                    "matrix_V": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
                    "expected_rank_U": 1,
                    "expected_rank_V": 4,
                },
                {
                    "cut_label": "z_cut",
                    "matrix_U": [[1, 0], [0, 0]],
                    "matrix_V": [[1, 0], [0, 0]],
                    "expected_rank_U": 1,
                    "expected_rank_V": 1,
                },
            ],
            "four_point_test": {
                "d_AB": [3, 1],
                "d_BC": [4, 1],
                "d_CD": [3, 1],
                "d_DA": [4, 1],
                "d_AC": [5, 1],
                "d_BD": [5, 1],
                "expect_cat0_violation": False,
            },
        },
    ))

    # C04: Proportional-cut control (additive equality case over F_5, Euclidean 3x4 rectangle passes Eq. 4.5)
    certs.append((
        "C04_proportional_cuts_additive_F5.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C04",
            "mode": "flattening_cut_nonadditivity",
            "claim_ref": "Sec. 2.1 proportional flattening control over F_5: U=(2,2), V=(3,3) preserves additivity",
            "field_char": 5,
            "certify_nonadditivity": False,
            "cuts": [
                {
                    "cut_label": "cut_1",
                    "matrix_U": [[1, 2], [3, 4]],
                    "matrix_V": [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
                    "expected_rank_U": 2,
                    "expected_rank_V": 3,
                },
                {
                    "cut_label": "cut_2",
                    "matrix_U": [[2, 1], [1, 2]],
                    "matrix_V": [[1, 1, 0], [0, 1, 1], [0, 0, 1]],
                    "expected_rank_U": 2,
                    "expected_rank_V": 3,
                },
            ],
            "four_point_test": {
                "d_AB": [3, 1],
                "d_BC": [4, 1],
                "d_CD": [3, 1],
                "d_DA": [4, 1],
                "d_AC": [5, 1],
                "d_BD": [5, 1],
                "expect_cat0_violation": False,
            },
        },
    ))

    # ---- MODE 2: Exterior algebra Eq. (3.2) & Theorem 4.1 Shifted Crossing / Unimodularity (C05-C08) ----
    # C05: Eq. (3.2) exterior crossing on Lambda^1 Z^2 (Q=2, b=1 -> R=2, R+=1, R-=1)
    certs.append((
        "C05_exterior_crossing_Q2_b1.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C05",
            "mode": "exterior_crossing_fermion",
            "submode": "exterior_crossing",
            "claim_ref": "Eq. (3.2) & Eq. (3.3) on Lambda^1 Z^2: symmetric rank-1 + rank-1 crossing of E(2) x E(2)",
            "Q": 2,
            "b": 1,
            "expected_R": 2,
            "expected_R_plus": 1,
            "expected_R_minus": 1,
        },
    ))

    # C06: Eq. (3.2) shifted asymmetric crossing on Lambda^1 Z^4 (Q=4, b=1 -> R=4, R+=6, R-=1)
    certs.append((
        "C06_exterior_shifted_Q4_b1.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C06",
            "mode": "exterior_crossing_fermion",
            "submode": "exterior_crossing",
            "claim_ref": "Eq. (3.2) & Theorem 4.1 asymmetric shift on Lambda^1 Z^4: R=4, R+=6, R-=1 over Z",
            "Q": 4,
            "b": 1,
            "expected_R": 4,
            "expected_R_plus": 6,
            "expected_R_minus": 1,
        },
    ))

    # C07: Proof of Theorem 4.1 (p. 19) Vandermonde & Cauchy-Binet unimodularity on energy downset H=9/2, m=2
    # m=2 subsets of {0,1,2,3,4} with sum(2i+1) <= 9:
    # {0,1}: 1+3=4, {0,2}: 1+5=6, {0,3}: 1+7=8, {1,2}: 3+5=8 -> 4 states!
    certs.append((
        "C07_vandermonde_unimodular_m2_H9half.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C07",
            "mode": "exterior_crossing_fermion",
            "submode": "vandermonde_unimodular",
            "claim_ref": "Theorem 4.1 proof (p. 19): unit lower triangular P^[2] and det(C^[2])=+-1 for H=9/2 (4 states)",
            "two_H": 9,
            "block_size_m": 2,
            "expected_num_states": 4,
        },
    ))

    # C08: Explicit integer crossing factor matrices for d_X=3, d_Y=3 with r+=3, r-=1 (from Lambda^1 Z^3)
    dx, dy, rp, rm, ap, bp, am, bm = _build_explicit_crossing(3, 1)
    certs.append((
        "C08_explicit_crossing_factors_3x3.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C08",
            "mode": "exterior_crossing_fermion",
            "submode": "explicit_crossing_factors",
            "claim_ref": "Theorem 4.1 explicit integer factor matrices A+,B+,A-,B- realizing E(3) x E(3) = T+ + T- (r+=3, r-=1)",
            "d_X": dx,
            "d_Y": dy,
            "r_plus": rp,
            "r_minus": rm,
            "A_plus": ap,
            "B_plus": bp,
            "A_minus": am,
            "B_minus": bm,
        },
    ))

    # ---- MODE 3: CAT(0), Norm Rigidity, Mean-Width Law & Prop. 5.2 Copying (C09-C12) ----
    # C09: Euclidean 3x4 rectangle + center midpoint (5 sites: (0,0),(3,0),(3,4),(0,4),(3/2,2))
    # All pairwise distances are exact rationals! (sides 3,4; diags 5; corner-to-mid 5/2!)
    d_rect = [
        [[0, 1], [3, 1], [5, 1], [4, 1], [5, 2]],
        [[3, 1], [0, 1], [4, 1], [5, 1], [5, 2]],
        [[5, 1], [4, 1], [0, 1], [3, 1], [5, 2]],
        [[4, 1], [5, 1], [3, 1], [0, 1], [5, 2]],
        [[5, 2], [5, 2], [5, 2], [5, 2], [0, 1]],
    ]
    certs.append((
        "C09_cat0_euclidean_rectangle_equality.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C09",
            "mode": "cat0_flat_mean_width",
            "submode": "cat0_quadrilateral_midpoint",
            "claim_ref": "Lemma 4.4 & Sec. 4.4: 3x4 Euclidean rectangle + midpoint achieves exact CAT(0) Eq. (4.4) and Eq. (4.5) equality",
            "num_sites": 5,
            "distance_matrix": d_rect,
            "geodesic_checks": [
                {"x": 0, "y": 2, "m": 4, "z": 1, "s": [1, 2], "expect_equality": True},
                {"x": 1, "y": 3, "m": 4, "z": 0, "s": [1, 2], "expect_equality": True},
            ],
        },
    ))

    # C10: Strict CAT(0) metric tree (3-ray tripod with center 0, leaves 1,2,3 at dist 2, midpoint 4 of [1,2] is 0!)
    # Wait: for tripod with center 0 and arms of length 2 to leaves 1, 2, 3:
    # d(1,2)=4, midpoint of [1,2] is 0 (d(1,0)=2, d(0,2)=2).
    # Eq. (4.4) for z=3: d(3,0)^2 = 4 < 1/2*16 + 1/2*16 - 1/4*16 = 12 (strictly negatively curved!)
    d_tree = [
        [[0, 1], [2, 1], [2, 1], [2, 1]],
        [[2, 1], [0, 1], [4, 1], [4, 1]],
        [[2, 1], [4, 1], [0, 1], [4, 1]],
        [[2, 1], [4, 1], [4, 1], [0, 1]],
    ]
    certs.append((
        "C10_cat0_tripod_strict_curvature.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C10",
            "mode": "cat0_flat_mean_width",
            "submode": "cat0_quadrilateral_midpoint",
            "claim_ref": "Theorem 4.3 & Eq. (4.4): metric tripod satisfies strict CAT(0) curvature (4 < 12) and Eq. (4.5)",
            "num_sites": 4,
            "distance_matrix": d_tree,
            "geodesic_checks": [
                {"x": 1, "y": 2, "m": 0, "z": 3, "s": [1, 2], "expect_equality": False},
                {"x": 1, "y": 3, "m": 0, "z": 2, "s": [1, 2], "expect_equality": False},
            ],
        },
    ))

    # C11: Lemma 4.2 Norm rigidity violation for l_1 plane (dual norm l_inf) with f=(1,0), g=(0,1), t=1, u=1 -> S=2, A+=1, A-=1? Wait:
    # In l_inf dual norm: f=(1, 1), g=(1, -1) are on l_inf unit sphere (max(|1|,|1|)=1)!
    # With t=1, u=1: S=2, tf+ug = (2, 0) -> A_+=2, A_+^2=4; tf-ug = (0, 2) -> A_-=2, A_-^2=4!
    # A_+^2 + A_-^2 = 8 > 4 = 2*S! And S - A_-^2 = -2 < delta = 0 < 2 = A_+^2 - S!
    certs.append((
        "C11_norm_rigidity_l1_violation.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C11",
            "mode": "cat0_flat_mean_width",
            "submode": "norm_rigidity_violation",
            "claim_ref": "Lemma 4.2 Eq. (4.2): l_1 plane dual norm (l_inf) violates parallelogram law 4+4=8 > 4 with shift delta=0",
            "dual_norm": "l_inf",
            "t": [1, 1],
            "u": [1, 1],
            "f": [[1, 1], [1, 1]],
            "g": [[1, 1], [-1, 1]],
            "delta": [0, 1],
        },
    ))

    # C12: Theorem 5.3 planar polygon mean width + Prop. 5.2 submodularity + convex-hull zero-cost insertion
    # 3x4 rectangle: verts (0,0), (3,0), (3,4), (0,4). Edges 3,4,3,4 -> per=14 -> W=7.
    # Interior z=(3/2, 2) with weights (1/4, 1/4, 1/4, 1/4) preserves delta=7.
    # Submodularity: A={0,1,2} (3-4-5 triangle, W=6), B={0,2,3} (3-4-5 triangle, W=6),
    # A cap B={0,2} (diagonal len 5, W=5), A cup B={0,1,2,3} (W=7): 7 + 5 = 12 <= 6 + 6 = 12!
    certs.append((
        "C12_mean_width_rectangle_submodularity.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C12",
            "mode": "cat0_flat_mean_width",
            "submode": "polygon_mean_width_copying",
            "claim_ref": "Theorem 5.3 & Prop. 5.2: 3x4 rectangle W=7, interior copy z=(3/2,2), tight submodularity 7+5 <= 6+6",
            "polygon_vertices": [
                [[0, 1], [0, 1]],
                [[3, 1], [0, 1]],
                [[3, 1], [4, 1]],
                [[0, 1], [4, 1]],
            ],
            "edge_lengths": [[3, 1], [4, 1], [3, 1], [4, 1]],
            "expected_mean_width": [7, 1],
            "interior_point": {
                "z": [[3, 2], [2, 1]],
                "convex_weights": [[1, 4], [1, 4], [1, 4], [1, 4]],
                "delta_with_z": [7, 1],
            },
            "submodularity_check": {
                "indices_A": [0, 1, 2],
                "indices_B": [0, 2, 3],
                "delta_A": [6, 1],
                "delta_B": [6, 1],
                "delta_A_inter_B": [5, 1],
                "delta_A_union_B": [7, 1],
            },
        },
    ))

    # ---- MODE 4: Prop. 5.2 Copying Degeneration + Parts II-III Catalyst/Routing/Diagonal (C13-C16) ----
    # C13: Prop. 5.2 quadratic-weight two-stage degeneration (r=4, q=5, h=2, v=[2,1,2,3] -> g=2 survives stage 1, killed at stage 2 with w2=2)
    certs.append((
        "C13_prop52_quadratic_degeneration.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C13",
            "mode": "cat0_flat_mean_width",
            "submode": "copying_quadratic_degeneration",
            "claim_ref": "Prop. 5.2 (p. 24): stage-1 survival (g=h=2) and stage-2 positive quadratic weight sum (v_j - v_1)^2 = 2",
            "r": 4,
            "q": 5,
            "h": 2,
            "v": [2, 1, 2, 3],
            "g": 2,
            "expected_L": 33,
            "survives_first_degeneration": True,
            "survives_second_degeneration": False,
        },
    ))

    # C14: Def. 8.1 Tag(2, 4) permutation witness + Lemma 8.2 substitution gap
    certs.append((
        "C14_tag_perm_and_catalyst_substitution.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C14",
            "mode": "catalyst_gap_routing_budget",
            "submode": "tag_perm_and_catalyst_gap",
            "claim_ref": "Def. 8.1 Tag(2,4) indicator = 1 & Lemma 8.2 substitution gap D_sub = 30 - 2*3*2 = 18 on 3-flat",
            "h": 2,
            "w": 4,
            "payloads": [11, 5, 14, 2],
            "records": [
                {"address": 2, "payload": 14},
                {"address": 0, "payload": 11},
                {"address": 3, "payload": 2},
                {"address": 1, "payload": 5},
            ],
            "expected_tag_entry": 1,
            "flat_dim_k": 3,
            "catalyst_pairs": [10, 12, 8],
            "expected_frame_price_hX": 60,
            "base_gap_D": [36, 1],
            "occurrence_gaps": [[30, 1], [34, 1]],
            "expected_substituted_gap": [18, 1],
        },
    ))

    # C15: Def. 9.1 Eq. (9.1) routing budget & Prop. 9.3 Eq. (9.2) adaptive probing Bregman decrease
    # n=4, w=8, M=5, N_c=64 -> budget = 40 + 2*5*64/8 = 40 + 80 = 120.
    # Def. 9.1 item 3: inter-layer links (lambda=5) have scale zero (t_pow_p=0); positive-scale links have lambda=1!
    # Links: 3 links totaling 5*10 + 0 + 5*10 + 0 + 1*12 + 1/2 = 225/2 = 112.5 <= 120.
    # Prop. 9.3: s=1/2, mu=2, a=[2,2,2,2], z=[1,2,3,4] -> x_i = min(2, z_i / 1) = [1,2,2,2].
    # Probed x' = [0,1,1,1] -> diff = [-1,-1,-1,-1], sum sq = 4, breg = (1/2)*2*(1/4)*4 = 1.
    # (s/2)*breg - P = (1/4)*1 - 1/8 = 1/8 > 0!
    certs.append((
        "C15_routing_budget_and_bregman_probing.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C15",
            "mode": "catalyst_gap_routing_budget",
            "submode": "routing_budget_and_bregman",
            "claim_ref": "Def. 9.1 Eq. (9.1) budget 225/2 <= 120 & Prop. 9.3 Eq. (9.2) Bregman margin 1/8 > 0",
            "n": 4,
            "w": 8,
            "M": 5,
            "N_c": 64,
            "links": [
                {"lambda": [5, 1], "y": 10, "t_pow_p": [0, 1]},
                {"lambda": [5, 1], "y": 10, "t_pow_p": [0, 1]},
                {"lambda": [1, 1], "y": 12, "t_pow_p": [1, 2]},
            ],
            "s": [1, 2],
            "mu": [2, 1],
            "probe_fee_P": [1, 8],
            "box_a": [[2, 1], [2, 1], [2, 1], [2, 1]],
            "inventory_z": [[1, 1], [2, 1], [3, 1], [4, 1]],
            "optimizer_x": [[1, 1], [2, 1], [2, 1], [2, 1]],
            "probed_x_prime": [[0, 1], [1, 1], [1, 1], [1, 1]],
        },
    ))

    # C16: Sec. 16.1 Eq. (16.1) in-regime diagonal contradiction gap (0 < upper < H(I))
    # k=2, g=4, R=64 (>= 8*g=32), E=1000, D=50000, C_k=2, ell=50
    # H(I) = 2*64*1000 = 128000
    # grid_err = 2 * (1000 * 16) / 64 = 500
    # sel_err = 2 * 2 * 64 * (50 + 16) = 256 * 66 = 16896
    # upper = 128000 - 50000 + 500 + 16896 = 95396 (0 < 95396 < 128000, net saving 32604!)
    certs.append((
        "C16_sec16_diagonal_contradiction_gap.json",
        {
            "schema": SCHEMA_VERSION,
            "certificate_id": "C16",
            "mode": "catalyst_gap_routing_budget",
            "submode": "diagonal_contradiction_gap",
            "claim_ref": "Sec. 16.1 Eq. (16.1): in-regime gap 0 < 95396 < 128000 = H(I) (saving 32604) contradicts Z - Z_hat = I",
            "flat_dim_k": 2,
            "grid_cells_g": 4,
            "scale_R": 64,
            "external_pairs_E": 1000,
            "cell_gap_D": [50000, 1],
            "grid_const_Ck": [2, 1],
            "selector_len_ell": 50,
        },
    ))

    return certs


def main() -> int:
    out_dir = ROOT / "certificates"
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.json"):
        old.unlink()
    certs = build_certificates()
    for fname, cert in certs:
        res = verify_certificate(cert)
        (out_dir / fname).write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8")
        print(f"Generated & verified {fname}: {res['mode']} ({res.get('submode', 'core')})")
    print(f"Total certificates generated: {len(certs)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
