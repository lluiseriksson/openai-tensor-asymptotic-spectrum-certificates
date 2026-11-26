"""Comprehensive pytest & adversarial forgery suite for TENSOR-SPEC-CERT-1.0."""

from __future__ import annotations

import ast
import copy
import json
from pathlib import Path
import sys
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.tensor_spectrum_verifier import (
    MAX_JSON_BYTES,
    SCHEMA_VERSION,
    VerificationError,
    verify_certificate,
    verify_file,
)

CERT_DIR = ROOT / "certificates"


def _load(name: str) -> dict:
    return json.loads((CERT_DIR / name).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 1. Static AST Hygiene & All 16 Canonical Certificates
# ---------------------------------------------------------------------------
def test_zero_float_and_no_div_ast_hygiene() -> None:
    src_path = ROOT / "src" / "tensor_spectrum_verifier.py"
    tree = ast.parse(src_path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        assert not isinstance(node, ast.Div), "True division (/) forbidden in exact verifier"
        if isinstance(node, ast.AugAssign):
            assert not isinstance(node.op, ast.Div), "/= forbidden in exact verifier"
        if isinstance(node, ast.Constant):
            assert not isinstance(node.value, float), "Float literal forbidden"
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert alias.name != "math", "math module forbidden"
        if isinstance(node, ast.ImportFrom):
            assert node.module != "math", "math module forbidden"
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            assert node.func.id != "float", "float() call forbidden"


@pytest.mark.parametrize("cert_path", sorted(CERT_DIR.glob("*.json")), ids=lambda p: p.name)
def test_all_canonical_certificates_verify(cert_path: Path) -> None:
    res = verify_file(cert_path)
    assert res["status"] == "VERIFIED"
    assert res["schema"] == SCHEMA_VERSION


def test_exact_certificate_count() -> None:
    assert len(list(CERT_DIR.glob("*.json"))) == 16


# ---------------------------------------------------------------------------
# 2. Structural / JSON / Rational / DoS Forgeries
# ---------------------------------------------------------------------------
def test_reject_wrong_schema() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["schema"] = "WRONG-1.0"
    with pytest.raises(VerificationError, match="Unsupported schema"):
        verify_certificate(c)


def test_reject_unknown_mode() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["mode"] = "bogus_mode"
    with pytest.raises(VerificationError, match="Unknown mode"):
        verify_certificate(c)


def test_reject_duplicate_json_key(tmp_path: Path) -> None:
    p = tmp_path / "dup.json"
    p.write_text('{"schema":"TENSOR-SPEC-CERT-1.0","schema":"TENSOR-SPEC-CERT-1.0"}', encoding="utf-8")
    with pytest.raises(VerificationError, match="Duplicate JSON key"):
        verify_file(p)


def test_reject_oversized_json(tmp_path: Path) -> None:
    p = tmp_path / "big.json"
    p.write_bytes(b" " * (MAX_JSON_BYTES + 1))
    with pytest.raises(VerificationError, match="MAX_JSON_BYTES"):
        verify_file(p)


def test_reject_non_lowest_terms_rational() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"]["d_AB"] = [2, 2]
    with pytest.raises(VerificationError, match="lowest terms"):
        verify_certificate(c)


def test_reject_zero_denominator() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"]["d_AB"] = [1, 0]
    with pytest.raises(VerificationError, match="< 1"):
        verify_certificate(c)


def test_reject_negative_denominator() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"]["d_AB"] = [1, -2]
    with pytest.raises(VerificationError, match="< 1"):
        verify_certificate(c)


def test_reject_bool_as_int() -> None:
    c = _load("C05_exterior_crossing_Q2_b1.json")
    c["Q"] = True
    with pytest.raises(VerificationError, match="must be an integer"):
        verify_certificate(c)


def test_reject_whitespace_string() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["claim_ref"] = "   "
    with pytest.raises(VerificationError, match="non-empty string"):
        verify_certificate(c)


def test_reject_huge_rational_numerator() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"]["d_AB"] = [10**13, 1]
    with pytest.raises(VerificationError, match="> 1000000000000"):
        verify_certificate(c)


# ---------------------------------------------------------------------------
# 3. Mode 1 Forgeries (Flattening Cuts & Sec. 1.1 Non-Additivity)
# ---------------------------------------------------------------------------
def test_mode1_reject_composite_char() -> None:
    c = _load("C02_sec11_square_nonadditivity_F2.json")
    c["field_char"] = 4
    with pytest.raises(VerificationError, match="prime"):
        verify_certificate(c)


def test_mode1_reject_single_cut() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["cuts"] = c["cuts"][:1]
    with pytest.raises(VerificationError, match="at least 2"):
        verify_certificate(c)


def test_mode1_reject_spoofed_rank_U() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["cuts"][0]["expected_rank_U"] = 1
    with pytest.raises(VerificationError, match="rank\\(U\\)"):
        verify_certificate(c)


def test_mode1_reject_spoofed_rank_V() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["cuts"][1]["expected_rank_V"] = 1
    with pytest.raises(VerificationError, match="rank\\(V\\)"):
        verify_certificate(c)


def test_mode1_reject_false_nonadditivity_claim() -> None:
    c = _load("C04_proportional_cuts_additive_F5.json")
    c["certify_nonadditivity"] = True
    with pytest.raises(VerificationError, match="non-additivity failed"):
        verify_certificate(c)


def test_mode1_reject_false_additivity_claim() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["certify_nonadditivity"] = False
    with pytest.raises(VerificationError, match="Expected additive"):
        verify_certificate(c)


def test_mode1_reject_bad_3cut_bracket() -> None:
    c = _load("C03_3cut_cube_nonadditivity_Q.json")
    c["geom_upper_U"] = [1, 1]
    with pytest.raises(VerificationError, match="do not dominate"):
        verify_certificate(c)


def test_mode1_reject_loose_3cut_bracket() -> None:
    c = _load("C03_3cut_cube_nonadditivity_Q.json")
    c["geom_upper_U"] = [3, 1]
    c["geom_upper_V"] = [3, 1]
    with pytest.raises(VerificationError, match="not strictly less"):
        verify_certificate(c)


def test_mode1_reject_triangle_violation_in_four_point() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"]["d_AC"] = [5, 1]
    with pytest.raises(VerificationError, match="triangle inequality"):
        verify_certificate(c)


def test_mode1_reject_spoofed_cat0_violation_flag() -> None:
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"]["expect_cat0_violation"] = False
    with pytest.raises(VerificationError, match="Unexpected Eq. \\(4.5\\)"):
        verify_certificate(c)


def test_mode1_reject_f2_rank_drop() -> None:
    c = _load("C02_sec11_square_nonadditivity_F2.json")
    c["cuts"][0]["matrix_U"] = [[2, 0], [0, 2]]
    with pytest.raises(VerificationError, match="rank\\(U\\)=0"):
        verify_certificate(c)


# ---------------------------------------------------------------------------
# 4. Mode 2 Forgeries (Exterior Algebra, Vandermonde, Crossing Factors)
# ---------------------------------------------------------------------------
def test_mode2_reject_bad_submode() -> None:
    c = _load("C05_exterior_crossing_Q2_b1.json")
    c["submode"] = "invalid"
    with pytest.raises(VerificationError, match="Unknown Mode 2 submode"):
        verify_certificate(c)


def test_mode2_reject_Q_too_large() -> None:
    c = _load("C05_exterior_crossing_Q2_b1.json")
    c["Q"] = 9
    with pytest.raises(VerificationError, match="> 6"):
        verify_certificate(c)


def test_mode2_reject_wrong_expected_R() -> None:
    c = _load("C06_exterior_shifted_Q4_b1.json")
    c["expected_R_plus"] = 5
    with pytest.raises(VerificationError, match="Exterior dimensions"):
        verify_certificate(c)


def test_mode2_reject_b_out_of_range() -> None:
    c = _load("C05_exterior_crossing_Q2_b1.json")
    c["b"] = 2
    with pytest.raises(VerificationError, match="> 1"):
        verify_certificate(c)


def test_mode2_reject_wrong_downset_size() -> None:
    c = _load("C07_vandermonde_unimodular_m2_H9half.json")
    c["expected_num_states"] = 5
    with pytest.raises(VerificationError, match="Energy downset size"):
        verify_certificate(c)


def test_mode2_reject_corrupted_crossing_factor() -> None:
    c = _load("C08_explicit_crossing_factors_3x3.json")
    c["A_plus"][0][0] = 7
    with pytest.raises(VerificationError, match="Crossing sum mismatch"):
        verify_certificate(c)


def test_mode2_reject_wrong_factor_dimensions() -> None:
    c = _load("C08_explicit_crossing_factors_3x3.json")
    c["A_minus"].pop()
    with pytest.raises(VerificationError, match="row count"):
        verify_certificate(c)


def test_mode2_reject_vacuous_zero_summand() -> None:
    # Identity on cut (x+,y+)|(x-,y-) has rank 4 -> rejected as trivial crossing split r_plus >= dim
    eye4 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
    c = {
        "schema": SCHEMA_VERSION,
        "certificate_id": "VAC",
        "mode": "exterior_crossing_fermion",
        "submode": "explicit_crossing_factors",
        "claim_ref": "vacuous test",
        "d_X": 2,
        "d_Y": 2,
        "r_plus": 4,
        "r_minus": 1,
        "A_plus": eye4,
        "B_plus": eye4,
        "A_minus": [[0], [0], [0], [0]],
        "B_minus": [[0, 0, 0, 0]],
    }
    with pytest.raises(VerificationError, match="> 3"):
        verify_certificate(c)


def test_mode2_reject_trivial_rank_sum_split() -> None:
    c = {
        "schema": SCHEMA_VERSION,
        "certificate_id": "TRIV",
        "mode": "exterior_crossing_fermion",
        "submode": "explicit_crossing_factors",
        "claim_ref": "trivial split r_+ + r_- = dim",
        "d_X": 2,
        "d_Y": 2,
        "r_plus": 2,
        "r_minus": 2,
        "A_plus": [[1, 0], [0, 1], [0, 0], [0, 0]],
        "B_plus": [[1, 0, 0, 0], [0, 1, 0, 0]],
        "A_minus": [[0, 0], [0, 0], [1, 0], [0, 1]],
        "B_minus": [[0, 0, 1, 0], [0, 0, 0, 1]],
    }
    with pytest.raises(VerificationError, match="Trivial crossing split"):
        verify_certificate(c)


def test_mode2_exterior_crossing_Q5_b2_verified_exact() -> None:
    # B5 check: Q=5, b=2 has dim_flat=50 > 32, must genuinely compute exact crossing ranks [10, 5]
    c = {
        "schema": SCHEMA_VERSION,
        "certificate_id": "Q5B2",
        "mode": "exterior_crossing_fermion",
        "submode": "exterior_crossing",
        "claim_ref": "exact crossing rank computation at dim_flat=50",
        "Q": 5,
        "b": 2,
        "expected_R": 10,
        "expected_R_plus": 10,
        "expected_R_minus": 5,
    }
    res = verify_certificate(c)
    assert res["verified_crossing_ranks"] == [10, 5]


def test_mode1_reject_side_triangle_inequality_violation() -> None:
    # B4 check: d_AB=100 violates d_AB <= d_AC + d_BC = 2
    c = _load("C01_sec11_square_nonadditivity_Q.json")
    c["four_point_test"] = {
        "d_AB": [100, 1],
        "d_BC": [1, 1],
        "d_CD": [1, 1],
        "d_DA": [1, 1],
        "d_AC": [1, 1],
        "d_BD": [1, 1],
        "expect_cat0_violation": False,
    }
    with pytest.raises(VerificationError, match="triangle inequality"):
        verify_certificate(c)


def test_mode1_accept_3cut_proportional_additive() -> None:
    c = _load("C03_3cut_cube_nonadditivity_Q.json")
    c["certify_nonadditivity"] = False
    eye2 = [[1, 0], [0, 1]]
    for cut in c["cuts"]:
        cut["matrix_U"] = eye2
        cut["matrix_V"] = eye2
        cut["expected_rank_U"] = 2
        cut["expected_rank_V"] = 2
    assert verify_certificate(c)["status"] == "VERIFIED"


def test_mode3_reject_double_wound_polygon() -> None:
    # B3 exploit: double-wound 3x4 rectangle must be rejected!
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["polygon_vertices"] = c["polygon_vertices"] + c["polygon_vertices"]
    c["edge_lengths"] = c["edge_lengths"] + c["edge_lengths"]
    c["expected_mean_width"] = [14, 1]
    c["interior_point"] = None
    c["submodularity_check"] = None
    with pytest.raises(VerificationError, match="globally strictly CCW convex"):
        verify_certificate(c)


def test_mode3_reject_three_site_cat0() -> None:
    # P2 check: num_sites < 4 must be rejected
    c = _load("C10_cat0_tripod_strict_curvature.json")
    c["num_sites"] = 3
    c["distance_matrix"] = [row[:3] for row in c["distance_matrix"][:3]]
    with pytest.raises(VerificationError, match="< 4"):
        verify_certificate(c)


def test_mode3_reject_degenerate_geodesic_z_eq_m() -> None:
    c = _load("C09_cat0_euclidean_rectangle_equality.json")
    c["geodesic_checks"][0]["z"] = 4
    with pytest.raises(VerificationError, match="degenerate observer"):
        verify_certificate(c)


def test_mode3_reject_cat0_distance_matrix_triangle_violation() -> None:
    c = _load("C10_cat0_tripod_strict_curvature.json")
    c["distance_matrix"][1][2] = [10, 1]
    c["distance_matrix"][2][1] = [10, 1]
    with pytest.raises(VerificationError, match="Triangle inequality failed"):
        verify_certificate(c)


def test_mode3_norm_rigidity_l1_and_bad_delta_interval() -> None:
    # P7: dual_norm = l_1 with f=(1,0), g=(0,1), t=u=1 -> A+=2, A-=2, S=2 -> 8 > 4 passes!
    c = {
        "schema": SCHEMA_VERSION,
        "certificate_id": "NR_L1",
        "mode": "cat0_flat_mean_width",
        "submode": "norm_rigidity_violation",
        "claim_ref": "l_inf plane dual norm l_1 violation",
        "dual_norm": "l_1",
        "t": [1, 1],
        "u": [1, 1],
        "f": [[1, 1], [0, 1]],
        "g": [[0, 1], [1, 1]],
        "delta": [0, 1],
    }
    assert verify_certificate(c)["status"] == "VERIFIED"
    # Non-violating pair f=g=(1,0) at t=u=1: A+=2, A-=0 -> 4+0 = 4 == 2S
    c2 = dict(c, g=[[1, 1], [0, 1]])
    with pytest.raises(VerificationError, match="not violated"):
        verify_certificate(c2)
    # Valid violation at t=2, u=1 (S=5, f=(1,0), g=(0,1): A+=3->9, A-=3->9, 18>10, interval (-4, 4))
    # Picking delta=9/2 lies in (-5, 5) but outside (-4, 4)
    c3 = dict(c, t=[2, 1], u=[1, 1], delta=[9, 2])
    with pytest.raises(VerificationError, match="strict saving interval"):
        verify_certificate(c3)


def test_mode3_polygon_extra_forgeries() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["interior_point"]["convex_weights"] = [[3, 4], [1, 2], [-1, 4], [0, 1]]
    with pytest.raises(VerificationError, match="< 0"):
        verify_certificate(c)
    c2 = _load("C12_mean_width_rectangle_submodularity.json")
    c2["interior_point"]["delta_with_z"] = [8, 1]
    with pytest.raises(VerificationError, match="must preserve delta"):
        verify_certificate(c2)
    c3 = _load("C12_mean_width_rectangle_submodularity.json")
    c3["submodularity_check"]["indices_A"] = [0, 1, 2, 2]
    with pytest.raises(VerificationError, match="Duplicate vertex indices"):
        verify_certificate(c3)
    # Non-square diagonal in sub-polygon: 1x1 square has diagonal sqrt(2) -> _exact_isqrt rejects
    c4 = _load("C12_mean_width_rectangle_submodularity.json")
    c4["polygon_vertices"] = [[[0, 1], [0, 1]], [[1, 1], [0, 1]], [[1, 1], [1, 1]], [[0, 1], [1, 1]]]
    c4["edge_lengths"] = [[1, 1], [1, 1], [1, 1], [1, 1]]
    c4["expected_mean_width"] = [2, 1]
    c4["interior_point"] = None
    c4["submodularity_check"]["delta_A_union_B"] = [2, 1]
    with pytest.raises(VerificationError, match="Non-square integer"):
        verify_certificate(c4)


def test_mode3_prop52_w1_nonzero_and_survives_both() -> None:
    # w1 != 0 with survives_second_degeneration=True must fail (P5)
    c = _load("C13_prop52_quadratic_degeneration.json")
    c["v"] = [1, 1, 2, 3]
    c["g"] = 5  # wait q=5 -> g in 0..4, use h=1, v=[1,2,2,1] -> g=1+5-3=3 != h=1
    c["h"] = 1
    c["v"] = [1, 2, 2, 1]
    c["g"] = 3
    c["survives_first_degeneration"] = False
    c["survives_second_degeneration"] = True
    with pytest.raises(VerificationError, match="cannot be True when w1 != 0"):
        verify_certificate(c)
    c["survives_second_degeneration"] = False
    assert verify_certificate(c)["survives_both"] is False
    # survives_both = True when all v_j = v_1 and g = h
    c_ok = _load("C13_prop52_quadratic_degeneration.json")
    c_ok["v"] = [2, 2, 2, 2]
    c_ok["survives_second_degeneration"] = True
    assert verify_certificate(c_ok)["survives_both"] is True


def test_mode4_reject_positive_scale_nonunit_lambda() -> None:
    # B8 check: Def. 9.1 requires lambda = 1 whenever t_pow_p > 0
    c = _load("C15_routing_budget_and_bregman_probing.json")
    c["links"][0]["t_pow_p"] = [1, 4]
    with pytest.raises(VerificationError, match="positive scale requires lambda = 1"):
        verify_certificate(c)


def test_mode4_reject_negative_diagonal_upper_bound() -> None:
    # B7 check: out-of-regime D > H(I) producing negative upper bound must be rejected
    c = _load("C16_sec16_diagonal_contradiction_gap.json")
    c["external_pairs_E"] = 100
    with pytest.raises(VerificationError, match="out of spectral price regime"):
        verify_certificate(c)


# ---------------------------------------------------------------------------
# 5. Mode 3 Forgeries (CAT(0), Norm Rigidity, Mean Width, Copying Degeneration)
# ---------------------------------------------------------------------------
def test_mode3_reject_asymmetric_distance() -> None:
    c = _load("C09_cat0_euclidean_rectangle_equality.json")
    c["distance_matrix"][0][1] = [4, 1]
    with pytest.raises(VerificationError, match="Invalid distance"):
        verify_certificate(c)


def test_mode3_reject_nonzero_diagonal_distance() -> None:
    c = _load("C09_cat0_euclidean_rectangle_equality.json")
    c["distance_matrix"][0][0] = [1, 1]
    with pytest.raises(VerificationError, match="d\\[0\\]\\[0\\] != 0"):
        verify_certificate(c)


def test_mode3_reject_hamming_in_cat0_mode() -> None:
    c = _load("C10_cat0_tripod_strict_curvature.json")
    c["distance_matrix"] = [
        [[0, 1], [1, 1], [2, 1], [1, 1]],
        [[1, 1], [0, 1], [1, 1], [2, 1]],
        [[2, 1], [1, 1], [0, 1], [1, 1]],
        [[1, 1], [2, 1], [1, 1], [0, 1]],
    ]
    with pytest.raises(VerificationError, match="CAT\\(0\\) Eq. \\(4.5\\) failed"):
        verify_certificate(c)


def test_mode3_reject_wrong_midpoint_label() -> None:
    c = _load("C09_cat0_euclidean_rectangle_equality.json")
    c["geodesic_checks"][0]["m"] = 3
    with pytest.raises(VerificationError, match="is not \\[x,y\\]_s"):
        verify_certificate(c)


def test_mode3_reject_false_strict_curvature_claim() -> None:
    c = _load("C09_cat0_euclidean_rectangle_equality.json")
    c["geodesic_checks"][0]["expect_equality"] = False
    with pytest.raises(VerificationError, match="Expected strict negative curvature"):
        verify_certificate(c)


def test_mode3_reject_false_equality_claim_on_tripod() -> None:
    c = _load("C10_cat0_tripod_strict_curvature.json")
    c["geodesic_checks"][0]["expect_equality"] = True
    with pytest.raises(VerificationError, match="Expected CAT\\(0\\) equality"):
        verify_certificate(c)


def test_mode3_reject_non_unit_dual_vector() -> None:
    c = _load("C11_norm_rigidity_l1_violation.json")
    c["f"] = [[2, 1], [0, 1]]
    with pytest.raises(VerificationError, match="dual unit sphere"):
        verify_certificate(c)


def test_mode3_reject_out_of_range_shift_delta() -> None:
    c = _load("C11_norm_rigidity_l1_violation.json")
    c["delta"] = [2, 1]
    with pytest.raises(VerificationError, match="-S < delta < S"):
        verify_certificate(c)


def test_mode3_reject_non_convex_polygon() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    # Reverse vertices -> CW instead of CCW
    c["polygon_vertices"] = list(reversed(c["polygon_vertices"]))
    with pytest.raises(VerificationError, match="strictly CCW convex"):
        verify_certificate(c)


def test_mode3_reject_spoofed_edge_length() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["edge_lengths"][0] = [4, 1]
    with pytest.raises(VerificationError, match="edge length mismatch"):
        verify_certificate(c)


def test_mode3_reject_spoofed_mean_width() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["expected_mean_width"] = [6, 1]
    with pytest.raises(VerificationError, match="Mean width"):
        verify_certificate(c)


def test_mode3_reject_bad_convex_weights() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["interior_point"]["convex_weights"] = [[1, 2], [1, 4], [1, 4], [1, 4]]
    with pytest.raises(VerificationError, match="sum to 1"):
        verify_certificate(c)


def test_mode3_reject_exterior_point_as_interior() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["interior_point"]["z"] = [[5, 1], [2, 1]]
    with pytest.raises(VerificationError, match="does not equal z"):
        verify_certificate(c)


def test_mode3_reject_disjoint_submodularity_sets() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["submodularity_check"]["indices_A"] = [0, 1]
    c["submodularity_check"]["indices_B"] = [2, 3]
    with pytest.raises(VerificationError, match="A intersect B nonempty"):
        verify_certificate(c)


def test_mode3_reject_spoofed_submod_delta_A() -> None:
    c = _load("C12_mean_width_rectangle_submodularity.json")
    c["submodularity_check"]["delta_A"] = [7, 1]
    with pytest.raises(VerificationError, match="delta_A mismatch"):
        verify_certificate(c)


def test_mode3_reject_prop52_broken_linear_eq() -> None:
    c = _load("C13_prop52_quadratic_degeneration.json")
    c["g"] = 3
    with pytest.raises(VerificationError, match="linear constraint failed"):
        verify_certificate(c)


def test_mode3_reject_prop52_wrong_L() -> None:
    c = _load("C13_prop52_quadratic_degeneration.json")
    c["expected_L"] = 32
    with pytest.raises(VerificationError, match="expected_L"):
        verify_certificate(c)


def test_mode3_reject_prop52_wrong_stage2_flag() -> None:
    c = _load("C13_prop52_quadratic_degeneration.json")
    c["survives_second_degeneration"] = True
    with pytest.raises(VerificationError, match="survives_second_degeneration"):
        verify_certificate(c)


# ---------------------------------------------------------------------------
# 6. Mode 4 Forgeries (Tag(h,w), Catalyst Gap, Routing Budget, Bregman, Eq. 16.1)
# ---------------------------------------------------------------------------
def test_mode4_reject_bad_tag_entry_on_duplicate_addr() -> None:
    c = _load("C14_tag_perm_and_catalyst_substitution.json")
    c["records"][0]["address"] = 0
    c["records"][0]["payload"] = 11
    with pytest.raises(VerificationError, match="Tag\\(h,w\\) entry 0 != expected 1"):
        verify_certificate(c)


def test_mode4_reject_payload_mismatch_in_tag() -> None:
    c = _load("C14_tag_perm_and_catalyst_substitution.json")
    c["records"][0]["payload"] = 13
    with pytest.raises(VerificationError, match="Tag\\(h,w\\) entry 0 != expected 1"):
        verify_certificate(c)


def test_mode4_reject_wrong_frame_price() -> None:
    c = _load("C14_tag_perm_and_catalyst_substitution.json")
    c["expected_frame_price_hX"] = 58
    with pytest.raises(VerificationError, match="h_X=60"):
        verify_certificate(c)


def test_mode4_reject_spoofed_substituted_gap() -> None:
    c = _load("C14_tag_perm_and_catalyst_substitution.json")
    c["expected_substituted_gap"] = [20, 1]
    with pytest.raises(VerificationError, match="substituted gap"):
        verify_certificate(c)


def test_mode4_reject_vacuous_nonpositive_gap() -> None:
    c = _load("C14_tag_perm_and_catalyst_substitution.json")
    c["occurrence_gaps"] = [[10, 1], [12, 1]]
    c["expected_substituted_gap"] = [-2, 1]
    with pytest.raises(VerificationError, match="strictly positive"):
        verify_certificate(c)


def test_mode4_reject_routing_budget_exceeded() -> None:
    c = _load("C15_routing_budget_and_bregman_probing.json")
    c["links"][2]["y"] = 25
    with pytest.raises(VerificationError, match="Eq. \\(9.1\\) budget exceeded"):
        verify_certificate(c)


def test_mode4_reject_non_optimal_kkt_x() -> None:
    c = _load("C15_routing_budget_and_bregman_probing.json")
    c["optimizer_x"][0] = [2, 1]
    with pytest.raises(VerificationError, match="exact KKT argmax"):
        verify_certificate(c)


def test_mode4_reject_nonpositive_bregman_margin() -> None:
    c = _load("C15_routing_budget_and_bregman_probing.json")
    c["probe_fee_P"] = [1, 2]
    with pytest.raises(VerificationError, match="Bregman margin not strictly positive"):
        verify_certificate(c)


def test_mode4_reject_diagonal_gap_failure() -> None:
    c = _load("C16_sec16_diagonal_contradiction_gap.json")
    c["cell_gap_D"] = [1000, 1]
    with pytest.raises(VerificationError, match="Eq. \\(16.1\\) strict contradiction bound failed"):
        verify_certificate(c)


def test_mode4_reject_R_below_8g() -> None:
    c = _load("C16_sec16_diagonal_contradiction_gap.json")
    c["scale_R"] = 16
    with pytest.raises(VerificationError, match="< 32"):
        verify_certificate(c)
