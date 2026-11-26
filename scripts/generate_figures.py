"""Generate 3 publication-quality figures for Paper 68 (TENSOR-SPEC-CERT-1.0)."""

from __future__ import annotations

from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 300,
})


def make_fig1() -> None:
    """Fig. 1: Sec. 1.1 Square Non-Additivity & Eq. (4.5) CAT(0) Contrast."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.85))
    fig.subplots_adjust(wspace=0.34, left=0.08, right=0.98, bottom=0.19, top=0.84)

    # Panel (a): Geometric mean G(U)+G(V) = 2*sqrt(2) vs G(U+V) = 3 (C01-C02)
    labels = ["$G(U)$\n$(2,1)$", "$G(V)$\n$(1,2)$", "$G(U)+G(V)$\n(additive)", "$G(U\\oplus V)$\n$(3,3)$"]
    vals = [np.sqrt(2), np.sqrt(2), 2 * np.sqrt(2), 3.0]
    colors = ["#4C78A8", "#72B7B2", "#F58518", "#E45756"]
    bars = ax1.bar(range(len(labels)), vals, color=colors, edgecolor="black", linewidth=0.7, width=0.58)
    ax1.axhline(2 * np.sqrt(2), color="#F58518", linestyle="--", linewidth=0.9, alpha=0.8)
    for idx, (bar, v) in enumerate(zip(bars, vals)):
        txt = f"{v:.3f}" if idx != 3 else "3.000"
        ax1.text(bar.get_x() + bar.get_width() / 2, v + 0.08, txt, ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax1.annotate(
        "Strict gap: $9 > 8$",
        xy=(2.95, 3.02),
        xytext=(0.15, 3.18),
        arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "#E45756", "connectionstyle": "arc3,rad=-0.12"},
        fontsize=7.5,
        color="#B22222",
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.24", "fc": "#FFF5F5", "ec": "#E45756", "lw": 0.7},
    )
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, fontsize=7.5)
    ax1.set_ylabel("Midline geometric mean $G(\\cdot)$")
    ax1.set_ylim(0, 3.78)
    ax1.set_title("(a) Sec. 1.1 Direct-Sum Obstruction (C01)")
    ax1.grid(axis="y", linestyle=":", alpha=0.5)

    # Panel (b): Eq. (4.5) 4-point CAT(0) ratio (d_AC^2 + d_BD^2) / sum(side^2)
    configs = ["Tripod\n(C10)", "3x4 Rect.\n(C09)", "Hamming\nSq. (C01)"]
    ratios = [0.50, 1.00, 2.00]
    c2 = ["#54A24B", "#4C78A8", "#E45756"]
    b2 = ax2.bar(range(len(configs)), ratios, color=c2, edgecolor="black", linewidth=0.7, width=0.52)
    ax2.axhline(1.0, color="#222222", linestyle="--", linewidth=1.0, label="CAT(0) ceiling $= 1$")
    for bar, r in zip(b2, ratios):
        ax2.text(bar.get_x() + bar.get_width() / 2, r + 0.06, f"{r:.2f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold")
    ax2.set_xticks(range(len(configs)))
    ax2.set_xticklabels(configs, fontsize=7.8)
    ax2.set_ylabel("$(d_{AC}^2+d_{BD}^2)/\\sum d_{\\mathrm{side}}^2$")
    ax2.set_ylim(0, 2.42)
    ax2.set_title("(b) Eq. (4.5) Four-Point Curvature Test")
    ax2.legend(loc="upper left", frameon=True)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)

    fig.savefig(FIG_DIR / "fig1_nonadditivity_and_cat0.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig1_nonadditivity_and_cat0.png", bbox_inches="tight")
    plt.close(fig)


def make_fig2() -> None:
    """Fig. 2: Exterior Crossing Tradeoff (Mode 2) & Planar Mean-Width / Submodularity (Mode 3)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    fig.subplots_adjust(wspace=0.34, left=0.08, right=0.98, bottom=0.18, top=0.84)

    # Panel (a): Honest exact discrete exterior-algebra crossing ranks (Sec. 3.2 Eq. 3.2-3.3)
    # Compare trivial single-cut rank d_X*d_Y vs certified crossing ranks (r_+, r_-) for C05, C08, C06!
    inst = ["C05\n$(Q,b)=(2,1)$", "C08\n$(3,1)$", "C06\n$(4,1)$"]
    x_pos = np.arange(len(inst))
    w_bar = 0.25
    triv_rk = [4, 9, 16]
    r_plus = [1, 3, 6]
    r_minus = [1, 1, 1]
    b1 = ax1.bar(x_pos - w_bar, triv_rk, width=w_bar, color="#BAB0AC", edgecolor="black", linewidth=0.7, label="Trivial $d_X d_Y$")
    b2 = ax1.bar(x_pos, r_plus, width=w_bar, color="#4C78A8", edgecolor="black", linewidth=0.7, label="Cut $X^+Y^+\\mid X^-Y^-$ ($r_+$)")
    b3 = ax1.bar(x_pos + w_bar, r_minus, width=w_bar, color="#E45756", edgecolor="black", linewidth=0.7, label="Cut $X^+Y^-\\mid X^-Y^+$ ($r_-$)")
    for bars in (b1, b2, b3):
        for bar in bars:
            h = int(bar.get_height())
            ax1.text(bar.get_x() + bar.get_width() / 2, h + 0.35, str(h), ha="center", va="bottom", fontsize=7.4, fontweight="bold")
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(inst, fontsize=7.6)
    ax1.set_ylabel("Exact $\\mathbb{Z}$-flattening rank")
    ax1.set_ylim(0, 19.5)
    ax1.set_title("(a) Sec. 3.2 Exterior Crossing Ranks")
    ax1.legend(loc="upper left", fontsize=6.5, frameon=True)
    ax1.grid(axis="y", linestyle=":", alpha=0.5)

    # Panel (b): C12 3x4 Rectangle mean width W=7, interior copy z=(3/2,2), triangles A & B
    poly = np.array([[0, 0], [3, 0], [3, 4], [0, 4], [0, 0]])
    ax2.fill([0, 3, 3], [0, 0, 4], color="#4C78A8", alpha=0.16, label="$A$: $W(A)=6$")
    ax2.fill([0, 3, 0], [0, 4, 4], color="#F58518", alpha=0.16, label="$B$: $W(B)=6$")
    ax2.plot(poly[:, 0], poly[:, 1], "k-", linewidth=1.4)
    ax2.plot([0, 3], [0, 4], color="#54A24B", linestyle="--", linewidth=1.6, label="$A\\cap B$: $d=5$")
    ax2.scatter([0, 3, 3, 0], [0, 0, 4, 4], color="#222222", s=28, zorder=5)
    ax2.scatter([1.5], [2.0], color="#E45756", edgecolor="black", s=55, marker="*", zorder=6, label="$z=(\\frac{3}{2},2)$ ($\\Delta=0$)")
    ax2.text(1.5, -0.45, "$d=3$", ha="center", fontsize=7.8)
    ax2.text(3.22, 2.0, "$d=4$", va="center", fontsize=7.8)
    ax2.text(1.62, 1.62, "$z$", color="#B22222", fontsize=8.5, fontweight="bold")
    ax2.text(0.15, 3.45, "$W(K)=\\frac{14}{2}=7$", fontsize=7.8, fontweight="bold",
             bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "#666666", "lw": 0.6})
    ax2.set_xlim(-0.55, 3.95)
    ax2.set_ylim(-0.75, 4.65)
    ax2.set_xlabel("$x_1$")
    ax2.set_ylabel("$x_2$")
    ax2.set_title("(b) Thm. 5.3 Mean Width & Submodularity (C12)")
    ax2.legend(loc="lower right", fontsize=6.6, frameon=True)
    ax2.grid(linestyle=":", alpha=0.4)

    fig.savefig(FIG_DIR / "fig2_crossing_and_mean_width.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig2_crossing_and_mean_width.png", bbox_inches="tight")
    plt.close(fig)


def make_fig3() -> None:
    """Fig. 3: Parts II-III Catalyst Substitution, Routing Budget & Eq. (16.1) Diagonal Gap."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8))
    fig.subplots_adjust(wspace=0.35, left=0.09, right=0.98, bottom=0.19, top=0.84)

    # Panel (a): C14 Lemma 8.2 catalyst frame price vs gap budget & C15 routing budget utilization
    cats = ["Frame $h_X$\n(C14)", "Base $D$\n(C14)", "Subst. $D_{\\mathrm{sub}}$\n(C14)", "Budget\n(C15)", "Charge\n(C15)"]
    vals = [60, 36, 18, 120, 112.5]
    cols = ["#4C78A8", "#72B7B2", "#54A24B", "#BAB0AC", "#F58518"]
    bars = ax1.bar(range(len(cats)), vals, color=cols, edgecolor="black", linewidth=0.7, width=0.58)
    for b, v in zip(bars, vals):
        lbl_v = f"{v:.1f}" if v != int(v) else str(int(v))
        ax1.text(b.get_x() + b.get_width() / 2, v + 2.5, lbl_v, ha="center", va="bottom", fontsize=7.6, fontweight="bold")
    ax1.set_xticks(range(len(cats)))
    ax1.set_xticklabels(cats, fontsize=7.2)
    ax1.set_ylabel("Exact rational units")
    ax1.set_ylim(0, 138)
    ax1.set_title("(a) Catalyst Gaps (C14) & Routing Budget (C15)")
    ax1.grid(axis="y", linestyle=":", alpha=0.5)

    # Panel (b): C16 Eq. (16.1) in-regime diagonal contradiction breakdown (units of 10^3)
    steps = ["$H(I)=2RE$", "$-D$ saving", "$+\\mathrm{err}_{\\mathrm{grid}}$", "$+\\mathrm{err}_{\\mathrm{sel}}$", "Upper bnd"]
    values_k = [128.000, -50.000, 0.500, 16.896, 95.396]
    col_w = ["#4C78A8", "#54A24B", "#F58518", "#F58518", "#E45756"]
    b2 = ax2.bar(range(len(steps)), values_k, color=col_w, edgecolor="black", linewidth=0.7, width=0.56)
    ax2.axhline(128.000, color="#4C78A8", linestyle="--", linewidth=0.95, label="$H(I)=128{,}000$")
    ax2.axhline(0.0, color="black", linewidth=0.7)
    for bar, v in zip(b2, values_k):
        if v >= 0:
            ax2.text(bar.get_x() + bar.get_width() / 2, v + 3.2, f"{v:+.2f}k", ha="center", va="bottom", fontsize=6.9, fontweight="bold")
        else:
            ax2.text(bar.get_x() + bar.get_width() / 2, v - 3.5, f"{v:+.2f}k", ha="center", va="top", fontsize=6.9, fontweight="bold")
    ax2.set_xticks(range(len(steps)))
    ax2.set_xticklabels(steps, fontsize=7.0)
    ax2.set_ylabel("Price contribution ($\\times 10^3$)")
    ax2.set_ylim(-68, 152)
    ax2.set_title("(b) Eq. (16.1) In-Regime Gap Contradiction (C16)")
    ax2.legend(loc="lower right", frameon=True)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)

    fig.savefig(FIG_DIR / "fig3_catalyst_and_diagonal.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig3_catalyst_and_diagonal.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    make_fig1()
    make_fig2()
    make_fig3()
    print("Generated 3 figures in paper/figures/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
