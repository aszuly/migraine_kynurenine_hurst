#!/usr/bin/env python
"""Figure 1: schematic of the kynurenine pathway.

Writes figures/Figure1_kynurenine_pathway.pdf and .png (600 dpi).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures", "Figure1_kynurenine_pathway")

BLUE = "#2F5597"      # measured in plasma
ORANGE = "#C55A11"    # L-kynurenine, the primary metabolite of the study
INK = "#222222"
LINE = "#555555"

XL, XC, XR = 1.9, 6.0, 10.1          # column centres
ROW = [10.3, 8.6, 6.9, 5.2, 3.5, 1.8]  # row centres, top to bottom
W, H = 3.1, 1.0

# name: (label, x, y, kind)
BOX = {
    "TRP":   ("Tryptophan\n(TRP)", XC, ROW[0], "measured"),
    "LKYN":  ("L-Kynurenine\n(L-KYN)", XC, ROW[1], "primary"),
    "ANA":   ("Anthranilic acid\n(ANA)", XL, ROW[1], "measured"),
    "KYNA":  ("Kynurenic acid\n(KYNA)", XR, ROW[1], "measured"),
    "HK":    ("3-Hydroxykynurenine\n(3-HK)", XC, ROW[2], "measured"),
    "XA":    ("Xanthurenic acid\n(XA)", XR, ROW[2], "measured"),
    "HANA":  ("3-Hydroxyanthranilic acid\n(3-HANA)", XC, ROW[3], "measured"),
    "ACMS":  ("2-Amino-3-carboxy-\nmuconate semialdehyde\n(ACMS)", XC, ROW[4], "other"),
    "PICA":  ("Picolinic acid\n(PICA)", XR, ROW[4], "measured"),
    "QUINA": ("Quinolinic acid\n(QUINA)", XC, ROW[5], "measured"),
    "NAD":   ("NAD$^{+}$", XR, ROW[5], "other"),
}

STYLE = {
    "measured": dict(fc=BLUE, ec=BLUE, tc="white", lw=1.0),
    "primary":  dict(fc=ORANGE, ec=ORANGE, tc="white", lw=1.0),
    "other":    dict(fc="white", ec=LINE, tc=INK, lw=1.2),
}


def box(ax, label, x, y, kind, w=W, h=H, fs=10.5):
    s = STYLE[kind]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=0.12",
                                fc=s["fc"], ec=s["ec"], lw=s["lw"], zorder=2))
    ax.text(x, y, label, ha="center", va="center", color=s["tc"], fontsize=fs,
            fontweight="bold" if kind != "other" else "normal", linespacing=1.25, zorder=3)


def arrow(ax, p0, p1):
    ax.annotate("", xy=p1, xytext=p0, zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.3,
                                shrinkA=0, shrinkB=0, mutation_scale=13))


def down(ax, a, b, label=None):
    xa, ya = BOX[a][1], BOX[a][2]
    yb = BOX[b][2]
    arrow(ax, (xa, ya - H / 2 - 0.05), (xa, yb + H / 2 + 0.05))
    if label:
        ax.text(xa + 0.18, (ya + yb) / 2, label, ha="left", va="center", fontsize=9.5, color=INK)


def side(ax, a, b, label=None):
    xa, y = BOX[a][1], BOX[a][2]
    xb = BOX[b][1]
    sgn = 1 if xb > xa else -1
    arrow(ax, (xa + sgn * (W / 2 + 0.05), y), (xb - sgn * (W / 2 + 0.05), y))
    if label:
        ax.text((xa + xb) / 2, y + 0.14, label, ha="center", va="bottom", fontsize=9.5, color=INK)


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
fig, ax = plt.subplots(figsize=(7.2, 6.0))
ax.set_xlim(0, 12)
ax.set_ylim(1.05, 11.05)
ax.set_aspect("equal")
ax.axis("off")

for name, (label, x, y, kind) in BOX.items():
    box(ax, label, x, y, kind, fs=9.3 if name == "ACMS" else 10.5)

down(ax, "TRP", "LKYN", "TDO, IDOs")
side(ax, "LKYN", "KYNA", "KATs")
side(ax, "LKYN", "ANA", "KYNU")
down(ax, "LKYN", "HK", "KMO")
side(ax, "HK", "XA", "KATs")
down(ax, "HK", "HANA", "KYNU")
down(ax, "HANA", "ACMS", "3-HAO")
side(ax, "ACMS", "PICA", "ACMSD")
down(ax, "ACMS", "QUINA", "non-enzymatic")
side(ax, "QUINA", "NAD", "QPRT")

# anthranilic acid -> 3-hydroxyanthranilic acid (elbow)
ax.plot([XL, XL], [ROW[1] - H / 2 - 0.05, ROW[3]], color=LINE, lw=1.3, zorder=1, solid_capstyle="projecting")
arrow(ax, (XL, ROW[3]), (XC - W / 2 - 0.05, ROW[3]))

# key
kx, ky, kw, kh = 0.45, ROW[4] + 0.15, 0.55, 0.36
for i, (kind, text) in enumerate([("primary", "primary metabolite of the study"),
                                  ("measured", "measured in plasma"),
                                  ("other", "not measured")]):
    y = ky - i * 0.62
    s = STYLE[kind]
    ax.add_patch(FancyBboxPatch((kx, y - kh / 2), kw, kh, boxstyle="round,pad=0,rounding_size=0.06",
                                fc=s["fc"], ec=s["ec"], lw=s["lw"]))
    ax.text(kx + kw + 0.2, y, text, ha="left", va="center", fontsize=9, color=INK)

fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
fig.savefig(OUT + ".pdf")
fig.savefig(OUT + ".png", dpi=600)
print("written:", OUT + ".pdf/.png")
