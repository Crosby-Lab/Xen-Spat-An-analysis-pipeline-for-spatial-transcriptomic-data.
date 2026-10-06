"""
Dice overlap between each niche in df1 and the 'Luminal BC Epithelium'
cells in df2, with a spatial (X, Y) plot for every comparison.

Cells are matched by cell ID, not by raw X/Y: the two files store the
coordinates with different precision (e.g. 1554.8275 vs 1554.827533), so
exact coordinate matching finds almost no shared cells.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

# ---------------- settings ----------------
FILE_NICHE = "XeniumdataNiche_100_10_25_91225.csv"   # df1: numeric niches
FILE_LABEL = "100metadataUniform32125_Dice.csv"          # df2: uniform cell types
REFERENCE = "Tumor cells"                  # df2 group to compare against
EXCLUDE_NICHES = []                                 # niches to skip; use [] to include all
MATCH_BY = "id"                                      # "id" (recommended) or "xy"
XY_DECIMALS = 3                                      # only used when MATCH_BY = "xy"
POINT_SIZE = 0.3
OUT_PREFIX = "dice_100"

# image output
SAVE_INDIVIDUAL = True          # True = one image per niche, in addition to the combined figure
SAVE_COMBINED = False            # True = all niches in one panel figure
SAVE_BARPLOT = True             # True = bar chart of Dice values
IMAGE_FORMATS = ["png"]         # any of "png", "pdf", "svg", "tiff"; e.g. ["png", "svg"] for Inkscape
DPI = 300
OUT_DIR = "dice_plots_100"          # folder for the images (created if missing)
# ------------------------------------------


os.makedirs(OUT_DIR, exist_ok=True)


def save(fig, name):
    """Save a figure in every format listed in IMAGE_FORMATS."""
    for ext in IMAGE_FORMATS:
        path = os.path.join(OUT_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        print(f"  saved {path}")


def dice(set1, set2):
    """Dice = 2|A ∩ B| / (|A| + |B|)."""
    total = len(set1) + len(set2)
    return 1.0 if total == 0 else 2.0 * len(set1 & set2) / total


def keys(df, id_col):
    """Return a set of cell keys (IDs or rounded XY pairs)."""
    if MATCH_BY == "id":
        return set(df[id_col])
    return set(zip(df["X"].round(XY_DECIMALS), df["Y"].round(XY_DECIMALS)))


# --- load ---
df1 = pd.read_csv(FILE_NICHE)
df1 = df1[df1["leiden"] != "notassigned"].copy()
df2 = pd.read_csv(FILE_LABEL)

ref = df2[df2["niches"] == REFERENCE]
ref_keys = keys(ref, "cell_id")

niches = sorted(n for n in df1["niches"].unique() if n not in EXCLUDE_NICHES)
results = []

# fixed axis limits so every panel is directly comparable
xlim = (min(df1["X"].min(), df2["X"].min()), max(df1["X"].max(), df2["X"].max()))
ylim = (min(df1["Y"].min(), df2["Y"].min()), max(df1["Y"].max(), df2["Y"].max()))

ncol = min(3, len(niches))
nrow = -(-len(niches) // ncol)
if SAVE_COMBINED:
    fig_all, axes_all = plt.subplots(nrow, ncol, figsize=(6 * ncol, 3.6 * nrow), squeeze=False)


def draw(ax, niche_df, overlap_mask_niche, ref_only, n, d):
    """Background = all cells; blue = niche only; orange = reference only; red = both."""
    ax.scatter(df1["X"], df1["Y"], s=POINT_SIZE, c="#DDDDDD", linewidths=0, rasterized=True)
    ax.scatter(niche_df.loc[~overlap_mask_niche, "X"], niche_df.loc[~overlap_mask_niche, "Y"],
               s=POINT_SIZE, c="#3B75AF", linewidths=0, rasterized=True, label=f"Niche {n} only")
    ax.scatter(ref_only["X"], ref_only["Y"],
               s=POINT_SIZE, c="#EF8636", linewidths=0, rasterized=True, label=f"{REFERENCE} only")
    ax.scatter(niche_df.loc[overlap_mask_niche, "X"], niche_df.loc[overlap_mask_niche, "Y"],
               s=POINT_SIZE, c="#C3423F", linewidths=0, rasterized=True, label="Both (overlap)")
    ax.set_title(f"Niche {n} vs {REFERENCE}\nDice = {d:.3f}", fontsize=11)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.invert_yaxis()                 # image convention: Y increases downwards
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.legend(markerscale=12, fontsize=7, loc="upper right", frameon=True)


for i, n in enumerate(niches):
    niche_df = df1[df1["niches"] == n]
    niche_keys = keys(niche_df, "Cell_ID")
    d = dice(niche_keys, ref_keys)
    inter = niche_keys & ref_keys

    results.append({"niche": n, "reference": REFERENCE,
                    "n_niche": len(niche_keys), "n_reference": len(ref_keys),
                    "n_overlap": len(inter), "dice": round(d, 4)})
    print(f"Niche {n}: Dice = {d:.4f} (overlap {len(inter)} / niche {len(niche_keys)} / ref {len(ref_keys)})")

    # which points are in the overlap / reference-only
    if MATCH_BY == "id":
        in_both = niche_df["Cell_ID"].isin(inter)
        ref_only = ref[~ref["cell_id"].isin(inter)]
    else:
        nk = list(zip(niche_df["X"].round(XY_DECIMALS), niche_df["Y"].round(XY_DECIMALS)))
        rk = list(zip(ref["X"].round(XY_DECIMALS), ref["Y"].round(XY_DECIMALS)))
        in_both = pd.Series([k in inter for k in nk], index=niche_df.index)
        ref_only = ref[[k not in inter for k in rk]]

    # one image per niche
    if SAVE_INDIVIDUAL:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        draw(ax, niche_df, in_both, ref_only, n, d)
        fig.tight_layout()
        save(fig, f"{OUT_PREFIX}_niche{n}")
        plt.close(fig)

    # same plot in the combined figure
    if SAVE_COMBINED:
        draw(axes_all[i // ncol][i % ncol], niche_df, in_both, ref_only, n, d)

if SAVE_COMBINED:
    for j in range(len(niches), nrow * ncol):      # hide unused panels
        axes_all[j // ncol][j % ncol].axis("off")
    fig_all.tight_layout()
    save(fig_all, f"{OUT_PREFIX}_all_niches")
    plt.close(fig_all)

# summary table + bar chart
res = pd.DataFrame(results).sort_values("dice", ascending=False)
best = res.iloc[0]
res["best_match"] = res["niche"] == best["niche"]

# file name carries the best niche and the reference it matched,
# e.g. dice_046_summary_niche6_Luminal_BC_Epithelium.csv
ref_tag = REFERENCE.replace(" ", "_")
summary_file = f"{OUT_PREFIX}_summary_niche{int(best['niche'])}_{ref_tag}.csv"
res.to_csv(summary_file, index=False)
print(res.to_string(index=False))
print(f"\nHighest Dice: niche {int(best['niche'])} = {best['dice']:.4f} -> saved as {summary_file}")
res = res.sort_values("niche")                      # back to niche order for the bar chart

if SAVE_BARPLOT:
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(res["niche"].astype(str), res["dice"], color="#3B75AF")
    for x, v in zip(res["niche"].astype(str), res["dice"]):
        ax.text(x, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_xlabel("Niche")
    ax.set_ylabel("Dice coefficient")
    ax.set_title(f"Overlap with {REFERENCE}")
    ax.set_ylim(0, max(1.0, res["dice"].max() + 0.1))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, f"{OUT_PREFIX}_barplot")
    plt.close(fig)
