#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ==================================================
# USER SETTINGS
# ==================================================
# Run this script from the directory that contains these folders:
#   ragtag_B73D1_v3/, ragtag_B73D2_v3/, ...
SAMPLES = [
    "ragtag_B73D1_v3",
    "ragtag_B73D2_v3",
    "ragtag_B73D4_v3",
    "ragtag_B73D7B_v3",
    "ragtag_B73D7c_v3",
    "ragtag_B73D10_v3",
    "ragtag_B73all_v3",
    "ragtag_J702_v3",
]

SUBDIR = "repeat_blast"
CHR = "chr7"

# Draw AGP contig blocks below the backbone if an AGP exists
DRAW_AGP_CONTIGS = True

# ==================================================
# COLORS
# ==================================================
COLORS = {
    "TR1": "#ff3b9d",
    "knob180": "#33cc33",
    "CentC": "#f4c430",
    "Contig": "lightgray",
}
BACKBONE_COLOR = "#1673A3"

# ==================================================
# HELPERS
# ==================================================
def hap_name_from_dir(sample_dir: str) -> str:
    """Label shown on plot: ragtag_B73D1_v3 -> B73D1_v3"""
    base = os.path.basename(os.path.normpath(sample_dir))
    return base.replace("ragtag_", "")

def file_prefix_from_dir(sample_dir: str) -> str:
    """Prefix used in filenames: ragtag_B73D1_v3 -> B73D1"""
    base = os.path.basename(os.path.normpath(sample_dir)).replace("ragtag_", "")
    return base.split("_")[0]

def first_glob(folder: str, patterns):
    """Return first match for any of patterns, else None."""
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(folder, pat)))
        if hits:
            return hits[0]
    return None

def norm_chr(x: str) -> str:
    """
    Normalize chromosome strings to 'chrN'
    Handles: chr07, Chr7, 7, chr7_RagTag, chr7-something, etc.
    """
    s = str(x).strip()
    s = s.replace("Chr", "chr").replace("CHR", "chr")

    # If it starts with chr, keep only chr + leading digits
    if s.lower().startswith("chr"):
        tail = s[3:]
        digits = ""
        for ch in tail:
            if ch.isdigit():
                digits += ch
            else:
                break
        if digits:
            return f"chr{int(digits)}"  # chr07_RagTag -> chr7
        return s

    # If it starts with digits, use them
    digits = ""
    for ch in s:
        if ch.isdigit():
            digits += ch
        else:
            break
    if digits:
        return f"chr{int(digits)}"

    return s

def read_chr_length_mb(folder: str, chrom: str) -> float | None:
    """
    Use <chrom>.txt (e.g., chr7.txt) to infer chromosome length: max(position) / 1e6.
    Assumes <chrom>.txt has at least 2 columns where col2 is position.
    """
    f = os.path.join(folder, f"{chrom}.txt")
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return None

    df = pd.read_csv(f, sep=r"\s+|\t", engine="python", header=None)
    if df.shape[1] < 2:
        return None

    pos = pd.to_numeric(df.iloc[:, 1], errors="coerce").dropna()
    if pos.empty:
        return None

    return float(pos.max()) / 1e6

def read_bed(path: str) -> pd.DataFrame:
    """Read a BED-like file with at least 3 columns: chr start end (extra columns OK)."""
    cols = ["chr", "start", "end"]
    if (not path) or (not os.path.exists(path)) or (os.path.getsize(path) == 0):
        return pd.DataFrame(columns=cols)

    df = pd.read_csv(path, sep=r"\s+|\t", engine="python", header=None)
    if df.shape[1] < 3:
        return pd.DataFrame(columns=cols)

    df = df.iloc[:, :3].copy()
    df.columns = cols

    df["chr"] = df["chr"].map(norm_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    return df

def read_agp(folder: str, chrom: str) -> pd.DataFrame:
    """
    Read an AGP file for chromosome contig blocks. Tries common names.
    AGP columns: object, object_beg, object_end, part_number, component_type, component_id/gap_length, ...
    """
    agp_path = first_glob(folder, [f"{chrom}.agp", f"*{chrom}*.agp", "*.agp"])
    if not agp_path:
        return pd.DataFrame(columns=["chr", "start", "end", "ctype"])

    agp = pd.read_csv(agp_path, sep="\t", header=None, comment="#")
    if agp.shape[1] < 6:
        return pd.DataFrame(columns=["chr", "start", "end", "ctype"])

    agp = agp.iloc[:, [0, 1, 2, 4]].copy()
    agp.columns = ["chr", "start", "end", "ctype"]
    agp["chr"] = agp["chr"].map(norm_chr)
    agp["start"] = pd.to_numeric(agp["start"], errors="coerce")
    agp["end"]   = pd.to_numeric(agp["end"], errors="coerce")
    agp = agp.dropna(subset=["start", "end"])
    return agp

def read_repeat_tracks(folder: str, file_prefix: str):
    """
    Your naming convention:
      intervals_B73D1_TR1.bed
      intervals_B73D1_knob180.bed
      intervals_B73D1_CentC.bed
    """
    tr1_path   = os.path.join(folder, f"intervals_{file_prefix}_TR1.bed")
    knob_path  = os.path.join(folder, f"intervals_{file_prefix}_knob180.bed")
    centc_path = os.path.join(folder, f"intervals_{file_prefix}_CentC.bed")

    tr1 = read_bed(tr1_path)
    knob = read_bed(knob_path)
    centc = read_bed(centc_path)

    # filter to CHR
    tr1 = tr1[tr1["chr"] == CHR] if not tr1.empty else tr1
    knob = knob[knob["chr"] == CHR] if not knob.empty else knob
    centc = centc[centc["chr"] == CHR] if not centc.empty else centc

    return tr1, knob, centc

def add_track(ax, df, color, y, h, z=10, alpha=0.95):
    """Add rectangles for intervals in df (bp -> Mb)."""
    if df.empty:
        return
    for _, r in df.iterrows():
        s = float(r["start"])
        e = float(r["end"])
        if e <= s:
            continue
        x0 = s / 1e6
        w  = (e - s) / 1e6
        ax.add_patch(Rectangle((x0, y), w, h,
                               facecolor=color, edgecolor="none",
                               alpha=alpha, zorder=z))

# ==================================================
# LOAD ALL SAMPLES
# ==================================================
sample_data = []
for sdir in SAMPLES:
    folder = os.path.join(sdir, SUBDIR)
    if not os.path.isdir(folder):
        print(f"WARNING: missing folder: {folder}  (skipping)")
        continue

    label = hap_name_from_dir(sdir)         # e.g., B73D1_v3 (plot label)
    prefix = file_prefix_from_dir(sdir)     # e.g., B73D1    (filename prefix)

    total_mb = read_chr_length_mb(folder, CHR)
    if total_mb is None:
        print(f"WARNING: missing/invalid {CHR}.txt in: {folder}  (skipping)")
        continue

    tr1_df, knob_df, centc_df = read_repeat_tracks(folder, prefix)

    agp_df = read_agp(folder, CHR) if DRAW_AGP_CONTIGS else pd.DataFrame()
    if not agp_df.empty:
        agp_df = agp_df[agp_df["chr"] == CHR]

    print(f"{label}: prefix={prefix}  TR1={len(tr1_df)} knob180={len(knob_df)} CentC={len(centc_df)}")

    sample_data.append({
        "label": label,
        "total_mb": total_mb,
        "tr1": tr1_df,
        "knob": knob_df,
        "centc": centc_df,
        "agp": agp_df,
    })

if not sample_data:
    raise RuntimeError("No samples loaded. Check SAMPLES paths and repeat_blast contents.")

max_len = max(d["total_mb"] for d in sample_data)

# ==================================================
# PLOT (STACKED)
# ==================================================
n = len(sample_data)
fig_h = max(3.5, 1.15 * n)
fig, ax = plt.subplots(figsize=(14, fig_h))

row_gap = 1.0
backbone_h = 0.40
repeat_h = 0.38
contig_h = backbone_h * 0.35
pad = 0.05

for i, d in enumerate(sample_data):
    # plot top -> bottom in the same order as SAMPLES
    y_center = (n - 1 - i) * row_gap
    y_backbone = y_center - backbone_h / 2
    y_repeat = y_center - repeat_h / 2
    y_contig = y_backbone - contig_h - pad

    # Backbone
    ax.add_patch(Rectangle(
        (0, y_backbone),
        d["total_mb"],
        backbone_h,
        facecolor=BACKBONE_COLOR,
        edgecolor="none",
        zorder=1
    ))

    # AGP contigs (below backbone)
    if DRAW_AGP_CONTIGS and (d["agp"] is not None) and (not d["agp"].empty):
        agp_comp = d["agp"]
        # If you want only contigs (not gaps), uncomment:
        # agp_comp = agp_comp[agp_comp["ctype"] == "W"]
        for _, r in agp_comp.iterrows():
            s = float(r["start"]) / 1e6
            e = float(r["end"]) / 1e6
            if e <= s:
                continue
            ax.add_patch(Rectangle(
                (s, y_contig),
                e - s,
                contig_h,
                facecolor=COLORS["Contig"],
                edgecolor="black",
                linewidth=0.5,
                zorder=2
            ))

    # Repeat tracks
    add_track(ax, d["tr1"],   COLORS["TR1"],     y_repeat, repeat_h, z=10)
    add_track(ax, d["knob"],  COLORS["knob180"], y_repeat, repeat_h, z=10)
    add_track(ax, d["centc"], COLORS["CentC"],   y_repeat, repeat_h, z=10)

    # Label
    ax.text(-0.01 * max_len, y_center, d["label"], va="center", ha="right", fontsize=10)

# Axes formatting
ax.set_xlim(0, max_len * 1.02)
ax.set_ylim(-0.8, (n - 1) * row_gap + 0.8)
ax.set_yticks([])
ax.set_xlabel(f"Position on {CHR} (Mb)")
ax.set_title(f"{CHR} repeat landscape (stacked)")

# ticks every 10 Mb
ax.set_xticks(np.arange(0, max_len + 10, 10))

# Legend
legend_handles = [
    Rectangle((0, 0), 1, 1, facecolor=COLORS["TR1"], edgecolor="none", label="TR1"),
    Rectangle((0, 0), 1, 1, facecolor=COLORS["knob180"], edgecolor="none", label="knob180"),
    Rectangle((0, 0), 1, 1, facecolor=COLORS["CentC"], edgecolor="none", label="CentC"),
]
if DRAW_AGP_CONTIGS:
    legend_handles.append(Rectangle((0, 0), 1, 1, facecolor=COLORS["Contig"], edgecolor="black", label="AGP contig"))

ax.legend(handles=legend_handles, loc="upper right", frameon=False)

fig.tight_layout()
outfile = f"{CHR}_stacked_landscape.png"
fig.savefig(outfile, dpi=300, bbox_inches="tight")
print("Saved:", outfile)

# On HPC this may warn; safe to keep or remove
# plt.show()
