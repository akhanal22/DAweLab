import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle  # <-- add this

base_dir = "repeat_blast"

COLORS = {
    "TR1": "#ff3b9d",
    "knob180": "#33cc33",
    "CentC": "#f4c430",
    "Contig": "grey",
}

def read_bed(file):
    cols = ["chr", "start", "end"]
    if (not os.path.exists(file)) or (os.path.getsize(file) == 0):
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(file, sep="\t", header=None, usecols=[0,1,2])
    df.columns = cols
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start","end"])
    return df

tr1_bed     = read_bed(os.path.join(base_dir, "intervals_B73D1_TR1.bed"))
knob180_bed = read_bed(os.path.join(base_dir, "intervals_B73D1_knob180.bed"))
centc_bed   = read_bed(os.path.join(base_dir, "intervals_B73D1_CentC.bed"))  # fixed

length_df = pd.read_csv(os.path.join(base_dir, "length.txt"), sep="\t", header=None)
total_length_mb = int(length_df.iloc[0, 1]) / 1e6

agp = pd.read_csv(os.path.join(base_dir, "chr7.agp"), sep="\t", header=None, usecols=[0,1,2,5])
agp.columns = ["chr","start","end","ID"]
agp["start"] = pd.to_numeric(agp["start"])
agp["end"]   = pd.to_numeric(agp["end"])

fig, ax = plt.subplots(figsize=(12,3))

bar_height = 0.4
y_center = 0.0
y_bottom = y_center - bar_height/2

contig_h = bar_height * 0.35
repeat_h = 0.38
gap = 0.05
contig_y = y_bottom - contig_h - gap
repeat_y = y_center - repeat_h/2

ax.add_patch(Rectangle((0,y_bottom), total_length_mb, bar_height, facecolor="#1673A3", edgecolor="none", zorder=1))

for _, row in agp.iterrows():
    x0 = row["start"]/1e6
    w  = (row["end"]-row["start"])/1e6
    ax.add_patch(Rectangle((x0,contig_y), w, contig_h, facecolor=COLORS["Contig"], edgecolor="black", linewidth=0.5, zorder=2))

def add_track(df, color, z=10):
    if df.empty: return
    for _, r in df.iterrows():
        x0 = r["start"]/1e6
        w  = (r["end"]-r["start"])/1e6
        if w > 0:
            ax.add_patch(Rectangle((x0,repeat_y), w, repeat_h, facecolor=color, edgecolor="none", alpha=0.95, zorder=z))

add_track(tr1_bed, COLORS["TR1"], z=10)
add_track(knob180_bed, COLORS["knob180"], z=10)
add_track(centc_bed, COLORS["CentC"], z=10)

ax.set_xlim(0, total_length_mb*1.02)
ax.set_ylim(contig_y-0.15, y_bottom+bar_height+0.15)
ax.set_xlabel("Position (Mb)")
ax.set_yticks([])
ax.set_xticks(np.arange(0, total_length_mb+10, 10))

output_pdf = os.path.join(base_dir, "Chr7_agp_landscape.pdf")
fig.tight_layout()
fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
print("Saved to:", output_pdf)
plt.show()
