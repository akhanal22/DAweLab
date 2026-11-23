#!/usr/bin/env python3
"""J702 region enrichment / density with chr7 + chr8 knob180 blocks broken at
bad bins, bad bins excluded everywhere, and event length bias removed.

Changes relative to 07-region_enrichment_three_shuffles.py
  * Bad bins: union of the five lines' SyRI bad 500 kb bins (the bins used to
    filter the SyRI VCF that the validated SVs come from; 31.0 Mb). They are
    used (i) as knob180 breakpoints and (ii) as bedtools -excl for every
    shuffle (with the AGP U gaps). No events are removed by bad bins: the
    validated SVs are already outside their own line's bins.
  * Knob180: raw EDTA 'Name=knob180' hits on chr7 and chr8 only. A hit that
    overlaps a bad bin or AGP gap is dropped; the remaining hits are merged when
    <= --knob-merge-kb apart, but never across a bad bin / gap (these act as
    breakpoints). The resulting blocks are the knob180 region.
  * Centromere: tentative J702 spans (AGP gaps removed).
  * Length bias: every event (SNV, indel, SV) is represented by its midpoint
    (1 bp), so a 5 kb deletion is no more likely to land in a region than an SNV.
    Regions are an exclusive partition on the midpoint, priority
    Centromere > Knob180 > Genic > Intergenic (genic not plotted).
  * Shuffles (bedtools shuffle -chrom -noOverlapping -excl gaps+bad bins -seed i):
      knob180 blocks, centromere pieces, intergenic intervals (largest first,
      -maxTries 10000). Shuffled knob180 counts exclude real-centromere
      midpoints (same rule as real); region bp is measured per shuffle.

Outputs (prefix j702_region_badbin_knob):
  .fold.{png,pdf,svg}     fold vs SNVs, six bars per class, no y = 1 line
  .density.{png,pdf,svg}  three panels: SNVs/Mb, indels/Mb, SVs/Mb
  .tsv, .per_line.tsv, .null.tsv, regions in --work-dir

Statistics
  fold  = (class in R / class total) / (SNV in R / SNV total)
  density = events per Mb of region per MA line (pooled count / 5 / Mb)
  Real vs shuffled: two-sided empirical p around the null mean,
      p = (1 + #{|null_k - m| >= |obs - m|}) / (N + 1)
  Real fold vs SNV: Fisher exact (class vs SNV, in vs out), BH over 6 -> '*'.
"""
import argparse
import glob
import importlib.util
import multiprocessing as mp
import os
import subprocess
from bisect import bisect_left
from collections import defaultdict
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('s06', HERE / '06-region_enrichment_raw_knob180_shuffle.py')
s06 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s06)
reg, perm, LENS, merge, SAMPLES, bh = s06.reg, s06.perm, s06.LENS, s06.merge, s06.SAMPLES, s06.bh
CENTROMERES, raw_knob, as_arrays, hit = s06.CENTROMERES, s06.raw_knob, s06.as_arrays, s06.hit

KNOB_CHROMS = ['chr7_RagTag', 'chr8_RagTag']
REAL = ['Knob180', 'Centromere', 'Intergenic']
SHUF = {r: f'Shuffled {r.lower()}' for r in REAL}
REGIONS = [x for r in REAL for x in (r, SHUF[r])]
CLASSES = ['Indels', 'Validated SVs']
ALL = ['SNVs'] + CLASSES
COLOR = {'Knob180': '#8c510a', 'Shuffled knob180': '#dfc27d',
         'Centromere': '#444444', 'Shuffled centromere': '#b0b0b0',
         'Intergenic': '#5f7a95', 'Shuffled intergenic': '#c3cfdb'}
NULL_LABEL = lambda r, N: f'Shuffled {r} ({N:,}×)'   # legend text for null bars
EXTRA = {'Knob180': ['-noOverlapping'], 'Centromere': ['-noOverlapping'],
         'Intergenic': ['-noOverlapping', '-maxTries', '10000']}


# ---------------- interval helpers ----------------
def subtract(A, B):
    """merged A minus merged B."""
    out, j = [], 0
    for a, b in A:
        cur = a
        while j < len(B) and B[j][1] <= cur:
            j += 1
        k = j
        while k < len(B) and B[k][0] < b:
            if B[k][0] > cur:
                out.append((cur, B[k][0]))
            cur = max(cur, B[k][1])
            k += 1
        if cur < b:
            out.append((cur, b))
    return out


def bp(iv):
    return sum(b - a for a, b in iv)


def bad_bins(pattern):
    d = defaultdict(list)
    for f in sorted(glob.glob(pattern)):
        for line in open(f):
            x = line.split()
            if len(x) >= 3 and x[0] in LENS:
                d[x[0]].append((int(x[1]), int(x[2])))
    return {c: merge(d[c]) for c in LENS}


def knob_blocks(raw, brk, gap):
    """hits on one chromosome -> blocks merged <= gap apart, broken at brk."""
    starts = [x for x, _ in brk]
    blocks = []
    for a, b in sorted(raw):
        j = bisect_left(starts, b) - 1          # last breakpoint starting before b
        if j >= 0 and brk[j][1] > a:
            continue                             # hit overlaps bad bin / gap: drop
        seg = bisect_left(starts, a)             # breakpoint segment index
        if blocks and a - blocks[-1][1] <= gap and blocks[-1][2] == seg:
            blocks[-1][1] = max(blocks[-1][1], b)
        else:
            blocks.append([a, b, seg])
    return [(a, b) for a, b, _ in blocks]


def ci95(v):
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return np.nan, np.nan
    m = v.mean(); h = t.ppf(.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
    return max(m - h, 0), m + h


def emp_p(obs, null):
    null = np.asarray(null, float); null = null[np.isfinite(null)]
    m = null.mean()
    return (1 + np.sum(np.abs(null - m) >= abs(obs - m))) / (len(null) + 1)


# ---------------- shuffle worker ----------------
_W = {}


def one_shuffle(i):
    w = _W
    row = {'seed': i}
    for r in REAL:
        p = subprocess.run([w['bedtools'], 'shuffle', '-i', str(w['beds'][r]), '-g', str(w['wd'] / 'genome.txt'),
                            '-excl', str(w['wd'] / 'excluded_gaps_badbins.bed'), '-chrom', '-seed', str(i)] + EXTRA[r],
                           check=True, capture_output=True, text=True)
        if i <= 3:
            (w['sdir'] / f'{r.lower()}_shuffle_seed{i}.bed').write_text(p.stdout)
        df = pd.read_csv(StringIO(p.stdout), sep='\t', header=None, usecols=[0, 1, 2]) if p.stdout else \
            pd.DataFrame(columns=[0, 1, 2])
        row[f'{r}|placed'] = len(df)
        tot_bp = 0
        for cls in ALL:
            for l in SAMPLES:
                row[f'{r}|{cls}|{l}'] = 0
        for c, g in df.groupby(0):
            iv = merge(zip(g[1].to_numpy(), g[2].to_numpy()))
            if r == 'Knob180':                    # same rule as real: not in a real centromere
                iv = subtract(iv, w['cen'][c])
            tot_bp += bp(iv)
            arr = as_arrays(iv)
            for cls in ALL:
                m, smp = w['mid'][cls][c]
                h = hit(m, m + 1, arr)
                for l in SAMPLES:
                    row[f'{r}|{cls}|{l}'] += int((h & (smp == l)).sum())
        row[f'{r}|bp'] = tot_bp
    return row


# ---------------- plotting ----------------
def setup_rc():
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': 8, 'axes.linewidth': 0.7, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})


def draw_group(ax, x0, res, idf, cls, metric, N, label_fmt, star_col=None):
    """six bars (real/shuffled pairs) centred at x0; returns bar x and tops."""
    step = width = 0.12
    xs, tops, txt = {}, {}, {}
    for j, r in enumerate(REGIONS):
        x = x0 + (j - (len(REGIONS) - 1) / 2) * step
        row = res[(res.category == cls) & (res.region == r)].iloc[0]
        y = float(row[metric])
        dots = idf[(idf.category == cls) & (idf.region == r)].set_index('MA_line').loc[SAMPLES][metric].to_numpy(float)
        ax.bar(x, y, width=width, color=COLOR[r], edgecolor='none', zorder=2)
        top = y
        lo, hi = ci95(dots)
        if np.isfinite(hi):
            ax.plot([x, x], [lo, hi], color='#202020', lw=0.8, zorder=4)
            for yy in (lo, hi):
                ax.plot([x - .03, x + .03], [yy, yy], color='#202020', lw=0.8, zorder=4)
            top = max(top, hi)
        ok = np.isfinite(dots)
        if ok.any():
            ax.scatter(x + np.linspace(-.04, .04, 5)[ok], dots[ok], s=6, color='#1a1a1a',
                       edgecolors='white', linewidths=0.2, zorder=5)
            top = max(top, dots[ok].max())
        star = '*' if (star_col and r in REAL and row[star_col] < 0.05) else ''
        xs[r], tops[r] = x, top
        txt[r] = ax.text(x, top, label_fmt(y) + star, ha='center', va='bottom', rotation=90, fontsize=6.5, zorder=10)
    return xs, tops, txt


def finish_brackets(fig, ax, groups, res, pcol, N, ymax):
    """place labels/brackets above everything; groups = list of (cls, xs, tops, txt)."""
    ylim = ymax * 1.8
    for _ in range(8):
        ax.set_ylim(0, ylim)
        pad = ylim * 0.012
        for _, xs, tops, txt in groups:
            for r, tx in txt.items():
                tx.set_y(tops[r] + pad)
        fig.canvas.draw()
        inv = ax.transData.inverted()
        need = 0
        for _, xs, tops, txt in groups:
            lt = {r: inv.transform((0, tx.get_window_extent().y1))[1] for r, tx in txt.items()}
            need = max(need, max(max(lt[r], lt[SHUF[r]]) for r in REAL))
        need = (need + ylim * 0.09) / 0.82
        if abs(need - ylim) / ylim < 0.02:
            break
        ylim = need
    inv = ax.transData.inverted()
    for cls, xs, tops, txt in groups:
        lt = {r: inv.transform((0, tx.get_window_extent().y1))[1] for r, tx in txt.items()}
        for r in REAL:
            p = res[(res.category == cls) & (res.region == SHUF[r])].iloc[0][pcol]
            thr = 1 / (N + 1)
            ptxt = (f'p < {thr:.3f}' if thr >= 5e-4 else f'p < {thr:.0e}') if p <= thr + 1e-12 \
                else (f'p = {p:.3f}' if p >= 1e-3 else f'p = {p:.1e}')
            yb = max(lt[r], lt[SHUF[r]]) + ylim * 0.025
            x1, x2 = xs[r], xs[SHUF[r]]
            ax.plot([x1, x1, x2, x2], [yb - ylim * .012, yb, yb, yb - ylim * .012], color='#202020', lw=0.7)
            ax.text((x1 + x2) / 2, yb + ylim * .005, ptxt, ha='center', va='bottom', fontsize=6.5)
    ax.set_ylim(0, ylim)


def legend(ax, N, knob_mb, ncol=4, fs=7):
    handles = [Patch(facecolor=COLOR[r], edgecolor='none') for r in REGIONS] + [
        Line2D([], [], marker='o', ls='none', ms=3, mfc='#1a1a1a', mec='white', mew=0.2),
        Line2D([], [], color='#202020', lw=0.8, marker='_', ms=5)]
    names = [f'Knob180 chr7+8 ({knob_mb:.1f} Mb)', NULL_LABEL('knob180', N), 'Centromere (tentative)',
             NULL_LABEL('centromere', N), 'Intergenic', NULL_LABEL('intergenic', N), 'MA line', '95% CI']
    ax.legend(handles, names, loc='upper left', ncol=ncol, frameon=False, fontsize=fs, handlelength=1.0,
              handletextpad=0.35, columnspacing=0.7, bbox_to_anchor=(0.01, 1.0), borderaxespad=0)


def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_axisbelow(True); ax.grid(axis='y', linewidth=0.5, alpha=0.25)


def plot_fold(res, idf, N, knob_mb, root, prefix):
    setup_rc()
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    group_w = len(REGIONS) * 0.12
    centers = np.arange(len(CLASSES)) * (group_w + 0.10)
    groups, ymax = [], 0
    for i, cls in enumerate(CLASSES):
        xs, tops, txt = draw_group(ax, centers[i], res, idf, cls, 'fold_vs_SNV', N,
                                   lambda y: f'{y:.2f}', star_col='q_BH')
        groups.append((cls, xs, tops, txt)); ymax = max(ymax, max(tops.values()))
    finish_brackets(fig, ax, groups, res, 'p_fold_real_vs_shuffle', N, ymax)
    style(ax)
    n = {c: int(res.loc[res.category == c, 'n_class'].iloc[0]) for c in CLASSES}
    ax.set_xticks(centers, [f'Indels 1–50 bp\n(n = {n["Indels"]:,})', f'SVs >50 bp\n(n = {n["Validated SVs"]:,})'],
                  fontsize=9)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(centers[0] - group_w / 2 - .04, centers[-1] + group_w / 2 + .04)
    ax.set_ylabel('Fold enrichment relative to SNVs\n(class proportion / SNV proportion)', fontsize=9)
    legend(ax, N, knob_mb)
    fig.tight_layout(pad=0.3)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(root / f'{prefix}.fold.{ext}', dpi=600 if ext == 'png' else None,
                    bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)


def plot_density(res, idf, N, knob_mb, root, prefix):
    setup_rc()
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.0))
    names = {'SNVs': 'SNVs', 'Indels': 'Indels 1–50 bp', 'Validated SVs': 'SVs >50 bp'}
    for ax, cls in zip(axes, ALL):
        xs, tops, txt = draw_group(ax, 0, res, idf, cls, 'per_Mb_per_line', N,
                                   lambda y: f'{y:.3g}' if y < 10 else f'{y:.1f}')
        finish_brackets(fig, ax, [(cls, xs, tops, txt)], res, 'p_density_real_vs_shuffle', N, max(tops.values()))
        style(ax)
        ax.set_xticks([])
        ax.set_xlim(-0.4, 0.4)
        n = int(res.loc[res.category == cls, 'n_class'].iloc[0])
        ax.set_title(f'{names[cls]} (n = {n:,})', fontsize=9)
        ax.set_ylabel(f'{names[cls].split()[0]} per Mb per MA line', fontsize=8.5)
    handles = [Patch(facecolor=COLOR[r], edgecolor='none') for r in REGIONS] + [
        Line2D([], [], marker='o', ls='none', ms=3, mfc='#1a1a1a', mec='white', mew=0.2),
        Line2D([], [], color='#202020', lw=0.8, marker='_', ms=5)]
    labels = [f'Knob180 chr7+8 ({knob_mb:.1f} Mb)', NULL_LABEL('knob180', N), 'Centromere (tentative)',
              NULL_LABEL('centromere', N), 'Intergenic', NULL_LABEL('intergenic', N), 'MA line', '95% CI']
    fig.legend(handles, labels, loc='lower center', ncol=8, frameon=False, fontsize=7.5,
               handlelength=1.0, handletextpad=0.35, columnspacing=0.9, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=[0, 0.06, 1, 1], pad=0.4)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(root / f'{prefix}.density.{ext}', dpi=600 if ext == 'png' else None,
                    bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-dir', default=str(HERE))
    ap.add_argument('--knob-gff', default=str(HERE.parent / 'J702_flip.fa.mod.EDTA.TEanno.gff3'))
    ap.add_argument('--gene-gff', default=str(HERE / 'J702.liftoff.gff3'))
    ap.add_argument('--bad-bins', default='/scratch/ak18787/B73_project/vcf_comparison_re_after_normalization/'
                    'vcf_to_compare/three_callers_private_500kb/syri/bad_500kb_bins/*.bad_500kb_bins.bed')
    ap.add_argument('--knob-merge-kb', type=float, default=100.0)
    ap.add_argument('--n-shuffle', type=int, default=1000)
    ap.add_argument('--bedtools', default='bedtools')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 1)))
    ap.add_argument('--work-dir', default=str(HERE / 'badbin_knob'))
    ap.add_argument('--out-prefix', default='j702_region_badbin_knob')
    ap.add_argument('--plot-only', action='store_true')
    a = ap.parse_args()
    root, wd = Path(a.input_dir), Path(a.work_dir)
    if a.plot_only:
        res = pd.read_csv(root / f'{a.out_prefix}.tsv', sep='\t')
        idf = pd.read_csv(root / f'{a.out_prefix}.per_line.tsv', sep='\t')
        N = len(pd.read_csv(root / f'{a.out_prefix}.null.tsv', sep='\t'))
        kmb = float(res.loc[res.region == 'Knob180', 'region_bp'].iloc[0]) / 1e6
        plot_fold(res, idf, N, kmb, root, a.out_prefix); plot_density(res, idf, N, kmb, root, a.out_prefix)
        return
    wd.mkdir(exist_ok=True); sdir = wd / 'shuffles'; sdir.mkdir(exist_ok=True)

    # ---- regions ----
    bad = bad_bins(a.bad_bins)
    gaps = {c: perm.non_gap_intervals(c)[1] for c in LENS}
    excl = {c: merge(list(bad[c]) + list(gaps[c])) for c in LENS}
    kd = defaultdict(list)
    for c, s, e in raw_knob(a.knob_gff):
        kd[c].append((s, e))
    knob = {c: (knob_blocks(kd[c], excl[c], a.knob_merge_kb * 1e3) if c in KNOB_CHROMS else []) for c in LENS}
    gapm = {c: merge(list(gaps[c])) for c in LENS}
    cen = {c: subtract([CENTROMERES[c]], gapm[c]) for c in LENS}
    knob_part = {c: subtract(knob[c], cen[c]) for c in LENS}          # knob not in centromere
    gen = reg.genes(a.gene_gff)
    ig = {}
    for c in LENS:
        block = merge(list(gapm[c]) + cen[c] + knob[c] + list(gen[c]))
        ig[c] = subtract([(0, LENS[c])], block)
    region_bp = {'Knob180': sum(bp(knob_part[c]) for c in LENS), 'Centromere': sum(bp(cen[c]) for c in LENS),
                 'Intergenic': sum(bp(ig[c]) for c in LENS)}
    access_bp = sum(LENS[c] - bp(gapm[c]) for c in LENS)

    (wd / 'genome.txt').write_text(''.join(f'{c}\t{n}\n' for c, n in LENS.items()))
    (wd / 'excluded_gaps_badbins.bed').write_text(''.join(f'{c}\t{s}\t{e}\n' for c in LENS for s, e in excl[c]))
    beds = {'Knob180': wd / 'knob180_chr7_chr8_blocks.bed', 'Centromere': wd / 'centromere.bed',
            'Intergenic': wd / 'intergenic_largest_first.bed'}
    beds['Knob180'].write_text(''.join(f'{c}\t{s}\t{e}\n' for c in LENS for s, e in knob[c]))
    beds['Centromere'].write_text(''.join(f'{c}\t{s}\t{e}\n' for c in LENS for s, e in cen[c]))
    igl = sorted(((c, s, e) for c in LENS for s, e in ig[c]), key=lambda r: r[1] - r[2])
    beds['Intergenic'].write_text(''.join(f'{c}\t{s}\t{e}\n' for c, s, e in igl))

    feat = {'Centromere': {c: as_arrays(cen[c]) for c in LENS},
            'Knob180': {c: as_arrays(knob[c]) for c in LENS},
            'Genic': {c: as_arrays(gen[c]) for c in LENS}}
    excl_arr = {c: as_arrays(gapm[c]) for c in LENS}       # events: drop AGP-gap midpoints only
    ig_arr = {c: as_arrays(ig[c]) for c in LENS}

    # ---- events as midpoints; drop those in AGP gaps (none expected) ----
    events = reg.load_events(root)
    mid, dropped = {}, {}
    for cls in ALL:
        d = {c: ([], []) for c in LENS}
        for sample, c, s, e, _ in events[cls]:
            d[c][0].append((s + e) // 2 if e - s > 1 else s); d[c][1].append(sample)
        mid[cls], dropped[cls] = {}, 0
        for c in LENS:
            m, smp = np.array(d[c][0], np.int64), np.array(d[c][1])
            k = ~hit(m, m + 1, excl_arr[c])
            dropped[cls] += int((~k).sum())
            mid[cls][c] = (m[k], smp[k])

    counts, totals = defaultdict(int), defaultdict(int)
    for cls in ALL:
        for c in LENS:
            m, smp = mid[cls][c]
            left = np.ones(len(m), bool)
            for r in ['Centromere', 'Knob180', 'Genic']:
                h = hit(m, m + 1, feat[r][c]) & left
                for l in SAMPLES:
                    counts[(cls, l, r)] += int((h & (smp == l)).sum())
                left &= ~h
            for l in SAMPLES:
                counts[(cls, l, 'Intergenic')] += int((left & (smp == l)).sum())
                totals[(cls, l)] += int((smp == l).sum())
            assert (hit(m, m + 1, ig_arr[c]) == left).all(), 'intergenic partition mismatch'

    # ---- shuffles ----
    global _W
    _W = dict(bedtools=a.bedtools, wd=wd, sdir=sdir, beds=beds, mid=mid, cen=cen)
    rows = []
    with mp.get_context('fork').Pool(a.threads) as pool:
        for n, row in enumerate(pool.imap_unordered(one_shuffle, range(1, a.n_shuffle + 1)), 1):
            rows.append(row)
            if n % 100 == 0:
                print(f'shuffle {n}/{a.n_shuffle}', flush=True)
    null = pd.DataFrame(rows).sort_values('seed')
    null.to_csv(root / f'{a.out_prefix}.null.tsv', sep='\t', index=False)
    N = a.n_shuffle
    n_in = {r: sum(1 for _ in open(beds[r])) for r in REAL}

    # ---- statistics ----
    def T(cls):
        return sum(totals[(cls, l)] for l in SAMPLES)
    res, indiv = [], []
    ns = T('SNVs')
    for cls in ALL:
        n = T(cls)
        for r in REAL:
            Rbp = region_bp[r]
            k = sum(counts[(cls, l, r)] for l in SAMPLES)
            ks = sum(counts[('SNVs', l, r)] for l in SAMPLES)
            obs_fold = (k / n) / (ks / ns)
            obs_den = k / len(SAMPLES) / (Rbp / 1e6)
            p_f = fisher_exact([[k, n - k], [ks, ns - ks]])[1] if cls != 'SNVs' else np.nan
            res.append(dict(category=cls, region=r, n_class=n, class_in_region=k, n_SNV=ns, SNV_in_region=ks,
                            region_bp=Rbp, fold_vs_SNV=obs_fold, per_Mb_per_line=obs_den, p_fisher=p_f))
            for l in SAMPLES:
                k1, n1, ks1, ns1 = counts[(cls, l, r)], totals[(cls, l)], counts[('SNVs', l, r)], totals[('SNVs', l)]
                indiv.append(dict(category=cls, region=r, MA_line=l, class_in_region=k1, n_class=n1,
                                  fold_vs_SNV=(k1 / n1) / (ks1 / ns1) if ks1 and n1 else np.nan,
                                  per_Mb_per_line=k1 / (Rbp / 1e6)))
            # shuffled
            kc = sum(null[f'{r}|{cls}|{l}'] for l in SAMPLES)
            sc = sum(null[f'{r}|SNVs|{l}'] for l in SAMPLES)
            sbp = null[f'{r}|bp']
            fold_k = ((kc / n) / (sc / ns)).replace([np.inf], np.nan)
            den_k = kc / len(SAMPLES) / (sbp / 1e6)
            res.append(dict(category=cls, region=SHUF[r], n_class=n, class_in_region=kc.mean(), n_SNV=ns,
                            SNV_in_region=sc.mean(), region_bp=sbp.mean(),
                            fold_vs_SNV=(kc.mean() / n) / (sc.mean() / ns), per_Mb_per_line=den_k.mean(),
                            fold_null_q2_5=fold_k.quantile(.025), fold_null_q97_5=fold_k.quantile(.975),
                            density_null_q2_5=den_k.quantile(.025), density_null_q97_5=den_k.quantile(.975),
                            p_fold_real_vs_shuffle=emp_p(obs_fold, fold_k) if cls != 'SNVs' else np.nan,
                            p_density_real_vs_shuffle=emp_p(obs_den, den_k),
                            intervals_placed_min=int(null[f'{r}|placed'].min()), intervals_in=n_in[r]))
            for l in SAMPLES:
                k1 = (null[f'{r}|{cls}|{l}']).mean(); ks1 = null[f'{r}|SNVs|{l}'].mean()
                n1, ns1 = totals[(cls, l)], totals[('SNVs', l)]
                indiv.append(dict(category=cls, region=SHUF[r], MA_line=l, class_in_region=k1, n_class=n1,
                                  fold_vs_SNV=(k1 / n1) / (ks1 / ns1) if ks1 and n1 else np.nan,
                                  per_Mb_per_line=(null[f'{r}|{cls}|{l}'] / (sbp / 1e6)).mean()))
    res = pd.DataFrame(res)
    ok = res.p_fisher.notna()
    res.loc[ok, 'q_BH'] = bh(res.loc[ok, 'p_fisher'])
    idf = pd.DataFrame(indiv)
    res.to_csv(root / f'{a.out_prefix}.tsv', sep='\t', index=False)
    idf.to_csv(root / f'{a.out_prefix}.per_line.tsv', sep='\t', index=False)

    kmb = region_bp['Knob180'] / 1e6
    plot_fold(res, idf, N, kmb, root, a.out_prefix)
    plot_density(res, idf, N, kmb, root, a.out_prefix)

    pd.set_option('display.width', 250)
    print(f'SyRI bad bins (union) {sum(bp(bad[c]) for c in LENS)/1e6:.1f} Mb; non-gap genome {access_bp/1e6:,.1f} Mb')
    print(f'knob180 blocks chr7+8: {sum(len(knob[c]) for c in LENS)} ({sum(bp(knob[c]) for c in LENS)/1e6:.2f} Mb); '
          f'centromere {region_bp["Centromere"]/1e6:.2f} Mb; intergenic {region_bp["Intergenic"]/1e6:,.1f} Mb')
    print('events dropped (midpoint in AGP gap):', dropped)
    print('intervals in / min placed:', {r: (n_in[r], int(null[f"{r}|placed"].min())) for r in REAL})
    print(res.drop(columns=['n_SNV']).to_string(index=False))


if __name__ == '__main__':
    main()
