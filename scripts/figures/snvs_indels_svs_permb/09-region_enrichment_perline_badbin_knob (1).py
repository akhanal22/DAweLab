#!/usr/bin/env python3
"""Same analysis as 08-region_enrichment_badbin_knob.py, but bad bins are applied
PER MA LINE instead of as the five-line union.

For each line l (B73D1, B73D2, B73D4, B73D7c, B73D10):
  * bad bins = that line's own SyRI bad 500 kb bins (+ AGP U gaps)
  * knob180_l = chr7 + chr8 EDTA knob180 hits; hits overlapping l's bad bins / gaps
    dropped; remaining hits merged when <= --knob-merge-kb apart, never across one
    of l's bad bins / gaps
  * centromere = tentative spans (same for all lines)
  * intergenic_l = genome - gaps - centromere - knob180_l - genes
  * line l's events (midpoints) are partitioned on l's own regions
  * shuffles: for every seed i, l's knob180 blocks, the centromeres and l's
    intergenic intervals are shuffled (bedtools shuffle -chrom -noOverlapping
    -excl gaps + l's bad bins -seed i); only l's events are counted against them.
--genome-wide: the same, but every region (knob180, centromere, intergenic) is
  shuffled across all of chr1-10 (bedtools shuffle without -chrom), so all three
  controls are "same-size region placed anywhere in the genome".
--knob-piece-mb X: for the shuffle only, knob180 blocks are cut into consecutive
  pieces of <= X Mb (same total bp) so that large blocks are not kept away from
  chromosome ends by the need to fit whole; the real knob region is unchanged.
Pooling: for each seed, counts and region bp are summed over the five lines.
  density (pooled) = sum_l events_l / sum_l Mb_l      (events per Mb per line)
  fold             = (class in R / class total) / (SNV in R / SNV total)
  p (real vs shuffled) = two-sided empirical, (1 + #{|x_k - m| >= |obs - m|}) / (N + 1)
  '*' on real fold bars = Fisher exact class vs SNV (pooled), BH over 6 tests.
Dots = per-line values on that line's own regions.
"""
import argparse
import glob
import importlib.util
import multiprocessing as mp
import os
import subprocess
from collections import defaultdict
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('s08', HERE / '08-region_enrichment_badbin_knob.py')
s08 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s08)
LENS, merge, SAMPLES, bh, reg, perm = s08.LENS, s08.merge, s08.SAMPLES, s08.bh, s08.reg, s08.perm
CENTROMERES, raw_knob, as_arrays, hit = s08.CENTROMERES, s08.raw_knob, s08.as_arrays, s08.hit
subtract, bp, knob_blocks, emp_p = s08.subtract, s08.bp, s08.knob_blocks, s08.emp_p
REAL, SHUF, ALL, EXTRA, KNOB_CHROMS = s08.REAL, s08.SHUF, s08.ALL, s08.EXTRA, s08.KNOB_CHROMS

SYRI_BINS = ('/scratch/ak18787/B73_project/vcf_comparison_re_after_normalization/vcf_to_compare/'
             'three_callers_private_500kb/syri/bad_500kb_bins/{file_line}.bad_500kb_bins.bed')
FILE_LINE = {'B73D7c': 'B73D7c1'}          # bad-bin file name differs for this line

_W = {}


def one_shuffle(i):
    w = _W
    row = {'seed': i}
    for l in SAMPLES:
        for r in REAL:
            p = subprocess.run([w['bedtools'], 'shuffle', '-i', str(w['beds'][l][r]), '-g', str(w['wd'] / 'genome.txt'),
                                '-excl', str(w['excl_bed'][l]), '-seed', str(i)] + w['chrom_flag'] + EXTRA[r],
                               check=True, capture_output=True, text=True)
            if i <= 2:
                (w['wd'] / l / f'{r.lower()}_shuffle_seed{i}.bed').write_text(p.stdout)
            df = pd.read_csv(StringIO(p.stdout), sep='\t', header=None, usecols=[0, 1, 2]) if p.stdout else \
                pd.DataFrame(columns=[0, 1, 2])
            row[f'{r}|placed|{l}'] = len(df)
            tot = 0
            for cls in ALL:
                row[f'{r}|{cls}|{l}'] = 0
            for c, g in df.groupby(0):
                iv = merge(zip(g[1].to_numpy(), g[2].to_numpy()))
                if r == 'Knob180':
                    iv = subtract(iv, w['cen'][c])
                tot += bp(iv)
                arr = as_arrays(iv)
                for cls in ALL:
                    m = w['mid'][cls][l][c]
                    row[f'{r}|{cls}|{l}'] += int(hit(m, m + 1, arr).sum())
            row[f'{r}|bp|{l}'] = tot
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-dir', default=str(HERE))
    ap.add_argument('--knob-gff', default=str(HERE.parent / 'J702_flip.fa.mod.EDTA.TEanno.gff3'))
    ap.add_argument('--gene-gff', default=str(HERE / 'J702.liftoff.gff3'))
    ap.add_argument('--bad-bins', default=SYRI_BINS, help='{file_line} is substituted')
    ap.add_argument('--knob-merge-kb', type=float, default=100.0)
    ap.add_argument('--n-shuffle', type=int, default=1000)
    ap.add_argument('--bedtools', default='bedtools')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 1)))
    ap.add_argument('--work-dir', default=str(HERE / 'perline_badbin_knob'))
    ap.add_argument('--out-prefix', default='j702_region_perline_badbin_knob')
    ap.add_argument('--plot-only', action='store_true')
    ap.add_argument('--genome-wide', action='store_true',
                    help='shuffle every region across all of chr1-10 (no bedtools -chrom); '
                         'default output prefix/work dir get a _genomewide suffix')
    ap.add_argument('--knob-piece-mb', type=float, default=0,
                    help='if > 0, cut knob180 blocks into pieces of at most this many Mb '
                         'for shuffling only (real knob unchanged); adds _pieces<X>Mb suffix')
    a = ap.parse_args()
    if a.genome_wide:
        if a.out_prefix == 'j702_region_perline_badbin_knob':
            a.out_prefix += '_genomewide'
        if a.work_dir == str(HERE / 'perline_badbin_knob'):
            a.work_dir += '_genomewide'
    if a.knob_piece_mb > 0:
        suf = f'_pieces{a.knob_piece_mb:g}Mb'
        if not a.out_prefix.endswith(suf) and a.out_prefix.startswith('j702_region_perline_badbin_knob'):
            a.out_prefix += suf
        if a.work_dir.startswith(str(HERE / 'perline_badbin_knob')):
            a.work_dir += suf
    root, wd = Path(a.input_dir), Path(a.work_dir)

    def knob_label_mb(res):
        return float(res.loc[res.region == 'Knob180', 'region_bp'].iloc[0]) / len(SAMPLES) / 1e6

    if a.plot_only:
        res = pd.read_csv(root / f'{a.out_prefix}.tsv', sep='\t')
        idf = pd.read_csv(root / f'{a.out_prefix}.per_line.tsv', sep='\t')
        N = len(pd.read_csv(root / f'{a.out_prefix}.null.tsv', sep='\t'))
        s08.plot_fold(res, idf, N, knob_label_mb(res), root, a.out_prefix)
        s08.plot_density(res, idf, N, knob_label_mb(res), root, a.out_prefix)
        return
    wd.mkdir(exist_ok=True)

    # ---- shared pieces ----
    gaps = {c: merge(list(perm.non_gap_intervals(c)[1])) for c in LENS}
    cen = {c: subtract([CENTROMERES[c]], gaps[c]) for c in LENS}
    gen = reg.genes(a.gene_gff)
    kd = defaultdict(list)
    for c, s, e in raw_knob(a.knob_gff):
        kd[c].append((s, e))
    (wd / 'genome.txt').write_text(''.join(f'{c}\t{n}\n' for c, n in LENS.items()))
    cen_bed = wd / 'centromere.bed'
    cen_bed.write_text(''.join(f'{c}\t{s}\t{e}\n' for c in LENS for s, e in cen[c]))

    # ---- per-line regions ----
    bad, excl, knob, ig, rbp, beds, excl_bed = {}, {}, {}, {}, {}, {}, {}
    for l in SAMPLES:
        (wd / l).mkdir(exist_ok=True)
        d = defaultdict(list)
        f = a.bad_bins.format(file_line=FILE_LINE.get(l, l))
        for line in open(f):
            x = line.split()
            if len(x) >= 3 and x[0] in LENS:
                d[x[0]].append((int(x[1]), int(x[2])))
        bad[l] = {c: merge(d[c]) for c in LENS}
        excl[l] = {c: merge(bad[l][c] + gaps[c]) for c in LENS}
        knob[l] = {c: (knob_blocks(kd[c], excl[l][c], a.knob_merge_kb * 1e3) if c in KNOB_CHROMS else [])
                   for c in LENS}
        ig[l] = {c: subtract([(0, LENS[c])], merge(gaps[c] + cen[c] + knob[l][c] + list(gen[c]))) for c in LENS}
        rbp[l] = {'Knob180': sum(bp(subtract(knob[l][c], cen[c])) for c in LENS),
                  'Centromere': sum(bp(cen[c]) for c in LENS),
                  'Intergenic': sum(bp(ig[l][c]) for c in LENS)}
        excl_bed[l] = wd / l / 'excluded_gaps_badbins.bed'
        excl_bed[l].write_text(''.join(f'{c}\t{s}\t{e}\n' for c in LENS for s, e in excl[l][c]))
        kb = wd / l / 'knob180_chr7_chr8_blocks.bed'
        kb.write_text(''.join(f'{c}\t{s}\t{e}\n' for c in LENS for s, e in knob[l][c]))
        if a.knob_piece_mb > 0:                     # shuffle input: blocks cut into pieces
            step = int(a.knob_piece_mb * 1e6)
            pieces = [(c, x, min(x + step, e)) for c in LENS for s, e in knob[l][c] for x in range(s, e, step)]
            kb = wd / l / f'knob180_shuffle_pieces_{a.knob_piece_mb:g}Mb.bed'
            kb.write_text(''.join(f'{c}\t{s}\t{e}\n' for c, s, e in pieces))
        igb = wd / l / 'intergenic_largest_first.bed'
        igl = sorted(((c, s, e) for c in LENS for s, e in ig[l][c]), key=lambda r: r[1] - r[2])
        igb.write_text(''.join(f'{c}\t{s}\t{e}\n' for c, s, e in igl))
        beds[l] = {'Knob180': kb, 'Centromere': cen_bed, 'Intergenic': igb}

    # ---- events: midpoints per line ----
    events = reg.load_events(root)
    mid = {cls: {l: {c: [] for c in LENS} for l in SAMPLES} for cls in ALL}
    for cls in ALL:
        for sample, c, s, e, _ in events[cls]:
            mid[cls][sample][c].append((s + e) // 2 if e - s > 1 else s)
    gap_arr = {c: as_arrays(gaps[c]) for c in LENS}
    for cls in ALL:
        for l in SAMPLES:
            for c in LENS:
                m = np.array(mid[cls][l][c], np.int64)
                mid[cls][l][c] = m[~hit(m, m + 1, gap_arr[c])]

    # ---- observed partition, per line on its own regions ----
    counts, totals, own_bad = defaultdict(int), defaultdict(int), defaultdict(int)
    for l in SAMPLES:
        feat = {'Centromere': {c: as_arrays(cen[c]) for c in LENS},
                'Knob180': {c: as_arrays(knob[l][c]) for c in LENS},
                'Genic': {c: as_arrays(gen[c]) for c in LENS}}
        for cls in ALL:
            for c in LENS:
                m = mid[cls][l][c]
                own_bad[(cls, l)] += int(hit(m, m + 1, as_arrays(bad[l][c])).sum())
                left = np.ones(len(m), bool)
                for r in ['Centromere', 'Knob180', 'Genic']:
                    h = hit(m, m + 1, feat[r][c]) & left
                    counts[(cls, l, r)] += int(h.sum()); left &= ~h
                counts[(cls, l, 'Intergenic')] += int(left.sum())
                totals[(cls, l)] += len(m)
                assert (hit(m, m + 1, as_arrays(ig[l][c])) == left).all()

    # ---- shuffles ----
    global _W
    _W = dict(bedtools=a.bedtools, wd=wd, beds=beds, excl_bed=excl_bed, mid=mid, cen=cen,
              chrom_flag=[] if a.genome_wide else ['-chrom'])
    rows = []
    with mp.get_context('fork').Pool(a.threads) as pool:
        for n, row in enumerate(pool.imap_unordered(one_shuffle, range(1, a.n_shuffle + 1)), 1):
            rows.append(row)
            if n % 100 == 0:
                print(f'shuffle {n}/{a.n_shuffle}', flush=True)
    null = pd.DataFrame(rows).sort_values('seed')
    null.to_csv(root / f'{a.out_prefix}.null.tsv', sep='\t', index=False)
    N = a.n_shuffle

    # ---- statistics ----
    def T(cls):
        return sum(totals[(cls, l)] for l in SAMPLES)
    ns = T('SNVs')
    res, indiv = [], []
    for cls in ALL:
        n = T(cls)
        for r in REAL:
            k = sum(counts[(cls, l, r)] for l in SAMPLES)
            ks = sum(counts[('SNVs', l, r)] for l in SAMPLES)
            Rbp = sum(rbp[l][r] for l in SAMPLES)                  # line-bp summed over lines
            obs_fold = (k / n) / (ks / ns)
            obs_den = k / (Rbp / 1e6)
            p_f = fisher_exact([[k, n - k], [ks, ns - ks]])[1] if cls != 'SNVs' else np.nan
            res.append(dict(category=cls, region=r, n_class=n, class_in_region=k, n_SNV=ns, SNV_in_region=ks,
                            region_bp=Rbp, fold_vs_SNV=obs_fold, per_Mb_per_line=obs_den, p_fisher=p_f))
            for l in SAMPLES:
                k1, n1, ks1, ns1 = counts[(cls, l, r)], totals[(cls, l)], counts[('SNVs', l, r)], totals[('SNVs', l)]
                indiv.append(dict(category=cls, region=r, MA_line=l, class_in_region=k1, n_class=n1,
                                  region_bp=rbp[l][r], fold_vs_SNV=(k1 / n1) / (ks1 / ns1) if ks1 and n1 else np.nan,
                                  per_Mb_per_line=k1 / (rbp[l][r] / 1e6)))
            kc = sum(null[f'{r}|{cls}|{l}'] for l in SAMPLES)
            sc = sum(null[f'{r}|SNVs|{l}'] for l in SAMPLES)
            sbp = sum(null[f'{r}|bp|{l}'] for l in SAMPLES)
            fold_k = ((kc / n) / (sc / ns)).replace([np.inf], np.nan)
            den_k = kc / (sbp / 1e6)
            placed = min(int(null[f'{r}|placed|{l}'].min()) for l in SAMPLES)
            res.append(dict(category=cls, region=SHUF[r], n_class=n, class_in_region=kc.mean(), n_SNV=ns,
                            SNV_in_region=sc.mean(), region_bp=sbp.mean(),
                            fold_vs_SNV=(kc.mean() / n) / (sc.mean() / ns), per_Mb_per_line=den_k.mean(),
                            fold_null_q2_5=fold_k.quantile(.025), fold_null_q97_5=fold_k.quantile(.975),
                            density_null_q2_5=den_k.quantile(.025), density_null_q97_5=den_k.quantile(.975),
                            p_fold_real_vs_shuffle=emp_p(obs_fold, fold_k) if cls != 'SNVs' else np.nan,
                            p_density_real_vs_shuffle=emp_p(obs_den, den_k), intervals_placed_min_any_line=placed))
            for l in SAMPLES:
                k1 = null[f'{r}|{cls}|{l}']; ks1 = null[f'{r}|SNVs|{l}']
                n1, ns1 = totals[(cls, l)], totals[('SNVs', l)]
                indiv.append(dict(category=cls, region=SHUF[r], MA_line=l, class_in_region=k1.mean(), n_class=n1,
                                  region_bp=null[f'{r}|bp|{l}'].mean(),
                                  fold_vs_SNV=(k1.mean() / n1) / (ks1.mean() / ns1) if ks1.mean() and n1 else np.nan,
                                  per_Mb_per_line=(k1 / (null[f'{r}|bp|{l}'] / 1e6)).mean()))
    res = pd.DataFrame(res)
    ok = res.p_fisher.notna()
    res.loc[ok, 'q_BH'] = bh(res.loc[ok, 'p_fisher'])
    idf = pd.DataFrame(indiv)
    res.to_csv(root / f'{a.out_prefix}.tsv', sep='\t', index=False)
    idf.to_csv(root / f'{a.out_prefix}.per_line.tsv', sep='\t', index=False)
    s08.plot_fold(res, idf, N, knob_label_mb(res), root, a.out_prefix)
    s08.plot_density(res, idf, N, knob_label_mb(res), root, a.out_prefix)

    pd.set_option('display.width', 250)
    for l in SAMPLES:
        print(f'{l}: bad bins {sum(bp(bad[l][c]) for c in LENS)/1e6:.1f} Mb; knob180 blocks '
              f'{sum(len(knob[l][c]) for c in LENS)} ({rbp[l]["Knob180"]/1e6:.2f} Mb); '
              f'intergenic {rbp[l]["Intergenic"]/1e6:,.1f} Mb; events in own bad bins: '
              f'{ {cls: own_bad[(cls, l)] for cls in ALL} }')
    print(res.drop(columns=['n_SNV']).to_string(index=False))


if __name__ == '__main__':
    main()
