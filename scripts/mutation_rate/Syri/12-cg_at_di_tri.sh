#!/bin/bash
#SBATCH --job-name=syri_repeat_overlap
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
set -euo pipefail
module load BCFtools/1.23.1-GCC-13.3.0

# Run from the directory containing vcf_to_compare.
# Optional first positional argument overrides INPUT_ROOT.
REF="${REF:-/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa}"
INPUT_ROOT="${1:-${INPUT_ROOT:-vcf_to_compare/three_callers_private_500kb}}"
SYRI="${SYRI:-${INPUT_ROOT}/syri/syri.private.500kb.filtered.snps_indels.merged.vcf.gz}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MIN_HP=7
HP_FLANK=0 # No flank masks or flank-filter stage.
MIN_DINUC_COPIES=5
MIN_TRINUC_COPIES=5
MAX_REPEAT_INDEL_BP=50  # ALL repeat filters apply only to 1-50-bp indels.
# Remove every eligible repeat-overlapping indel; no percentage threshold.

HP_GROUPS=CG,AT # Filter CG first, then AT, then di and tri.
# Reuse unpadded BEDs generated against this REF. Optional second argument is BED_DIR.
BED_DIR="${2:-${BED_DIR:-${INPUT_ROOT}/repeat_filter_comparison/04_hp_flank1_di_tri}}"
# Explicit HP_BED (CG), AT_HP_BED, DI_BED, and TRI_BED overrides are supported.
choose_bed() {
    local prefix="$1" unit="$2" wanted="$3" path
    for path in "$BED_DIR/${prefix}${wanted}${unit}.bed" \
                "$BED_DIR/${prefix}$([[ "$unit" == copies ]] && echo 3 || echo 5)${unit}.bed"; do
        if [[ -f "$path" ]]; then printf '%s\n' "$path"; return; fi
    done
    echo "Missing core BED in $BED_DIR: ${prefix}${wanted}${unit}.bed (or older lower-threshold BED). Set BED_DIR or explicit BED paths." >&2
    return 1
}
HP_BED="${HP_BED:-$(choose_bed reference.CG_homopolymers.min '' "$MIN_HP")}" 
AT_HP_BED="${AT_HP_BED:-$(choose_bed reference.AT_homopolymers.min '' "$MIN_HP")}"
DI_BED="${DI_BED:-$(choose_bed reference.perfect_dinucleotide_repeats.min copies "$MIN_DINUC_COPIES")}" 
TRI_BED="${TRI_BED:-$(choose_bed reference.perfect_trinucleotide_repeats.min copies "$MIN_TRINUC_COPIES")}" 
for repeat_bed in "$HP_BED" "$AT_HP_BED" "$DI_BED" "$TRI_BED"; do
    [[ -f "$repeat_bed" ]] || { echo "Missing BED: $repeat_bed" >&2; exit 1; }
done
[[ -s "$REF" && -s "$SYRI" ]] || { echo 'Missing REF or SYRI input' >&2; exit 1; }
PARENT="${INPUT_ROOT}/repeat_filter_comparison"
mkdir -p "$PARENT"
OUT=$(mktemp -d "${PARENT}/syri_CGhp7_then_AThp7_di5_tri5_alloverlap_$(date +%Y%m%d_%H%M%S)_XXXXXX")
echo "Results: $OUT"
# Optional reuse of a normalized indel VCF; input BEDs are read only.
if [[ -n "${NORMALIZED_VCF:-}" ]]; then
    [[ -s "$NORMALIZED_VCF" ]] || { echo 'Missing NORMALIZED_VCF' >&2; exit 1; }
    cp --reflink=auto "$NORMALIZED_VCF" "$OUT/syri.indels.normalized.vcf.gz"
else
    # Preserve multiallelic rows during normalization, stop on REF mismatches.
    bcftools view --threads "$THREADS" -v indels -Ou "$SYRI" |
    bcftools norm --threads "$THREADS" -f "$REF" -c e -Ou |
    bcftools view -v indels -Ou |
    bcftools norm -d exact -Ou |
    bcftools sort -Oz -o "$OUT/syri.indels.normalized.vcf.gz"
fi
bcftools index "$OUT/syri.indels.normalized.vcf.gz"

python3 - "$REF" "$OUT/syri.indels.normalized.vcf.gz" "$OUT" "$MIN_HP" "$HP_FLANK" "$MIN_DINUC_COPIES" "$MIN_TRINUC_COPIES" "$MAX_REPEAT_INDEL_BP" "$HP_GROUPS" "$HP_BED" "$AT_HP_BED" "$DI_BED" "$TRI_BED" <<'PY'
import bisect, collections, csv, gzip, re, sys
from contextlib import ExitStack
from pathlib import Path

fasta, vcf, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
min_hp, flank, min_di, min_tri, max_di_tri = map(int, sys.argv[4:9])
hp_groups = set(sys.argv[9].split(','))
if not hp_groups or not hp_groups <= {'CG', 'AT'}:
    raise ValueError('HP_GROUPS must be CG, AT, or CG,AT')


def make_index(runs):
    # Prefix maxima support overlapping padded masks without scanning chromosome.
    runs.sort()
    starts, ends, maximum = [], [], -1
    for a, b, *_ in runs:
        starts.append(a)
        maximum = max(maximum, b)
        ends.append(maximum)
    return starts, ends, runs


def hits(index, start, end):
    # Insertions are points and touch both tract edges; deletions are half-open.
    starts, prefix_end, runs = index
    stop = bisect.bisect_right(starts, start) if start == end else bisect.bisect_left(starts, end)
    i = stop - 1
    result = []
    while i >= 0 and (prefix_end[i] >= start if start == end else prefix_end[i] > start):
        run = runs[i]
        if run[1] >= start if start == end else run[1] > start:
            result.append(run)
        i -= 1
    return result


def event(pos, ref, alt):
    a, r, t = 0, len(ref), len(alt)
    while a < min(r, t) and ref[a] == alt[a]:
        a += 1
    while r > a and t > a and ref[r - 1] == alt[t - 1]:
        r -= 1
        t -= 1
    return pos - 1 + a, pos - 1 + r, ref[a:r], alt[a:t]


# BED coordinates are 0-based half-open FULL, unpadded reference tracts.
# Older min5 HP / min3-copy di/tri BEDs are accepted: shorter rows are skipped.
# Never merge, split, pad, rewrite, or regenerate these tracts.
bed_paths = dict(zip(('hp_core', 'at_core', 'di', 'tri'), sys.argv[10:14]))
indexes = {key: {} for key in bed_paths}
repeat_counts = collections.Counter()
bed_skipped = collections.Counter()
for key, path in bed_paths.items():
    if re.search(r'plus\d+bp', Path(path).name):
        raise ValueError(f'Use an unpadded core BED, not a plus/flank BED: {path}')
    by_chrom = collections.defaultdict(list)
    period = {'hp_core': 1, 'at_core': 1, 'di': 2, 'tri': 3}[key]
    minimum = {'hp_core': min_hp, 'at_core': min_hp, 'di': 2*min_di, 'tri': 3*min_tri}[key]
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f'{path}:{line_number}: expected chrom start end motif; reuse original motif-bearing BED')
            chrom, a, b, motif = fields[0], int(fields[1]), int(fields[2]), fields[3].upper()
            if a < 0 or b <= a:
                raise ValueError(f'{path}:{line_number}: invalid BED interval')
            # Previous padded HP BEDs have 8 columns with core coordinates.
            if key in {'hp_core', 'at_core'} and len(fields) >= 8:
                raise ValueError(f'{path}:{line_number}: padded BED detected; use core BED')
            if len(motif) != period or any(c not in 'ACGT' for c in motif):
                raise ValueError(f'{path}:{line_number}: invalid {key} motif {motif}')
            if key in {'hp_core', 'at_core'} and motif not in ('CG' if key == 'hp_core' else 'AT'):
                bed_skipped[f'{key}_wrong_base_rows'] += 1
                continue
            if period > 1 and (len(set(motif)) == 1 or (b-a) % period):
                raise ValueError(f'{path}:{line_number}: expected whole copies of a di/tri motif')
            if b-a < minimum:
                bed_skipped[f'{key}_below_minimum_rows'] += 1
                continue
            by_chrom[chrom].append((a, b, a, b, motif))
    for chrom, runs in by_chrom.items():
        # Identical rows must not multiply audit hits; overlapping tracts stay separate.
        runs = list(set(runs))
        indexes[key][chrom] = make_index(runs)
        repeat_counts['CG' if key == 'hp_core' else 'AT' if key == 'at_core' else key] += len(runs)
    print(f'Loaded existing {key} BED: {path}; eligible tracts={repeat_counts["CG" if key == "hp_core" else "AT" if key == "at_core" else key]}', flush=True)
# A contig can legitimately have no repeat BED rows. When available, use the
# reference FASTA index, rather than BED membership, to check VCF contig names.
fai = Path(fasta + '.fai')
reference_contigs = {line.split('\t')[0] for line in fai.read_text().splitlines()} if fai.exists() else None

stages = ('input', 'after_HP_core', 'after_AT_core', 'after_di_1_50', 'after_tri_1_50')
count, bases = collections.Counter(), collections.Counter()
info_numbers = {}
info_header = re.compile(r'^##INFO=<ID=([^,>]+),Number=([^,>]+)')


def selected_info(info, number, total):
    kept = []
    for entry in info.split(';'):
        if entry == '.':
            continue
        name, sep, value = entry.partition('=')
        kind = info_numbers.get(name)
        if name in {'AC','AN','AF','MLEAC','MLEAF'} or kind == 'G' or (total > 1 and name == 'END'):
            continue
        if sep and kind in {'A','R'}:
            vals = value.split(',')
            if len(vals) != total + (kind == 'R'):
                continue
            value = vals[number-1] if kind == 'A' else vals[0]+','+vals[number]
        kept.append(name + ('=' + value if sep else ''))
    return ';'.join(kept) or '.'


def selected_gt(gt, number):
    return ''.join(x if x in {'/','|'} else '1' if x == str(number) else '0' if x == '0' else '.'
                   for x in re.split(r'([/|])', gt))


with ExitStack() as stack:
    source = stack.enter_context(gzip.open(vcf, 'rt'))
    kept = {s: stack.enter_context((out / f'syri.{s}.vcf').open('w')) for s in stages}
    excluded = {s: stack.enter_context((out / f'syri.removed_at_{s}.vcf').open('w')) for s in stages[1:]}
    audit = csv.writer(stack.enter_context((out / 'repeat_fraction_audit.tsv').open('w')), delimiter='\t')
    audit.writerow(['sample','CHROM','POS','ID','ALT_index','type','event_start_0based','event_end_0based',
                    'size_bp','repeat_class','repeat_start_0based','repeat_end_0based','motif',
                    'repeat_length_bp','remove_by_this_repeat_match','first_removed_stage'])
    decisions = csv.writer(stack.enter_context((out / 'variant_decisions.tsv').open('w')), delimiter='\t')
    decisions.writerow(['sample','CHROM','POS','ID','ALT_index','type','size_bp','assessment','first_removed_stage'])
    names = None
    skipped = collections.Counter()
    for raw in source:
        if raw.startswith('#'):
            match = info_header.match(raw)
            if match:
                info_numbers[match[1]] = match[2]
            if raw.startswith('#CHROM'):
                names = raw.rstrip('\n').split('\t')[9:]
                if not names:
                    raise ValueError('No sample columns')
            for h in [*kept.values(), *excluded.values()]:
                h.write(raw)
            continue
        if not raw.strip():
            continue
        if names is None:
            raise ValueError('No sample header')
        f = raw.rstrip('\n').split('\t')
        if len(f) != 9 + len(names):
            raise ValueError('Unexpected VCF column count')
        chrom, pos, ref = f[0], int(f[1]), f[3].upper()
        if reference_contigs is not None and chrom not in reference_contigs:
            raise ValueError(f'Reference contig missing: {chrom}')
        alts, fmt = f[4].split(','), f[8].split(':')
        if 'GT' not in fmt:
            raise ValueError('GT missing')
        gi = fmt.index('GT')
        # Process per sample and carried ALT, not once per overlapping repeat hit.
        for sample_i, (sample, value) in enumerate(zip(names, f[9:])):
            parts = value.split(':')
            gt = parts[gi] if gi < len(parts) else '.'
            numbers = sorted({int(a) for a in re.findall(r'\d+',gt) if 0 < int(a) <= len(alts)})
            for number in numbers:
                alt = alts[number-1].upper()
                if not re.fullmatch('[ACGTN]+', ref) or not re.fullmatch('[ACGTN]+', alt):
                    skipped['symbolic_or_nonsequence_calls'] += 1
                    continue
                if len(ref) == len(alt):
                    continue
                a, b, deleted, inserted = event(pos, ref, alt)
                assessment = 'complex_delins_not_filtered' if deleted and inserted else 'simple_indel_assessed'
                if deleted and inserted:
                    # Net length is not an inserted/deleted sequence length.
                    # Keep complex delins; they are not simple insertion/deletion events.
                    typ, size = ('INS' if len(alt)>len(ref) else 'DEL'), abs(len(alt)-len(ref))
                    skipped['complex_delins_retained_without_repeat_filter'] += 1
                    all_hits = {k: [] for k in indexes}
                else:
                    typ, size = ('INS', len(inserted)) if inserted else ('DEL', len(deleted))
                    all_hits = {k: hits(v.get(chrom, ([],[],[])),a,b) for k,v in indexes.items()}
                flags = []
                for k in ('hp_core','at_core','di','tri'):
                    eligible = 1 <= size <= max_di_tri
                    # Any overlap removes an eligible simple indel, regardless of its size/repeat ratio.
                    flags.append(eligible and bool(all_hits[k]))
                first = next((i+1 for i,flag in enumerate(flags) if flag), None)
                decisions.writerow([sample,chrom,pos,f[2],number,typ,size,assessment,stages[first] if first is not None else '.'])
                gts = ['./.']*len(names)
                gts[sample_i] = selected_gt(gt,number)
                row = '\t'.join([*f[:4],alts[number-1],f[5],f[6],selected_info(f[7],number,len(alts)), 'GT',*gts])+'\n'
                for i,s in enumerate(stages):
                    if first is not None and i >= first:
                        break
                    kept[s].write(row)
                    for size_group in ('all_sizes', '1_50' if size <= max_di_tri else 'gt50'):
                        count[sample,typ,s,size_group] += 1
                        bases[sample,typ,s,size_group] += size
                if first is not None:
                    excluded[stages[first]].write(row)
                for k in ('hp_core','at_core','di','tri'):
                    for _,_,rs,re_,motif in all_hits[k]:
                        eligible = 1 <= size <= max_di_tri
                        audit.writerow([sample,chrom,pos,f[2],number,typ,a,b,size,k,rs,re_,motif,re_-rs,
                                        int(eligible),
                                        stages[first] if first is not None else '.'])

if names is None:
    raise ValueError('No #CHROM header')
with (out/'indel_retention_by_stage.tsv').open('w') as h:
    w = csv.writer(h,delimiter='\t')
    w.writerow(['sample','type','size_group','stage','input_calls','calls_entering_step','removed_at_step','retained_calls',
                'bp_removed_at_step','retained_bp','percent_retained_from_input'])
    for sample in [*names,'ALL_SAMPLES_SUM']:
        ns = names if sample == 'ALL_SAMPLES_SUM' else [sample]
        for typ in ('INS','DEL','ALL_INDELS'):
            ts = ('INS','DEL') if typ == 'ALL_INDELS' else (typ,)
            for size_group in ('all_sizes','1_50','gt50'):
                vals = [sum(count[n,t,s,size_group] for n in ns for t in ts) for s in stages]
                bpvals = [sum(bases[n,t,s,size_group] for n in ns for t in ts) for s in stages]
                for i,s in enumerate(stages):
                    before = vals[i-1] if i else vals[0]
                    removed = before-vals[i]
                    bp_removed = bpvals[i-1]-bpvals[i] if i else 0
                    w.writerow([sample,typ,size_group,s,vals[0],before,removed,vals[i],bp_removed,bpvals[i],
                                f'{100*vals[i]/vals[0]:.2f}' if vals[0] else 'NA'])
                    if typ == 'ALL_INDELS':
                        print(f'{sample}\t{size_group}\t{s}\tremoved={removed}\tretained={vals[i]}')
with (out/'run_parameters.txt').open('w') as h:
    h.write(f'REFERENCE={fasta}\nINPUT_NORMALIZED={vcf}\nHP_MIN_BP={min_hp}\nHP_GROUPS={",".join(sorted(hp_groups))}\nFILTER_ORDER=CG,AT,di,tri\nHP_FLANK_BP={flank}\nALL_REPEAT_FILTERS_MAX_INDEL_BP={max_di_tri}\nDI_MIN_COPIES={min_di}\nTRI_MIN_COPIES={min_tri}\nDI_TRI_MAX_INDEL_BP={max_di_tri}\nREMOVE_RULE=any_reference_repeat_overlap_for_simple_indels_1_50bp\nN3_FILTER=not_applied\nVCF_UNIT=one_row_per_carried_sample_ALT\n')
    for key,path in bed_paths.items():
        h.write(f'EXISTING_BED_{key}={Path(path).resolve()}\n')
    h.write('BED_MODE=existing_unpadded_full_tracts_read_only\nSUMMARY_SIZE_GROUPS=all_sizes,1_50,gt50\n')
    for key,value in bed_skipped.items():
        h.write(f'BED_SKIPPED_{key}={value}\n')
    for key,value in repeat_counts.items():
        h.write(f'REFERENCE_REPEAT_COUNT_{key}={value}\n')
    for key,value in skipped.items():
        h.write(f'{key}={value}\n')
print(f'Summary: {out / "indel_retention_by_stage.tsv"}')
print(f'Audit: {out / "repeat_fraction_audit.tsv"}')
PY

for stage in input after_HP_core after_AT_core after_di_1_50 after_tri_1_50; do
    bcftools view --threads "$THREADS" -Oz -o "$OUT/syri.$stage.vcf.gz" "$OUT/syri.$stage.vcf"
    bcftools index "$OUT/syri.$stage.vcf.gz"
done
for stage in after_HP_core after_AT_core after_di_1_50 after_tri_1_50; do
    bcftools view --threads "$THREADS" -Oz -o "$OUT/syri.removed_at_$stage.vcf.gz" "$OUT/syri.removed_at_$stage.vcf"
    bcftools index "$OUT/syri.removed_at_$stage.vcf.gz"
done
echo "Final VCF: $OUT/syri.after_tri_1_50.vcf.gz"
echo "Step counts: $OUT/indel_retention_by_stage.tsv"
echo "Finished; original inputs and previous outputs unchanged."
