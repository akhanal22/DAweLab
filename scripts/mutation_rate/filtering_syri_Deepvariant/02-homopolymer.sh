#!/bin/bash
#SBATCH --job-name=01_hp_core
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load BCFtools/1.23.1-GCC-13.3.0
module load SAMtools/1.23.1-GCC-13.3.0

REF="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"
THREADS=16
MIN_HP=5

INPUT_ROOT=${1:-vcf_to_compare/three_callers_private_500kb}

SYRI="${INPUT_ROOT}/syri/syri.private.500kb.filtered.snps_indels.merged.vcf.gz"
ANCHORWAVE="${INPUT_ROOT}/anchorwave/anchorwave.private.500kb.filtered.snps_indels.merged.vcf.gz"
DEEPVARIANT="${INPUT_ROOT}/deepvariant/deepvariant.private.500kb.filtered.snps_indels.merged.vcf.gz"

OUT="${INPUT_ROOT}/repeat_filter_comparison/01_hp_core"
mkdir -p "${OUT}"

samtools faidx "${REF}"

###############################################################################
# 1. Verify inputs. All three have already been normalized to the J702 reference.
###############################################################################

for input_vcf in "${SYRI}" "${ANCHORWAVE}" "${DEEPVARIANT}"; do
    [[ -s "${input_vcf}" ]] || { echo "Missing VCF: ${input_vcf}" >&2; exit 1; }
done

###############################################################################
# 2. Keep multiallelic VCF rows intact; normalize and select carried indels
#
# -c e deliberately stops if REF does not agree with the supplied reference.
# Do not replace it with -c s just to force the files through.
###############################################################################

normalize_indels()
{
    local input_vcf=$1
    local output_vcf=$2

# Preserve original ALT lists and genotypes: no bcftools norm -m.
bcftools view --threads "${THREADS}" -v indels -Ou "${input_vcf}" |
bcftools norm --threads "${THREADS}" -f "${REF}" -c e -Ou |
bcftools view -v indels -Ov |
awk 'BEGIN { FS=OFS="\t" }
     /^#/ { print; next }
     {
         gt_col=0; n=split($9, fmt, ":")
         for (k=1; k<=n; k++) if (fmt[k]=="GT") gt_col=k
         if (!gt_col) { print "Missing GT: " $1 ":" $2 > "/dev/stderr"; exit 1 }
         carried=0
         for (i=10; i<=NF; i++) {
             split($i, f, ":"); a=split(f[gt_col], alleles, /[^0-9]+/)
             for (j=1; j<=a; j++)
                 if (alleles[j] ~ /^[0-9]+$/ && alleles[j]+0>0) carried=1
         }
         if (carried) print
     }' |
bcftools view -Ou - |
bcftools norm -d exact -Ou |
bcftools sort -Oz -o "${output_vcf}"

    bcftools index -f "${output_vcf}"
}

normalize_indels \
    "${SYRI}" \
    "${OUT}/syri.indels.normalized.vcf.gz"

normalize_indels \
    "${ANCHORWAVE}" \
    "${OUT}/anchorwave.indels.normalized.vcf.gz"

normalize_indels \
    "${DEEPVARIANT}" \
    "${OUT}/deepvariant.indels.normalized.vcf.gz"

###############################################################################
# 3. Build reference homopolymer BED files and classify overlapping indels
#
# Homopolymers are defined only from the unmodified reference genome. A run
# must already contain at least MIN_HP identical bases in the reference.
# Therefore:
#   reference CCCCC + insertion of A inside the run -> retained as C/G
#   reference CCCC  + insertion of C creating CCCCC -> not retained
#
# BED coordinates are 0-based, half-open. VCF POS remains 1-based. The Python
# code performs the coordinate conversion internally and writes the original
# 1-based VCF record to each output VCF.
#
# For deletions/complex indels, any overlap of the changed reference interval
# with a reference homopolymer is retained. For insertions, the insertion
# boundary must fall inside the run. A same-base insertion at a run boundary
# is also retained because left-normalization places such an insertion at the
# edge of the reference run.
###############################################################################

python3 - \
    "${REF}" "${MIN_HP}" \
    "${OUT}/reference.homopolymers.min${MIN_HP}.bed" \
    "${OUT}/reference.CG_homopolymers.min${MIN_HP}.bed" \
    "${OUT}/reference.AT_homopolymers.min${MIN_HP}.bed" \
    "${OUT}/syri.indels.normalized.vcf.gz" \
    "${OUT}/syri" \
    "${OUT}/anchorwave.indels.normalized.vcf.gz" \
    "${OUT}/anchorwave" \
    "${OUT}/deepvariant.indels.normalized.vcf.gz" \
    "${OUT}/deepvariant" <<'PY'
import gzip
import bisect
import re
import sys

fasta = sys.argv[1]
minimum = int(sys.argv[2])
all_bed = sys.argv[3]
cg_bed = sys.argv[4]
at_bed = sys.argv[5]

def read_fasta(path):
    sequences = {}
    name = None
    pieces = []
    with open(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(pieces).upper()
                name = raw[1:].split()[0]
                pieces = []
            else:
                pieces.append(raw.strip())
    if name is not None:
        sequences[name] = "".join(pieces).upper()
    return sequences

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def trim_shared(ref, alt):
    left = 0
    while left < len(ref) and left < len(alt) and ref[left] == alt[left]:
        left += 1

    ref_end = len(ref)
    alt_end = len(alt)
    while (
        ref_end > left
        and alt_end > left
        and ref[ref_end - 1] == alt[alt_end - 1]
    ):
        ref_end -= 1
        alt_end -= 1

    return left, ref[left:ref_end], alt[left:alt_end]

reference = read_fasta(fasta)

# Each run is stored as (start_0based, end_0based, base, group).
runs_by_chrom = {}
starts_by_chrom = {}
pattern = re.compile(
    "|".join(f"{base}{{{minimum},}}" for base in "ACGT")
)

with open(all_bed, "w") as all_out, \
     open(cg_bed, "w") as cg_out, \
     open(at_bed, "w") as at_out:

    for chrom, sequence in reference.items():
        chromosome_runs = []

        for match in pattern.finditer(sequence):
            start = match.start()
            end = match.end()
            base = sequence[start]
            group = "CG" if base in "CG" else "AT"
            run = (start, end, base, group)
            chromosome_runs.append(run)

            bed_line = (
                f"{chrom}\t{start}\t{end}\t{base}\t{end - start}\t{group}\n"
            )
            all_out.write(bed_line)
            if group == "CG":
                cg_out.write(bed_line)
            else:
                at_out.write(bed_line)

        runs_by_chrom[chrom] = chromosome_runs
        starts_by_chrom[chrom] = [run[0] for run in chromosome_runs]

def overlapping_interval_runs(chrom, start, end):
    """Return reference runs overlapping the 0-based half-open interval."""
    runs = runs_by_chrom.get(chrom, [])
    starts = starts_by_chrom.get(chrom, [])
    if not runs or start >= end:
        return []

    index = max(0, bisect.bisect_right(starts, start) - 1)
    hits = []
    while index < len(runs) and runs[index][0] < end:
        run = runs[index]
        if run[1] > start:
            hits.append(run)
        index += 1
    return hits

def insertion_boundary_runs(chrom, boundary, inserted):
    """Return runs containing an insertion boundary in reference space."""
    runs = runs_by_chrom.get(chrom, [])
    starts = starts_by_chrom.get(chrom, [])
    if not runs:
        return []

    center = bisect.bisect_right(starts, boundary)
    first = max(0, center - 2)
    last = min(len(runs), center + 1)
    hits = []

    for run in runs[first:last]:
        start, end, base, group = run

        # Any inserted sequence is retained when its boundary is strictly
        # inside an existing reference homopolymer.
        inside = start < boundary < end

        # A same-base insertion can be left-normalized to the edge of its run.
        same_base_at_edge = (
            boundary in (start, end)
            and len(inserted) > 0
            and all(character == base for character in inserted)
        )

        if inside or same_base_at_edge:
            hits.append(run)

    return hits

for argument in range(6, len(sys.argv), 2):
    input_vcf = sys.argv[argument]
    output_prefix = sys.argv[argument + 1]

    output_paths = {
        "CG": (
            output_prefix + ".CG_homopolymer_indels.vcf",
            output_prefix + ".CG_homopolymer_indels.tsv",
        ),
        "AT": (
            output_prefix + ".AT_homopolymer_indels.vcf",
            output_prefix + ".AT_homopolymer_indels.tsv",
        ),
    }

    destinations = {
        group: open(paths[0], "w") for group, paths in output_paths.items()
    }
    tables = {
        group: open(paths[1], "w") for group, paths in output_paths.items()
    }
    non_hp_destination = open(
        output_prefix + ".non_homopolymer_indels.vcf", "w"
    )

    table_header = (
        "CHROM\tPOS\tREF\tALT\tTYPE\tEVENT_START_0BASED\t"
        "EVENT_END_0BASED\tHP_START_0BASED\tHP_END_0BASED\t"
        "HP_BASE\tHP_LENGTH\tHP_GROUP\n"
    )
    for table in tables.values():
        table.write(table_header)

    try:
        with open_text(input_vcf) as source:

            for raw in source:
                if raw.startswith("#"):
                    for destination in destinations.values():
                        destination.write(raw)
                    non_hp_destination.write(raw)
                    continue

                fields = raw.rstrip("\n").split("\t")
                chrom, pos_text, ref, alt_field = fields[0], fields[1], fields[3], fields[4]
                if chrom not in reference:
                    raise RuntimeError(f"Reference contig missing: {chrom}")

                fmt = fields[8].split(":")
                if "GT" not in fmt:
                    raise RuntimeError(f"Missing GT: {chrom}:{pos_text}")
                gt_index = fmt.index("GT")
                alts = alt_field.split(",")
                carried_indices = set()
                for sample_field in fields[9:]:
                    values = sample_field.split(":")
                    gt = values[gt_index] if gt_index < len(values) else "."
                    carried_indices.update(
                        int(x) for x in re.findall(r"[0-9]+", gt)
                        if 0 < int(x) <= len(alts)
                    )

                pos0 = int(pos_text) - 1
                sequence = reference[chrom]
                groups = set()
                hp_hit = False
                saw_indel = False
                for alt_index in sorted(carried_indices):
                    alt = alts[alt_index - 1]
                    if not re.fullmatch(r"[ACGTNacgtn]+", alt) or len(alt)==len(ref):
                        continue
                    left_trim, deleted, inserted = trim_shared(ref.upper(), alt.upper())
                    if not deleted and not inserted:
                        continue
                    saw_indel = True
                    event_start = pos0 + left_trim
                    if deleted:
                        event_end = event_start + len(deleted)
                        if sequence[event_start:event_end] != deleted:
                            raise RuntimeError(
                                f"Reference mismatch at {chrom}:{pos_text}: "
                                f"expected {deleted}, found {sequence[event_start:event_end]}"
                            )
                        hits = overlapping_interval_runs(chrom, event_start, event_end)
                    else:
                        event_end = event_start
                        hits = insertion_boundary_runs(chrom, event_start, inserted)

                    if deleted and inserted:
                        event_type = "COMPLEX"
                    elif deleted:
                        event_type = "DEL"
                    else:
                        event_type = "INS"
                    hp_hit = hp_hit or bool(hits)
                    groups.update(run[3] for run in hits)
                    for start, end, base, group in hits:
                        tables[group].write(
                            f"{chrom}\t{pos_text}\t{ref}\t{alt}\t"
                            f"{event_type}\t{event_start}\t{event_end}\t"
                            f"{start}\t{end}\t{base}\t{end-start}\t"
                            f"{group}\n"
                        )
                if not saw_indel:
                    continue
                for group in groups:
                    destinations[group].write(raw)
                if not hp_hit:
                    non_hp_destination.write(raw)
    finally:
        for handle in list(destinations.values()) + list(tables.values()):
            handle.close()
        non_hp_destination.close()
PY

###############################################################################
# 4. Compress and index the precisely classified homopolymer indels
###############################################################################

for method in syri anchorwave deepvariant; do
    for hp_group in CG AT; do
        bcftools view --threads "${THREADS}" -Oz -o "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf.gz" "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf"

        bcftools index -f \
            "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf.gz"

        rm "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf"
    done

    bcftools view --threads "${THREADS}" -Oz -o "${OUT}/${method}.non_homopolymer_indels.vcf.gz" "${OUT}/${method}.non_homopolymer_indels.vcf"

    bcftools index -f \
        "${OUT}/${method}.non_homopolymer_indels.vcf.gz"

    rm "${OUT}/${method}.non_homopolymer_indels.vcf"
done

# One-sample VCF containing only genotypes that carry an ALT allele.
subset_alt_sample() {
    local source_vcf=$1 sample_name=$2 target_vcf=$3
    bcftools view -s "${sample_name}" -Ov "${source_vcf}" |
    awk 'BEGIN { FS=OFS="\t" }
         /^#/ { print; next }
         {
             gt_col=0; n=split($9, fmt, ":")
             for (k=1; k<=n; k++) if (fmt[k]=="GT") gt_col=k
             if (!gt_col) { print "Missing GT: " $1 ":" $2 > "/dev/stderr"; exit 1 }
             split($10, f, ":"); a=split(f[gt_col], alleles, /[^0-9]+/)
             for (j=1; j<=a; j++)
                 if (alleles[j] ~ /^[0-9]+$/ && alleles[j]+0>0) { print; break }
         }' |
    bcftools view -Oz -o "${target_vcf}" -
}

###############################################################################
# 5. Compare each line separately against DeepVariant
#
# DeepVariant uses B73D7c while SyRI/AnchorWave use B73D7c1.
###############################################################################

while read -r assembly_sample dv_sample; do

    for method in syri anchorwave; do
        for hp_group in CG AT; do

        assembly_vcf="${OUT}/${method}.${hp_group}_homopolymer_indels.vcf.gz"
        dv_vcf="${OUT}/deepvariant.${hp_group}_homopolymer_indels.vcf.gz"

        assembly_sample_vcf="${OUT}/${method}.${assembly_sample}.${hp_group}_homopolymer_indels.vcf.gz"
        dv_sample_vcf="${OUT}/deepvariant.${assembly_sample}.${hp_group}_homopolymer_indels.vcf.gz"

        # Handle the B73D7c/B73D7c1 naming difference if either method
        # used the shorter name in the new merged inputs.
        assembly_name="${assembly_sample}"
        if [[ "${assembly_sample}" == "B73D7c1" ]] &&
           ! bcftools query -l "${assembly_vcf}" | grep -Fxq "B73D7c1"; then
            assembly_name="B73D7c"
        fi
        dv_name="${dv_sample}"
        if [[ "${dv_sample}" == "B73D7c" ]] &&
           ! bcftools query -l "${dv_vcf}" | grep -Fxq "B73D7c"; then
            dv_name="B73D7c1"
        fi

        subset_alt_sample "${assembly_vcf}" "${assembly_name}" "${assembly_sample_vcf}"
        subset_alt_sample "${dv_vcf}" "${dv_name}" "${dv_sample_vcf}"

        bcftools index -f "${assembly_sample_vcf}"
        bcftools index -f "${dv_sample_vcf}"

        # Assembly-method indels absent from DeepVariant.
        bcftools isec \
            -c none \
            -n=1 \
            -w1 \
            -Oz \
            -o "${OUT}/${method}.${assembly_sample}.${hp_group}_homopolymer_indels.not_in_deepvariant.vcf.gz" \
            "${assembly_sample_vcf}" \
            "${dv_sample_vcf}"

        # Exact normalized indels shared with DeepVariant.
        bcftools isec \
            -c none \
            -n=2 \
            -w1 \
            -Oz \
            -o "${OUT}/${method}.${assembly_sample}.${hp_group}_homopolymer_indels.shared_with_deepvariant.vcf.gz" \
            "${assembly_sample_vcf}" \
            "${dv_sample_vcf}"

        bcftools index -f \
            "${OUT}/${method}.${assembly_sample}.${hp_group}_homopolymer_indels.not_in_deepvariant.vcf.gz"

        bcftools index -f \
            "${OUT}/${method}.${assembly_sample}.${hp_group}_homopolymer_indels.shared_with_deepvariant.vcf.gz"

        done
    done

done <<'SAMPLES'
B73D1 B73D1
B73D2 B73D2
B73D4 B73D4
B73D7c1 B73D7c
B73D10 B73D10
SAMPLES

###############################################################################
# 6. Print result counts
###############################################################################

printf "METHOD\tSAMPLE\tHP_GROUP\tHP_INDELS_NOT_IN_DV\tHP_INDELS_SHARED_WITH_DV\n"

for method in syri anchorwave; do
    for sample in B73D1 B73D2 B73D4 B73D7c1 B73D10; do
        for hp_group in CG AT; do

        unique_vcf="${OUT}/${method}.${sample}.${hp_group}_homopolymer_indels.not_in_deepvariant.vcf.gz"
        shared_vcf="${OUT}/${method}.${sample}.${hp_group}_homopolymer_indels.shared_with_deepvariant.vcf.gz"

        unique_count=$(bcftools view -H "${unique_vcf}" | wc -l)
        shared_count=$(bcftools view -H "${shared_vcf}" | wc -l)

        printf "%s\t%s\t%s\t%s\t%s\n" \
            "${method}" "${sample}" "${hp_group}" \
            "${unique_count}" "${shared_count}"
        done
    done
done


###############################################################################
# 7. Per-sample, per-ALT indel-size counts. VCF rows stay intact.
#    Plotting input: combine barplot_counts.tsv across the four stage dirs.
###############################################################################
python3 - "${OUT}" "01_hp_core" "non_homopolymer_indels" <<'PLOT_PY'
import collections
import gzip
import pathlib
import re
import sys

out = pathlib.Path(sys.argv[1])
stage = sys.argv[2]
retained_suffix = sys.argv[3]
counts = collections.Counter()

def add_counts(method, stage_label, path):
    with gzip.open(path, "rt") as handle:
        names = []
        for line in handle:
            if line.startswith("#CHROM"):
                names = line.rstrip("\n").split("\t")[9:]
            elif line.startswith("#"):
                continue
            else:
                f = line.rstrip("\n").split("\t")
                if not names:
                    raise ValueError(f"Sample header missing from {path}")
                alts = f[4].split(",")
                fmt = f[8].split(":")
                if "GT" not in fmt:
                    raise ValueError(f"GT missing at {f[0]}:{f[1]} in {path}")
                index = fmt.index("GT")
                for name, sample_field in zip(names, f[9:]):
                    fields = sample_field.split(":")
                    gt = fields[index] if index < len(fields) else "."
                    canonical = "B73D7c" if name == "B73D7c1" else name
                    for allele_index in {int(a) for a in re.findall(r"[0-9]+", gt)
                                         if 0 < int(a) <= len(alts)}:
                        allele = alts[allele_index - 1]
                        if not re.fullmatch(r"[ACGTNacgtn]+", allele):
                            continue
                        size = len(allele)-len(f[3])
                        if size == 0:
                            continue
                        typ = "INS" if size > 0 else "DEL"
                        counts[(stage_label, method, canonical, typ, size)] += 1
for caller in ("syri", "anchorwave", "deepvariant"):
    if stage == "01_hp_core":
        add_counts(caller, "00_input", out / f"{caller}.indels.normalized.vcf.gz")
    add_counts(caller, stage, out / f"{caller}.{retained_suffix}.vcf.gz")

with (out / "barplot_counts.tsv").open("w") as target:
    target.write("stage\tcaller\tsample\tindel_type\tsigned_size_bp\tcount\n")
    for (label, caller, sample, typ, size), count in sorted(counts.items()):
        target.write(f"{label}\t{caller}\t{sample}\t{typ}\t{size}\t{count}\n")
print(f"Barplot data: {out / 'barplot_counts.tsv'}")
PLOT_PY
