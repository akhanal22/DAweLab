#!/bin/bash
#SBATCH --job-name=04_hp_flank1_di_tri
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
HP_FLANK=1
MIN_DINUC_COPIES=3
MIN_TRINUC_COPIES=3

INPUT_ROOT=${1:-vcf_to_compare/three_callers_private_500kb}

SYRI="${INPUT_ROOT}/syri/syri.private.500kb.filtered.snps_indels.merged.vcf.gz"
ANCHORWAVE="${INPUT_ROOT}/anchorwave/anchorwave.private.500kb.filtered.snps_indels.merged.vcf.gz"
DEEPVARIANT="${INPUT_ROOT}/deepvariant/deepvariant.private.500kb.filtered.snps_indels.merged.vcf.gz"

OUT="${INPUT_ROOT}/repeat_filter_comparison/04_hp_flank1_di_tri"
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
# 3. Build reference homopolymer and perfect di-/trinucleotide BED files
#    and classify overlapping indels
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
# The conservative filtering mask expands every qualifying reference run by
# HP_FLANK bp on each side (clipped at chromosome boundaries). Thus, any-base
# insertion directly before, inside, or directly after the run is removed.
# Deletions and complex indels are removed when their changed reference
# interval overlaps either the core run or its flanking mask.
#
# A strict dinucleotide repeat has a minimal period of exactly 2 bp (the two
# motif bases differ) and at least MIN_DINUC_COPIES complete copies in the
# reference, for example ACACAC. Any indel with net size 1-50 bp is removed
# when its event overlaps one of these reference dinucleotide-repeat tracts.
#
# A strict trinucleotide repeat has a minimal period of exactly 3 bp and at
# least MIN_TRINUC_COPIES complete copies in the reference, for example
# AAGAAGAAG. Homopolymer motifs such as AAA are excluded. Any indel with net
# size 1-50 bp is removed when its event overlaps one of these reference
# trinucleotide-repeat tracts. Insertions at either repeat boundary count as
# touching the repeat and are removed.
###############################################################################

python3 - \
    "${REF}" "${MIN_HP}" "${HP_FLANK}" "${MIN_DINUC_COPIES}" \
    "${MIN_TRINUC_COPIES}" \
    "${OUT}/reference.homopolymers.min${MIN_HP}.bed" \
    "${OUT}/reference.CG_homopolymers.min${MIN_HP}.bed" \
    "${OUT}/reference.AT_homopolymers.min${MIN_HP}.bed" \
    "${OUT}/reference.homopolymers.min${MIN_HP}.plus${HP_FLANK}bp.bed" \
    "${OUT}/reference.CG_homopolymers.min${MIN_HP}.plus${HP_FLANK}bp.bed" \
    "${OUT}/reference.AT_homopolymers.min${MIN_HP}.plus${HP_FLANK}bp.bed" \
    "${OUT}/reference.perfect_dinucleotide_repeats.min${MIN_DINUC_COPIES}copies.bed" \
    "${OUT}/reference.perfect_trinucleotide_repeats.min${MIN_TRINUC_COPIES}copies.bed" \
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
flank = int(sys.argv[3])
minimum_dinucleotide_copies = int(sys.argv[4])
minimum_trinucleotide_copies = int(sys.argv[5])
all_bed = sys.argv[6]
cg_bed = sys.argv[7]
at_bed = sys.argv[8]
padded_all_bed = sys.argv[9]
padded_cg_bed = sys.argv[10]
padded_at_bed = sys.argv[11]
dinucleotide_bed = sys.argv[12]
trinucleotide_bed = sys.argv[13]

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
padded_starts_by_chrom = {}
pattern = re.compile(
    "|".join(f"{base}{{{minimum},}}" for base in "ACGT")
)

with open(all_bed, "w") as all_out, \
     open(cg_bed, "w") as cg_out, \
     open(at_bed, "w") as at_out, \
     open(padded_all_bed, "w") as padded_all_out, \
     open(padded_cg_bed, "w") as padded_cg_out, \
     open(padded_at_bed, "w") as padded_at_out:

    for chrom, sequence in reference.items():
        chromosome_runs = []

        for match in pattern.finditer(sequence):
            start = match.start()
            end = match.end()
            base = sequence[start]
            group = "CG" if base in "CG" else "AT"
            run = (start, end, base, group)
            chromosome_runs.append(run)
            padded_start = max(0, start - flank)
            padded_end = min(len(sequence), end + flank)

            bed_line = (
                f"{chrom}\t{start}\t{end}\t{base}\t{end - start}\t{group}\n"
            )
            all_out.write(bed_line)
            if group == "CG":
                cg_out.write(bed_line)
            else:
                at_out.write(bed_line)

            padded_bed_line = (
                f"{chrom}\t{padded_start}\t{padded_end}\t{base}\t"
                f"{end - start}\t{group}\t{start}\t{end}\n"
            )
            padded_all_out.write(padded_bed_line)
            if group == "CG":
                padded_cg_out.write(padded_bed_line)
            else:
                padded_at_out.write(padded_bed_line)

        runs_by_chrom[chrom] = chromosome_runs
        starts_by_chrom[chrom] = [run[0] for run in chromosome_runs]
        padded_starts_by_chrom[chrom] = [
            max(0, run[0] - flank) for run in chromosome_runs
        ]

def find_perfect_tandem_runs(sequence, motif_length, minimum_copies):
    """Find non-overlapping perfect repeats of minimal period 2 or 3."""
    runs = []
    minimum_length = motif_length * minimum_copies
    index = 0

    while index + minimum_length <= len(sequence):
        motif = sequence[index:index + motif_length]

        # A length-2 or length-3 motif has a smaller period only when all of
        # its bases are identical. Excluding those motifs prevents a
        # homopolymer from also being called a di-/trinucleotide repeat.
        if any(base not in "ACGT" for base in motif) or len(set(motif)) == 1:
            index += 1
            continue

        end = index + motif_length
        while (
            end < len(sequence)
            and sequence[end] in "ACGT"
            and sequence[end] == sequence[end - motif_length]
        ):
            end += 1

        # Retain only complete motif copies.
        complete_end = (
            index
            + motif_length * ((end - index) // motif_length)
        )
        if complete_end - index >= minimum_length:
            runs.append((index, complete_end, motif))
            index = complete_end
        else:
            index += 1

    return runs

dinucleotide_runs_by_chrom = {}
dinucleotide_starts_by_chrom = {}

with open(dinucleotide_bed, "w") as dinucleotide_out:
    for chrom, sequence in reference.items():
        chromosome_runs = find_perfect_tandem_runs(
            sequence, 2, minimum_dinucleotide_copies
        )
        dinucleotide_runs_by_chrom[chrom] = chromosome_runs
        dinucleotide_starts_by_chrom[chrom] = [
            run[0] for run in chromosome_runs
        ]

        for start, end, motif in chromosome_runs:
            dinucleotide_out.write(
                f"{chrom}\t{start}\t{end}\t{motif}\t"
                f"{(end - start) // 2}\n"
            )

trinucleotide_runs_by_chrom = {}
trinucleotide_starts_by_chrom = {}

with open(trinucleotide_bed, "w") as trinucleotide_out:
    for chrom, sequence in reference.items():
        chromosome_runs = find_perfect_tandem_runs(
            sequence, 3, minimum_trinucleotide_copies
        )
        trinucleotide_runs_by_chrom[chrom] = chromosome_runs
        trinucleotide_starts_by_chrom[chrom] = [
            run[0] for run in chromosome_runs
        ]

        for start, end, motif in chromosome_runs:
            trinucleotide_out.write(
                f"{chrom}\t{start}\t{end}\t{motif}\t"
                f"{(end - start) // 3}\n"
            )

def overlapping_interval_runs(chrom, start, end):
    """Return reference runs overlapping the 0-based half-open interval."""
    runs = runs_by_chrom.get(chrom, [])
    padded_starts = padded_starts_by_chrom.get(chrom, [])
    if not runs or start >= end:
        return []

    index = bisect.bisect_right(padded_starts, start)
    while index > 0 and runs[index - 1][1] + flank > start:
        index -= 1
    hits = []
    while index < len(runs) and padded_starts[index] < end:
        run = runs[index]
        padded_end = run[1] + flank
        if padded_end > start:
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

        # With a 1-bp flank, the two run-edge boundaries represent insertions
        # directly before or after the homopolymer. Any inserted base is
        # conservatively classified as homopolymer-associated.
        inside_or_adjacent = start <= boundary <= end

        if inside_or_adjacent:
            hits.append(run)

    return hits

def overlapping_tandem_runs(
    runs_by_chromosome, starts_by_chromosome, chrom, start, end
):
    """Return perfect tandem repeats overlapping a changed interval."""
    runs = runs_by_chromosome.get(chrom, [])
    starts = starts_by_chromosome.get(chrom, [])
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

def tandem_insertion_runs(
    runs_by_chromosome, starts_by_chromosome, chrom, boundary
):
    """Return tandem repeats containing or touching an insertion boundary."""
    runs = runs_by_chromosome.get(chrom, [])
    starts = starts_by_chromosome.get(chrom, [])
    if not runs:
        return []

    center = bisect.bisect_right(starts, boundary)
    first = max(0, center - 2)
    last = min(len(runs), center + 1)
    return [
        run for run in runs[first:last]
        if run[0] <= boundary <= run[1]
    ]

for argument in range(14, len(sys.argv), 2):
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
    non_hp_core_destination = open(
        output_prefix + ".non_homopolymer_core_indels.vcf", "w"
    )
    dinucleotide_destination = open(
        output_prefix + ".dinucleotide_repeat_indels.size1-50.vcf", "w"
    )
    dinucleotide_table = open(
        output_prefix + ".dinucleotide_repeat_indels.size1-50.tsv", "w"
    )
    trinucleotide_destination = open(
        output_prefix + ".trinucleotide_repeat_indels.size1-50.vcf", "w"
    )
    trinucleotide_table = open(
        output_prefix + ".trinucleotide_repeat_indels.size1-50.tsv", "w"
    )
    dinucleotide_clean_destination = open(
        output_prefix
        + ".non_homopolymer_non_dinucleotide_indels.vcf",
        "w",
    )
    combined_clean_destination = open(
        output_prefix
        + ".non_homopolymer_non_dinucleotide_non_trinucleotide_indels.vcf",
        "w",
    )

    table_header = (
        "CHROM\tPOS\tREF\tALT\tTYPE\tEVENT_START_0BASED\t"
        "EVENT_END_0BASED\tHP_START_0BASED\tHP_END_0BASED\t"
        "PADDED_START_0BASED\tPADDED_END_0BASED\tHP_BASE\tHP_LENGTH\t"
        "HP_GROUP\n"
    )
    for table in tables.values():
        table.write(table_header)

    dinucleotide_table.write(
        "CHROM\tPOS\tREF\tALT\tTYPE\tINDEL_SIZE\t"
        "EVENT_START_0BASED\tEVENT_END_0BASED\tTR_START_0BASED\t"
        "TR_END_0BASED\tMOTIF\tREPEAT_COPIES\n"
    )
    trinucleotide_table.write(
        "CHROM\tPOS\tREF\tALT\tTYPE\tINDEL_SIZE\t"
        "EVENT_START_0BASED\tEVENT_END_0BASED\tTR_START_0BASED\t"
        "TR_END_0BASED\tMOTIF\tREPEAT_COPIES\n"
    )

    try:
        with open_text(input_vcf) as source:

            for raw in source:
                if raw.startswith("#"):
                    for destination in destinations.values():
                        destination.write(raw)
                    non_hp_core_destination.write(raw)
                    non_hp_destination.write(raw)
                    dinucleotide_destination.write(raw)
                    trinucleotide_destination.write(raw)
                    dinucleotide_clean_destination.write(raw)
                    combined_clean_destination.write(raw)
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
                hp_core_hit = False
                hp_hit = False
                saw_indel = False
                di_hit = False
                tri_hit = False
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

                    # Compare the event against the original homopolymer run,
                    # before applying the 1-bp flank used by 'hits'.
                    for hp_start, hp_end, hp_base, _ in hits:
                        if deleted:
                            core_overlap = hp_start < event_end and hp_end > event_start
                        else:
                            core_overlap = hp_start < event_start < hp_end or (
                                event_start in (hp_start, hp_end)
                                and inserted
                                and all(base == hp_base for base in inserted)
                            )
                        hp_core_hit = hp_core_hit or core_overlap

                    if deleted and inserted:
                        event_type = "COMPLEX"
                    elif deleted:
                        event_type = "DEL"
                    else:
                        event_type = "INS"
                    hp_hit = hp_hit or bool(hits)
                    groups.update(run[3] for run in hits)
                    indel_size = abs(len(alt)-len(ref))
                    if 1 <= indel_size <= 50:
                        if deleted:
                            dinuc_hits = overlapping_tandem_runs(
                                dinucleotide_runs_by_chrom, dinucleotide_starts_by_chrom,
                                chrom, event_start, event_end
                            )
                            trinuc_hits = overlapping_tandem_runs(
                                trinucleotide_runs_by_chrom, trinucleotide_starts_by_chrom,
                                chrom, event_start, event_end
                            )
                        else:
                            dinuc_hits = tandem_insertion_runs(
                                dinucleotide_runs_by_chrom, dinucleotide_starts_by_chrom,
                                chrom, event_start
                            )
                            trinuc_hits = tandem_insertion_runs(
                                trinucleotide_runs_by_chrom, trinucleotide_starts_by_chrom,
                                chrom, event_start
                            )
                    else:
                        dinuc_hits = []
                        trinuc_hits = []
                    di_hit = di_hit or bool(dinuc_hits)
                    tri_hit = tri_hit or bool(trinuc_hits)
                    for start, end, motif in dinuc_hits:
                        dinucleotide_table.write(
                            f"{chrom}\t{pos_text}\t{ref}\t{alt}\t"
                            f"{event_type}\t{indel_size}\t{event_start}\t"
                            f"{event_end}\t{start}\t{end}\t{motif}\t"
                            f"{(end-start)//2}\n"
                        )
                    for start, end, motif in trinuc_hits:
                        trinucleotide_table.write(
                            f"{chrom}\t{pos_text}\t{ref}\t{alt}\t"
                            f"{event_type}\t{indel_size}\t{event_start}\t"
                            f"{event_end}\t{start}\t{end}\t{motif}\t"
                            f"{(end-start)//3}\n"
                        )
                    for start, end, base, group in hits:
                        tables[group].write(
                            f"{chrom}\t{pos_text}\t{ref}\t{alt}\t"
                            f"{event_type}\t{event_start}\t{event_end}\t"
                            f"{start}\t{end}\t{max(0,start-flank)}\t"
                            f"{min(len(sequence),end+flank)}\t{base}\t"
                            f"{end-start}\t{group}\n"
                        )
                if not saw_indel:
                    continue
                for group in groups:
                    destinations[group].write(raw)
                if not hp_core_hit:
                    non_hp_core_destination.write(raw)
                if not hp_hit:
                    non_hp_destination.write(raw)
                if di_hit:
                    dinucleotide_destination.write(raw)
                if tri_hit:
                    trinucleotide_destination.write(raw)
                if not hp_hit and not di_hit:
                    dinucleotide_clean_destination.write(raw)
                if not hp_hit and not di_hit and not tri_hit:
                    combined_clean_destination.write(raw)
    finally:
        for handle in list(destinations.values()) + list(tables.values()):
            handle.close()
        non_hp_core_destination.close()
        non_hp_destination.close()
        dinucleotide_destination.close()
        dinucleotide_table.close()
        trinucleotide_destination.close()
        trinucleotide_table.close()
        dinucleotide_clean_destination.close()
        combined_clean_destination.close()
PY

###############################################################################
# 4. Compress and index classified and filtered indel VCFs
###############################################################################

for method in syri anchorwave deepvariant; do
    for hp_group in CG AT; do
        bcftools view --threads "${THREADS}" -Oz -o "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf.gz" "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf"

        bcftools index -f \
            "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf.gz"

        rm "${OUT}/${method}.${hp_group}_homopolymer_indels.vcf"
    done

    bcftools view --threads "${THREADS}" -Oz -o "${OUT}/${method}.non_homopolymer_core_indels.vcf.gz" "${OUT}/${method}.non_homopolymer_core_indels.vcf"
    bcftools index -f "${OUT}/${method}.non_homopolymer_core_indels.vcf.gz"
    rm "${OUT}/${method}.non_homopolymer_core_indels.vcf"

    bcftools view --threads "${THREADS}" -Oz -o "${OUT}/${method}.non_homopolymer_indels.vcf.gz" "${OUT}/${method}.non_homopolymer_indels.vcf"

    bcftools index -f \
        "${OUT}/${method}.non_homopolymer_indels.vcf.gz"

    rm "${OUT}/${method}.non_homopolymer_indels.vcf"

    for category in \
        dinucleotide_repeat_indels.size1-50 \
        trinucleotide_repeat_indels.size1-50 \
        non_homopolymer_non_dinucleotide_indels \
        non_homopolymer_non_dinucleotide_non_trinucleotide_indels; do

        bcftools view --threads "${THREADS}" -Oz -o "${OUT}/${method}.${category}.vcf.gz" "${OUT}/${method}.${category}.vcf"

        bcftools index -f \
            "${OUT}/${method}.${category}.vcf.gz"

        rm "${OUT}/${method}.${category}.vcf"
    done
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
python3 - "${OUT}" "04_hp_flank1_di_tri" "non_homopolymer_non_dinucleotide_non_trinucleotide_indels" <<'PLOT_PY'
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
    add_counts(caller, "00_input", out / f"{caller}.indels.normalized.vcf.gz")
    add_counts(caller, "01_hp_core", out / f"{caller}.non_homopolymer_core_indels.vcf.gz")
    add_counts(caller, "02_hp_flank1", out / f"{caller}.non_homopolymer_indels.vcf.gz")
    add_counts(caller, "03_hp_flank1_di", out / f"{caller}.non_homopolymer_non_dinucleotide_indels.vcf.gz")
    add_counts(caller, stage, out / f"{caller}.{retained_suffix}.vcf.gz")

with (out / "barplot_counts.tsv").open("w") as target:
    target.write("stage\tcaller\tsample\tindel_type\tsigned_size_bp\tcount\n")
    for (label, caller, sample, typ, size), count in sorted(counts.items()):
        target.write(f"{label}\t{caller}\t{sample}\t{typ}\t{size}\t{count}\n")
print(f"Barplot data: {out / 'barplot_counts.tsv'}")
PLOT_PY

# C/G homopolymer + 1-bp flank retention per sample/caller (allele counts).
python3 - "${OUT}" <<'CG_COUNTS_PY'
#!/usr/bin/env python3
"""Count indel ALTs overlapping reference C/G homopolymers plus their 1-bp flanks.

Run automatically by 04-homopolymer_flanks_di_tri_repeat.sh.

Each CHROM,POS,REF,ALT allele is counted once per sample, even when its intact
VCF record contains other ALT alleles. No VCF is modified or split.
"""

import csv
import gzip
import re
import sys
from pathlib import Path


STAGE_DIR = Path(
    sys.argv[1]
    if len(sys.argv) > 1
    else "vcf_to_compare/three_callers_private_500kb/repeat_filter_comparison/04_hp_flank1_di_tri"
)
if len(sys.argv) > 2:
    raise SystemExit("Expected output directory from the current 04 script")


def cg_alleles(path):
    """CG TSV has one row per overlapping reference run; deduplicate hits."""
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"CHROM", "POS", "REF", "ALT"}
        if not rows.fieldnames or not required.issubset(rows.fieldnames):
            raise ValueError(f"Expected {required} in {path}")
        return {
            (row["CHROM"], row["POS"], row["REF"], row["ALT"])
            for row in rows
        }


def called_indels(path):
    """Sample -> set of carried sequence indels in a multisample VCF."""
    calls = {}
    with gzip.open(path, "rt") as source:
        names = None
        for line in source:
            if line.startswith("#CHROM"):
                names = line.rstrip("\n").split("\t")[9:]
                calls = {name: set() for name in names}
                continue
            if line.startswith("#"):
                continue
            if names is None:
                raise ValueError(f"Missing sample header in {path}")
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref = fields[0], fields[1], fields[3]
            alts = fields[4].split(",")
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise ValueError(f"Missing GT at {chrom}:{pos} in {path}")
            gt_index = fmt.index("GT")
            for name, value in zip(names, fields[9:]):
                sample_fields = value.split(":")
                gt = sample_fields[gt_index] if gt_index < len(sample_fields) else "."
                indices = {int(a) for a in re.findall(r"[0-9]+", gt)}
                for index in indices:
                    if not 1 <= index <= len(alts):
                        continue
                    alt = alts[index - 1]
                    if re.fullmatch(r"[ACGTNacgtn]+", alt) and len(alt) != len(ref):
                        calls[name].add((chrom, pos, ref, alt))
    return calls


def add_row(rows, caller, sample, calls, cg_hits):
    removed = calls & cg_hits
    total = len(calls)
    kept = total - len(removed)
    rows.append(
        (caller, sample, total, len(removed), kept, f"{100 * kept / total:.2f}" if total else "NA")
    )


rows = []
for caller in ("deepvariant", "anchorwave", "syri"):
    normalized = STAGE_DIR / f"{caller}.indels.normalized.vcf.gz"
    cg_table = STAGE_DIR / f"{caller}.CG_homopolymer_indels.tsv"
    for path in (normalized, cg_table):
        if not path.is_file():
            raise SystemExit(f"Missing input: {path}. Run the current 04 script first.")
    cg_hits = cg_alleles(cg_table)
    sample_calls = called_indels(normalized)
    all_calls = set()
    for original_name, calls in sample_calls.items():
        sample = "B73D7c" if original_name == "B73D7c1" else original_name
        add_row(rows, caller, sample, calls, cg_hits)
        all_calls.update(calls)
    add_row(rows, caller, "ALL_SAMPLES_UNIQUE", all_calls, cg_hits)

output = STAGE_DIR / "CG_plus_1bp_flank_retention_summary.tsv"
with output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        ["caller", "sample", "total_called_indel_alleles", "excluded_CG_plus_1bp_flank",
         "retained_after_CG_plus_1bp_flank", "percent_retained"]
    )
    writer.writerows(rows)

print(f"CG-plus-flank indel retention: {output}")
CG_COUNTS_PY

# Read the five stage totals from the per-size table and write one wide table
# for direct review. The input here is the VCF after the 500-kb filter.
python3 - "${OUT}" <<'STAGE_SUMMARY_PY'
import collections
import csv
import gzip
import pathlib
import sys

out = pathlib.Path(sys.argv[1])
stages = (
    "00_input",
    "01_hp_core",
    "02_hp_flank1",
    "03_hp_flank1_di",
    "04_hp_flank1_di_tri",
)
totals = collections.Counter()
with (out / "barplot_counts.tsv").open(newline="") as source:
    for row in csv.DictReader(source, delimiter="\t"):
        totals[(row["caller"], row["sample"], row["stage"])] += int(row["count"])

header = (
    "caller", "sample", "input_indels_after_500kb", "after_homopolymer_core",
    "after_homopolymer_plus_1bp_flank", "after_dinucleotide_1_50bp",
    "after_trinucleotide_1_50bp", "removed_by_homopolymer_core",
    "removed_by_flank", "removed_by_dinucleotide", "removed_by_trinucleotide",
    "percent_retained_final",
)
output = out / "indel_retention_by_stage.tsv"
with output.open("w", newline="") as target:
    writer = csv.writer(target, delimiter="\t")
    writer.writerow(header)
    for caller in ("deepvariant", "anchorwave", "syri"):
        with gzip.open(out / f"{caller}.indels.normalized.vcf.gz", "rt") as source:
            samples = []
            for line in source:
                if line.startswith("#CHROM"):
                    samples = ["B73D7c" if s == "B73D7c1" else s
                               for s in line.rstrip("\n").split("\t")[9:]]
                    break
        if not samples:
            raise ValueError(f"No sample columns in {caller} normalized VCF")
        all_counts = [0] * len(stages)
        for sample in samples:
            values = [totals[(caller, sample, stage)] for stage in stages]
            if any(after > before for before, after in zip(values, values[1:])):
                raise ValueError(f"Indel count increased between stages: {caller} {sample}")
            all_counts = [a + b for a, b in zip(all_counts, values)]
            removed = [before - after for before, after in zip(values, values[1:])]
            percent = 100 * values[-1] / values[0] if values[0] else 0.0
            writer.writerow([caller, sample, *values, *removed, f"{percent:.2f}"])
        removed = [before - after for before, after in zip(all_counts, all_counts[1:])]
        percent = 100 * all_counts[-1] / all_counts[0] if all_counts[0] else 0.0
        writer.writerow([caller, "ALL_SAMPLES_SUM", *all_counts, *removed, f"{percent:.2f}"])

print(f"Stage-by-stage retention: {output}")
STAGE_SUMMARY_PY
