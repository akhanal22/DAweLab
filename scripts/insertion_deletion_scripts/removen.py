python3 - <<'PY'
import gzip
import re

input_vcf = "syri.private.500kb.filtered.snps_indels.merged_againstb73all.vcf.gz"
output_vcf = "syri.private.500kb.filtered.snps_indels.merged_againstb73all.no_N3.vcf.gz"

n_run = re.compile(r"N{3,}", re.IGNORECASE)
total = retained = removed = 0

with gzip.open(input_vcf, "rt") as source, \
     gzip.open(output_vcf, "wt") as output:

    for line in source:
        if line.startswith("#"):
            output.write(line)
            continue

        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        ref = fields[3]
        alts = fields[4].split(",")
        total += 1

        if n_run.search(ref) or any(n_run.search(alt) for alt in alts):
            removed += 1
            continue

        output.write(line)
        retained += 1

print(f"Input variants: {total:,}")
print(f"Removed containing NNN: {removed:,}")
print(f"Retained: {retained:,}")
print(f"Output: {output_vcf}")
PY
