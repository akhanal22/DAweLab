python - <<'PY'
from pathlib import Path

source = Path("ragtag.scaffold.fasta")
output = Path("J702_flip.fa")
target = "chr2_RagTag"
flip_length = 1010498

complement = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn"
)

# Verify the chromosome before creating output.
with source.open() as f:
    matches = sum(
        line.startswith(">") and line[1:].split()[0] == target
        for line in f
    )
if matches != 1:
    raise SystemExit(f"Expected exactly one {target}; found {matches}.")
if output.exists():
    raise SystemExit(f"{output} already exists; choose a new output filename.")

def write_flipped(lines, out):
    seq = "".join(line.strip() for line in lines)
    if len(seq) < flip_length:
        raise ValueError("Chromosome is shorter than the requested segment.")

    seq = seq[:flip_length].translate(complement)[::-1] + seq[flip_length:]

    # Preserve the original sequence line lengths.
    offset = 0
    for line in lines:
        width = len(line.strip())
        out.write(seq[offset:offset + width] + "\n")
        offset += width

with source.open() as inp, output.open("w") as out:
    collecting = False
    chr2_lines = []

    for line in inp:
        if line.startswith(">"):
            if collecting:
                write_flipped(chr2_lines, out)
                chr2_lines = []
            collecting = line[1:].split()[0] == target
            out.write(line)
        elif collecting:
            chr2_lines.append(line)
        else:
            out.write(line)

    if collecting:
        write_flipped(chr2_lines, out)

print(f"Created {output}: flipped {target}:1-{flip_length}")
PY

samtools faidx J702_flip.fa
