import re
import sys

paf = sys.argv[1]
min_del = int(sys.argv[2]) if len(sys.argv) > 2 else 100

with open(paf) as f:
    for line in f:
        if line.startswith("#"):
            continue

        x = line.rstrip().split("\t")

        qname = x[0]
        strand = x[4]

        tname = x[5]
        tstart = int(x[7])

        cigar = None

        for field in x[12:]:
            if field.startswith("cg:Z:"):
                cigar = field[5:]
                break

        if cigar is None:
            continue

        # Current position on J702/reference
        ref_pos = tstart

        operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

        for length, op in operations:

            length = int(length)

            # consumes reference
            if op in ["M", "=", "X", "N"]:
                ref_pos += length

            # sequence present in reference/J702 but absent in query
            elif op == "D":

                if length >= min_del:

                    # BED coordinates: 0-based, half-open
                    start = ref_pos
                    end = ref_pos + length

                    print(
                        tname,
                        start,
                        end,
                        qname,
                        length,
                        strand,
                        sep="\t"
                    )

                ref_pos += length

            # insertion relative to reference does NOT move reference
            elif op == "I":
                pass
