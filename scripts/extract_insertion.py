#!/usr/bin/env python3

import re
import sys

paf = sys.argv[1]
min_ins = int(sys.argv[2]) if len(sys.argv) > 2 else 1000

print("\t".join([
    "chr",
    "start",
    "end",
    "insertion_size",
    "j702_chr",
    "j702_position",
    "strand"
]))

with open(paf) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip().split("\t")

        qname = fields[0]
        qstart = int(fields[2])
        strand = fields[4]

        tname = fields[5]
        tstart = int(fields[7])

        # Use primary alignments only
        tags = fields[12:]

        tp = None
        cigar = None

        for tag in tags:
            if tag.startswith("tp:A:"):
                tp = tag.split(":")[-1]
            elif tag.startswith("cg:Z:"):
                cigar = tag[5:]

        if tp != "P":
            continue

        if cigar is None:
            continue

        qpos = qstart
        tpos = tstart

        for length, op in re.findall(r"(\d+)([=XIDM])", cigar):
            length = int(length)

            if op in ("=", "X", "M"):
                qpos += length
                tpos += length

            elif op == "I":
                # Sequence present in B73D1 but absent from J702
                if length >= min_ins:
                    # BED-like 0-based coordinates
                    print("\t".join(map(str, [
                        qname,
                        qpos,
                        qpos + length,
                        length,
                        tname,
                        tpos,
                        strand
                    ])))

                qpos += length

            elif op == "D":
                # Sequence present in J702 but absent from B73D1
                tpos += length
