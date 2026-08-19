#!/usr/bin/env python3

import sys
import re

gff = sys.argv[1]


def parse_attrs(s):
    d = {}

    for item in s.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v

    return d


repeat_regions = {}
children = {}


# --------------------------------------------------
# Read GFF
# --------------------------------------------------

with open(gff) as f:

    for line in f:

        if not line.strip() or line.startswith("#"):
            continue

        p = line.rstrip().split("\t")

        if len(p) < 9:
            continue

        chrom = p[0]
        feature = p[2]
        start = int(p[3])
        end = int(p[4])
        strand = p[6]

        attrs = parse_attrs(p[8])

        # Keep chr1-chr10
        m = re.match(r"^(chr(?:[1-9]|10))_RagTag.*$", chrom)

        if not m:
            continue

        new_chr = m.group(1)

        fid = attrs.get("ID")
        parent = attrs.get("Parent")

        record = {
            "chr": new_chr,
            "feature": feature,
            "start": start,
            "end": end,
            "strand": strand,
            "attrs": attrs
        }

        # Main intact LTR element
        if feature == "repeat_region":

            classification = attrs.get("classification", "")

            if classification.startswith("LTR/") and fid:
                repeat_regions[fid] = record

        # Child features
        if parent:
            children.setdefault(parent, []).append(record)


# --------------------------------------------------
# Header
# --------------------------------------------------

header = [
    "#chr",
    "TE_start",
    "TE_end",
    "TE_size",
    "TE_name",
    "classification",
    "strand",
    "LTR_identity",

    "left_LTR_start",
    "left_LTR_end",
    "left_LTR_size",

    "internal_start",
    "internal_end",
    "internal_size",

    "right_LTR_start",
    "right_LTR_end",
    "right_LTR_size",

    "left_TSD_start",
    "left_TSD_end",
    "left_TSD_seq",

    "right_TSD_start",
    "right_TSD_end",
    "right_TSD_seq"
]

print("\t".join(header))


# --------------------------------------------------
# Process each intact LTR
# --------------------------------------------------

for rid, te in repeat_regions.items():

    kids = children.get(rid, [])

    ltrs = []
    tsds = []

    for x in kids:

        if x["feature"] == "long_terminal_repeat":
            ltrs.append(x)

        elif x["feature"] == "target_site_duplication":
            tsds.append(x)

    ltrs.sort(key=lambda x: x["start"])
    tsds.sort(key=lambda x: x["start"])

    # Require two LTRs
    if len(ltrs) < 2:
        continue

    left_ltr = ltrs[0]
    right_ltr = ltrs[-1]

    # Internal region
    internal_start = left_ltr["end"] + 1
    internal_end = right_ltr["start"] - 1

    if internal_end >= internal_start:
        internal_size = internal_end - internal_start + 1
    else:
        internal_size = 0

    te_size = te["end"] - te["start"] + 1

    left_ltr_size = (
        left_ltr["end"] - left_ltr["start"] + 1
    )

    right_ltr_size = (
        right_ltr["end"] - right_ltr["start"] + 1
    )

    # TSD information
    tsd_seq = te["attrs"].get("tsd", "NA")

    if len(tsds) >= 2:

        left_tsd = tsds[0]
        right_tsd = tsds[-1]

        left_tsd_start = left_tsd["start"]
        left_tsd_end = left_tsd["end"]

        right_tsd_start = right_tsd["start"]
        right_tsd_end = right_tsd["end"]

    else:

        left_tsd_start = "NA"
        left_tsd_end = "NA"

        right_tsd_start = "NA"
        right_tsd_end = "NA"

    row = [
        te["chr"],
        te["start"],
        te["end"],
        te_size,
        te["attrs"].get("Name", rid),
        te["attrs"].get("classification", "NA"),
        te["strand"],
        te["attrs"].get("ltr_identity", "NA"),

        left_ltr["start"],
        left_ltr["end"],
        left_ltr_size,

        internal_start,
        internal_end,
        internal_size,

        right_ltr["start"],
        right_ltr["end"],
        right_ltr_size,

        left_tsd_start,
        left_tsd_end,
        tsd_seq,

        right_tsd_start,
        right_tsd_end,
        tsd_seq
    ]

    print("\t".join(map(str, row)))
