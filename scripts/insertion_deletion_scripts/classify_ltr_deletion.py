#!/usr/bin/env python3

import sys

deletion_file = sys.argv[1]
ltr_file = sys.argv[2]

# Allow breakpoint differences up to 20 bp
tolerance = 20


# ---------------------------------------------------
# Helper functions
# ---------------------------------------------------

def close(a, b):
    return abs(a - b) <= tolerance


def overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1))


# ---------------------------------------------------
# Read deletion BED
#
# Expected format:
# chr  start  end  J702_chr  deletion_size  strand
# ---------------------------------------------------

deletions = []

with open(deletion_file) as f:
    for line in f:

        if not line.strip():
            continue

        x = line.rstrip().split("\t")

        # Skip header if present
        if x[0] in ("chr", "#chr"):
            continue

        if len(x) < 6:
            continue

        deletions.append({
            "chr": x[0],
            "start": int(x[1]),
            "end": int(x[2]),
            "j702_chr": x[3],
            "size": int(x[4]),
            "strand": x[5]
        })


# ---------------------------------------------------
# Read LTR structure table
# ---------------------------------------------------

ltrs = []

with open(ltr_file) as f:

    header = f.readline().rstrip().split("\t")

    # remove leading # from #chr
    header = [h.lstrip("#") for h in header]

    idx = {name: i for i, name in enumerate(header)}

    required = [
        "chr",
        "TE_start",
        "TE_end",
        "TE_name",
        "classification",
        "LTR_identity",
        "left_LTR_start",
        "left_LTR_end",
        "internal_start",
        "internal_end",
        "right_LTR_start",
        "right_LTR_end"
    ]

    for r in required:
        if r not in idx:
            sys.exit(f"ERROR: column {r} not found in LTR structure file")

    for line in f:

        if not line.strip():
            continue

        x = line.rstrip().split("\t")

        # GFF-derived coordinates are 1-based inclusive.
        # Convert starts to BED-style 0-based coordinates.
        ltrs.append({
            "chr": x[idx["chr"]],

            "TE_start": int(x[idx["TE_start"]]) - 1,
            "TE_end": int(x[idx["TE_end"]]),

            "TE_name": x[idx["TE_name"]],
            "class": x[idx["classification"]],
            "identity": x[idx["LTR_identity"]],

            "left_start": int(x[idx["left_LTR_start"]]) - 1,
            "left_end": int(x[idx["left_LTR_end"]]),

            "internal_start": int(x[idx["internal_start"]]) - 1,
            "internal_end": int(x[idx["internal_end"]]),

            "right_start": int(x[idx["right_LTR_start"]]) - 1,
            "right_end": int(x[idx["right_LTR_end"]])
        })


# ---------------------------------------------------
# Output header
# ---------------------------------------------------

print("\t".join([
    "chr",
    "deletion_start",
    "deletion_end",
    "deletion_size",
    "J702_chr",
    "deletion_strand",

    "TE_start",
    "TE_end",
    "TE_size",
    "TE_name",
    "TE_class",
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

    "classification",
    "start_diff",
    "end_diff",
    "expected_size",
    "size_diff"
]))


# ---------------------------------------------------
# Compare deletions against LTR structures
# ---------------------------------------------------

for d in deletions:

    for te in ltrs:

        if d["chr"] != te["chr"]:
            continue

        ov = overlap(
            d["start"],
            d["end"],
            te["TE_start"],
            te["TE_end"]
        )

        # Ignore deletions that do not overlap the TE at all
        if ov == 0:
            continue

        te_size = te["TE_end"] - te["TE_start"]

        left_size = te["left_end"] - te["left_start"]
        internal_size = te["internal_end"] - te["internal_start"]
        right_size = te["right_end"] - te["right_start"]

        category = "PARTIAL_OTHER"

        start_diff = "NA"
        end_diff = "NA"
        expected_size = "NA"
        size_diff = "NA"


        # ------------------------------------------------
        # Entire TE deleted
        # ------------------------------------------------

        if (
            close(d["start"], te["TE_start"]) and
            close(d["end"], te["TE_end"])
        ):

            category = "FULL_TE_DELETION"

            start_diff = abs(d["start"] - te["TE_start"])
            end_diff = abs(d["end"] - te["TE_end"])

            expected_size = te_size
            size_diff = abs(d["size"] - expected_size)


        # ------------------------------------------------
        # Left LTR only
        # ------------------------------------------------

        elif (
            close(d["start"], te["left_start"]) and
            close(d["end"], te["left_end"])
        ):

            category = "LEFT_LTR_ONLY"

            start_diff = abs(d["start"] - te["left_start"])
            end_diff = abs(d["end"] - te["left_end"])

            expected_size = left_size
            size_diff = abs(d["size"] - expected_size)


        # ------------------------------------------------
        # Right LTR only
        # ------------------------------------------------

        elif (
            close(d["start"], te["right_start"]) and
            close(d["end"], te["right_end"])
        ):

            category = "RIGHT_LTR_ONLY"

            start_diff = abs(d["start"] - te["right_start"])
            end_diff = abs(d["end"] - te["right_end"])

            expected_size = right_size
            size_diff = abs(d["size"] - expected_size)


        # ------------------------------------------------
        # Internal region only
        # ------------------------------------------------

        elif (
            close(d["start"], te["internal_start"]) and
            close(d["end"], te["internal_end"])
        ):

            category = "INTERNAL_ONLY"

            start_diff = abs(d["start"] - te["internal_start"])
            end_diff = abs(d["end"] - te["internal_end"])

            expected_size = internal_size
            size_diff = abs(d["size"] - expected_size)


        # ------------------------------------------------
        # Left LTR + internal deleted
        # Right LTR remains
        # ------------------------------------------------

        elif (
            close(d["start"], te["left_start"]) and
            close(d["end"], te["right_start"])
        ):

            category = "SOLO_LTR_CANDIDATE_RIGHT_REMAINS"

            start_diff = abs(d["start"] - te["left_start"])
            end_diff = abs(d["end"] - te["right_start"])

            expected_size = te["right_start"] - te["left_start"]
            size_diff = abs(d["size"] - expected_size)


        # ------------------------------------------------
        # Internal + right LTR deleted
        # Left LTR remains
        # ------------------------------------------------

        elif (
            close(d["start"], te["left_end"]) and
            close(d["end"], te["right_end"])
        ):

            category = "SOLO_LTR_CANDIDATE_LEFT_REMAINS"

            start_diff = abs(d["start"] - te["left_end"])
            end_diff = abs(d["end"] - te["right_end"])

            expected_size = te["right_end"] - te["left_end"]
            size_diff = abs(d["size"] - expected_size)


        print("\t".join(map(str, [
            d["chr"],
            d["start"],
            d["end"],
            d["size"],
            d["j702_chr"],
            d["strand"],

            te["TE_start"] + 1,
            te["TE_end"],
            te_size,
            te["TE_name"],
            te["class"],
            te["identity"],

            te["left_start"] + 1,
            te["left_end"],
            left_size,

            te["internal_start"] + 1,
            te["internal_end"],
            internal_size,

            te["right_start"] + 1,
            te["right_end"],
            right_size,

            category,
            start_diff,
            end_diff,
            expected_size,
            size_diff
        ])))
