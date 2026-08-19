BEGIN {
    FS="\t"
    OFS="\t"
}

# Skip comments
/^#/ {
    next
}

{
    chr=$1

    # Keep chr1-chr10 and rename
    # chr1_RagTag_R -> chr1
    # chr1_RagTag_  -> chr1
    if (chr ~ /^chr[0-9]+_RagTag/) {
        sub(/_RagTag.*$/, "", chr)
    } else {
        next
    }

    if (chr !~ /^chr([1-9]|10)$/)
        next

    feature=$3

    class="NA"
    name="NA"
    tsd="NA"
    identity="NA"

    n=split($9,a,";")

    for(i=1;i<=n;i++) {

        x=a[i]

        if(x ~ /^classification=/) {
            sub(/^classification=/,"",x)
            class=x
        }

        if(x ~ /^Name=/) {
            sub(/^Name=/,"",x)
            name=x
        }

        if(x ~ /^tsd=/) {
            sub(/^tsd=/,"",x)
            tsd=x
        }

        if(x ~ /^ltr_identity=/) {
            sub(/^ltr_identity=/,"",x)
            identity=x
        }
    }

    # LTRs:
    # keep ONLY repeat_region so an intact LTR
    # is represented once
    if (class ~ /^LTR\//) {

        if (feature=="repeat_region") {
            print chr,$4-1,$5,name,class,tsd,identity,$7,feature
        }

        next
    }

    # Remove subfeatures
    if (feature=="target_site_duplication" ||
        feature=="long_terminal_repeat")
        next

    # Keep non-LTR transposons, Helitrons,
    # and relevant repeat_region records
    if (feature ~ /transposon/ ||
        feature ~ /Helitron/ ||
        feature=="repeat_region") {

        print chr,$4-1,$5,name,class,tsd,identity,$7,feature
    }
}
