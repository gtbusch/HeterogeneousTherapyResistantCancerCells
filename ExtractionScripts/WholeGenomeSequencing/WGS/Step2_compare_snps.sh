#!/bin/bash

# Get list of all exonic SNP files
FILES=(*_exonic_snps.txt)

echo "Found ${#FILES[@]} files: ${FILES[@]}"

# Much faster awk-based approach
awk '
BEGIN {
    OFS="\t"
}

FNR==1 {
    filenum++
    split(FILENAME, parts, "_")
    clone[filenum] = parts[1]  # Adjust if your naming is different
}
{
    variant = $1":"$2":"$4":"$5
    presence[variant, filenum] = 1
    if (!(variant in seen)) {
        variants[++varcount] = variant
        seen[variant] = 1
    }
}
END {
    printf "Variant"
    for (i=1; i<=filenum; i++) {
        printf "\t%s", clone[i]
    }
    printf "\n"
    
    for (v=1; v<=varcount; v++) {
        variant = variants[v]
        printf "%s", variant
        for (i=1; i<=filenum; i++) {
            if (presence[variant, i]) {
                printf "\t1"
            } else {
                printf "\t0"
            }
        }
        printf "\n"
    }
}
' *_exonic_snps.txt > comparison_table.txt

echo "Done! Results in comparison_table.txt"

