#!/bin/bash
## from tx csv prepare sample_roi.tsv.gz

input_csv=$1
out_tsv_gz=$2
# X Y gene cell_id Count
# using bash
cat $input_csv | awk -F, 'NR==1 {OFS="\t"; print "X", "Y", "gene", "cell_id", "transcript_id", "Count"} NR>1 {OFS="\t"; print $5, $6, $4, $2, $1, 1}' | \
               (head -n 1 && tail -n +2 | sort -k1,1n) | \
               gzip > $out_tsv_gz