#!/bin/bash
#SBATCH --job-name=getminmax
#SBATCH --mem=26G

input_tsv=$1
out_tsv=$2

read xmin xmax ymin ymax <<< $(zcat "$input_tsv" | \
  awk -F'\t' 'NR>1 {
    if (NR==2) { xmin=$1; xmax=$1; ymin=$2; ymax=$2 }
    if ($1 < xmin) xmin=$1
    if ($1 > xmax) xmax=$1
    if ($2 < ymin) ymin=$2
    if ($2 > ymax) ymax=$2
  }
  END { print xmin, xmax, ymin, ymax }')

printf "xmin\t%s\n" "$xmin" >  "$out_tsv"
printf "xmax\t%s\n" "$xmax" >> "$out_tsv"
printf "ymin\t%s\n" "$ymin" >> "$out_tsv"
printf "ymax\t%s\n" "$ymax" >> "$out_tsv"
