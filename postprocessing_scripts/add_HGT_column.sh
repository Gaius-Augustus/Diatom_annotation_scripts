#!/bin/bash
# script: add_HGT_column.sh
# usage: ./add_HGT_column.sh input.tsv output.tsv

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

input_file="$1"
output_file="$2"

awk 'BEGIN { FS="\t"; OFS="\t" }
{
    if (NF == 1) {
        print $1, "HGT"
    } else {
        print $1, "HGT", $2
    }
}' "$input_file" > "$output_file"

echo "File with newky added HGT saved into the $output_file"

