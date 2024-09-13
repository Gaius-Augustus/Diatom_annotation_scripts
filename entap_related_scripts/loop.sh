#!/bin/bash
#SBATCH --job-name=script
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G

source activate datasets
# Define the directories
dir1="/core/labs/Wegrzyn/diatom_genomes/fixed_gff/gff" # path to GFF files
dir2="/core/labs/Wegrzyn/diatom_genomes/fixed_gff/entap" # path entap_results.tsv, which should be renamed to represent the species

# Loop through each *.gff file in the first directory
for gff_file in "$dir1"/*.gff; do
  # Extract the basename without extension
  basename=$(basename "$gff_file" .gff)

  # Look for the corresponding *.tsv file in the second directory
  tsv_file="$dir2/$basename.tsv"

  # Check if both files exist
  if [ -f "$gff_file" ] && [ -f "$tsv_file" ]; then
    echo "Found matching files: $gff_file and $tsv_file"

    # Run the Python script with the gff and tsv file as arguments
    python ncbi_gff.py "$tsv_file" "$gff_file" pdb_pfam_mapping.txt "$basename"_ncbi.gff
  else
    echo "No matching .tsv file found for $gff_file"
  fi
done

