#!/usr/bin/env python3

import pandas as pd
import os

# Read the TSV file into a DataFrame
df = pd.read_csv('Orthogroups.tsv', sep='\t')

# Create an output directory if it doesn't exist
output_dir = 'output_files'
os.makedirs(output_dir, exist_ok=True)

# Iterate over each column except the first one (Orthogroup)
for column in df.columns[1:]:
    # Prepare the output filename, replacing spaces with underscores
    species_name = column.replace(' ', '_')
    output_filename = f"{output_dir}/{species_name}_assigned.lst"
    
    # Get the non-empty transcript IDs in the current column
    non_empty_transcripts = df[column].dropna()
    
    # Flatten the list of transcripts if they are comma-separated
    all_transcripts = []
    for transcripts in non_empty_transcripts:
        all_transcripts.extend(transcripts.split(','))

    # Write to the output file, one transcript ID per line
    with open(output_filename, 'w') as outfile:
        for transcript in all_transcripts:
            outfile.write(f"{transcript.strip()}\n")

print("Files generated successfully.")
