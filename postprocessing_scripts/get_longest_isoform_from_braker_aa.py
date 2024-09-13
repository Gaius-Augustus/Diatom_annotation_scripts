#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def get_longest_isoforms(fasta_file, output_file):
    # Dictionary to store the longest sequence for each gene
    gene_dict = {}

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract gene name (everything before the first dot)
        gene_name = record.id.split('.')[0]

        # Check if this gene is already in the dictionary
        if gene_name in gene_dict:
            # If the current sequence is longer, replace the existing one
            if len(record.seq) > len(gene_dict[gene_name].seq):
                gene_dict[gene_name] = record
        else:
            # If the gene is not in the dictionary, add it
            gene_dict[gene_name] = record

    # Write the longest isoforms to the output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(gene_dict.values(), output_handle, "fasta")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract longest isoforms from protein FASTA.")
    
    # Define arguments
    parser.add_argument("-i", "--input", required=True, help="Input protein FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with longest isoforms")
    
    # Parse arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    get_longest_isoforms(args.input, args.output)

if __name__ == "__main__":
    main()
