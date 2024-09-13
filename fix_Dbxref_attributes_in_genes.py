#!/usr/bin/env python3

import re
import argparse
from collections import defaultdict

def parse_attributes(attributes_str):
    """Parse the attributes column of a GFF3 line into a dictionary."""
    attributes = {}
    for attribute in attributes_str.split(";"):
        if "=" in attribute:
            key, value = attribute.split("=")
            attributes[key] = value
    return attributes

def attributes_to_str(attributes):
    """Convert the attribute dictionary back to a GFF3 string."""
    return ";".join([f"{key}={value}" for key, value in attributes.items()])

def process_gff3(input_file, output_file):
    gene_id_to_geneid = defaultdict(list)  # Dict to store gene_id -> list of GeneID:...
    
    # Step 1 (a): Parse the input file to collect GeneIDs from CDS lines
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("#"):
                continue  # Skip comments

            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue  # Skip malformed lines

            feature_type = fields[2]
            if feature_type == "CDS":
                attributes = parse_attributes(fields[8])
                
                # Extract Dbxref and gene_id
                if "Dbxref" in attributes and "gene_id" in attributes:
                    dbxrefs = attributes["Dbxref"].split(",")
                    gene_id = attributes["gene_id"]
                    
                    # Filter only GeneID entries and add to dict
                    for dbxref in dbxrefs:
                        if dbxref.startswith("GeneID:"):
                            if dbxref not in gene_id_to_geneid[gene_id]:
                                gene_id_to_geneid[gene_id].append(dbxref)

    # Step 2 (b): Ensure lists in the dict are unique and convert to comma-separated strings
    for gene_id, geneid_list in gene_id_to_geneid.items():
        gene_id_to_geneid[gene_id] = ",".join(sorted(set(geneid_list)))  # Remove duplicates and sort

    # Step 3 (c): Re-read the input file and modify gene feature lines
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)  # Write comments as is
                continue

            fields = line.strip().split('\t')
            if len(fields) != 9:
                outfile.write(line)  # Write malformed lines as is
                continue

            feature_type = fields[2]
            attributes = parse_attributes(fields[8])
            
            if feature_type == "gene":
                gene_id = attributes.get("ID")
                
                # If the gene_id has corresponding GeneID entries
                if gene_id in gene_id_to_geneid and gene_id_to_geneid[gene_id].strip(","):
                    dbxref_value = gene_id_to_geneid[gene_id]
                    
                    if "Dbxref" in attributes:
                        # Update existing Dbxref
                        attributes["Dbxref"] = dbxref_value
                    else:
                        # Add new Dbxref
                        attributes["Dbxref"] = dbxref_value
                    
                    # Write the updated gene feature line
                    fields[8] = attributes_to_str(attributes)
                    outfile.write("\t".join(fields) + "\n")
                else:
                    # Write the gene feature line as is if no GeneID found
                    outfile.write(line)
            else:
                # Write non-gene lines as is
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Fix the Dbxref= field in GFF3 gene models using CDS feature data.")
    parser.add_argument('-i', '--input', required=True, help='Input GFF3 file')
    parser.add_argument('-o', '--output', required=True, help='Output GFF3 file')

    args = parser.parse_args()

    process_gff3(args.input, args.output)

if __name__ == "__main__":
    main()
