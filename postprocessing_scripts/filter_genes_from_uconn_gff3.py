#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description='Filter genes from a gff3 file based on a list of gene names')
parser.add_argument('-g', '--gff3', help='GFF3 file to filter', required=True)
parser.add_argument('-l', '--list', help='List of gene names to keep', required=True)
parser.add_argument('-o', '--output', help='Output file name', required=True)
args = parser.parse_args()

# open list file, chop trailing newline, store as key in dict (tx_dict), then split by dot and store first value in gene_dict
gene_dict = {}
tx_dict = {}
try:
    with open(args.list, 'r') as list_file:
        for line in list_file:
            tx_dict[line.rstrip()] = 1
            gene_dict[line.split('.')[0]] = 1
except IOError:
    print('Cannot open', args.list)

# open output file with try/except
try:
    with open(args.output, 'w') as output_file:
        try:
            with open(args.gff3, 'r') as gff3_file:
                for line in gff3_file:
                    if line.startswith('#'):
                        output_file.write(line)
                        continue
                    line = line.rstrip()
                    fields = line.split('\t')
                    if fields[2] == 'gene':
                        # last column has this format: ID=g1;Dbxref=GeneID:7204623;gene_biotype=protein_coding
                        # extract the value between ID= and ;
                        # Last column has format: ID=g1;Dbxref=GeneID:7204623;gene_biotype=protein_coding
                        attributes = fields[8]  # Last column is fields[8]
                        # Extract ID between ID= and ;
                        gene_id = None
                        for attribute in attributes.split(';'):
                            if attribute.startswith('ID='):
                                gene_id = attribute.split('=')[1]
                                break
                        if gene_id in gene_dict:
                            output_file.write(line + '\n')
                    elif fields[2] == 'mRNA':
                        # last column has this format: ID=g1.t1;Parent=g1;Dbxref=GeneID:7204623,PFAM:PF04055;product=hypothetical protein;note=EggNOG:Radical SAM superfamily
                        # extract the value between ID= and ;
                        attributes = fields[8]  # Last column is fields[8]
                        # Extract ID between ID= and ;
                        tx_id = None
                        for attribute in attributes.split(';'):
                            if attribute.startswith('ID='):
                                tx_id = attribute.split('=')[1]
                                break
                        if gene_id in tx_dict:
                            output_file.write(line + '\n')
                    elif fields[2] == 'exon' or fields[2] == 'CDS' or fields[2] == 'start_codon' or fields[2] == 'stop_codon' or fields[2] == 'intron' or fields[2] == 'five_prime_UTR' or fields[2] == 'three_prime_UTR':
                        # last colun has this format: ID=agat-cds-1;Parent=g1.t1;gene_id=g1;transcript_id=g1.t1;Dbxref=GeneID:7204623,PFAM:PF04055;product=hypothetical protein;note=EggNOG:Radical SAM superfamily
                        # extract the value between Parent= and ;
                        attributes = fields[8]
                        tx_id = None
                        for attribute in attributes.split(';'):
                            if attribute.startswith('Parent='):
                                tx_id = attribute.split('=')[1]
                                break
                        if tx_id in tx_dict:
                            output_file.write(line + '\n')
                    else:
                        print("Skipping line: ", line)
        except IOError:
            print('Cannot open', args.gff3)
except IOError:
    print('Cannot open', args.output)
