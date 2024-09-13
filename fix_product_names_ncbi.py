#!/usr/bin/env python3

import re
import argparse

def fix_product_name(product_name):
    # Fix uncharacterized protein
    if product_name.startswith("uncharacterized protein"):
        product_name = "putative protein"
    
    # if a product name ends in uncharacterised, remove that
    product_name = re.sub(r'uncharacterised$', '', product_name, flags=re.IGNORECASE)
    # remove trailing whitespace
    product_name = product_name.strip()

    # Fix hypothetical protein only if followed by additional text, keeping the rest of the text
    if product_name.startswith("hypothetical protein "):
        product_name = product_name.replace("hypothetical protein", "putative protein", 1)
    
    # Fix plural forms (simple heuristic)
    product_name = re.sub(r'proteins\b', 'protein', product_name)

    # Remove suspicious phrases
    product_name = re.sub(r'truncat(ed|ion)?', '', product_name)
    product_name = re.sub(r'Truncated ', 'Nonfunctional ', product_name)
    product_name = re.sub(r'Fragment', '', product_name, flags=re.IGNORECASE)
    product_name = re.sub(r'partial', '', product_name, flags=re.IGNORECASE)

    # Remove strings containing three or more consecutive digits (may be part of a larger string)
    product_name = re.sub(r'\S*\d{3,}\S*', '', product_name)

    product_name = product_name.replace('_', ' ')  # Replace underscores with spaces

    # Change homologs/paralogs to '-like protein'
    product_name = re.sub(r'\b(Homolog|paralog)\b', '-like protein', product_name)
    product_name = re.sub(r'\shomolog', '-like proteins', product_name)
    # Shorten long product names to 100 characters
    if len(product_name) > 100:
        # Find the last space before the 100th character
        last_space = product_name[:100].rfind(' ')
        # replace product_name with the first 100 characters up to the last space
        product_name = product_name[:last_space].rstrip()

    # Replace 'gene' with 'protein'
    product_name = re.sub(r'\bgene\b', 'protein', product_name)

    # Fix product names starting with "and"
    if re.match(r'\s?and\s', product_name):
        product_name = product_name[4:].strip()

    # Fix product names starting with "a"
    if re.match(r'\s?a\s', product_name):
        product_name = product_name[2:].strip()

    # Remove "FOG" from product names
    product_name = product_name.replace("FOG", "")

    # Replace "possibly", "residue", and "unnamed" with empty strings or replacements
    product_name = re.sub(r'\bpossibly\b', '', product_name, flags=re.IGNORECASE)
    product_name = re.sub(r'\bresidue\b', '', product_name, flags=re.IGNORECASE)
    product_name = re.sub(r'\bunnamed\b', 'putative protein', product_name, flags=re.IGNORECASE)

    # if product ends in family, replace with family member
    product_name = re.sub(r'family$', 'family member', product_name, flags=re.IGNORECASE)

    # residue Fragilaria_radians.dr # this should be ignored, it's valuable information
    # conserved protein Thalassiosira_delicatula.dr I think this should also be ignored

    # Convert all-caps product names to title case
    if product_name.isupper():
        product_name = product_name.title()

    # Collapse multiple spaces into a single space and strip leading/trailing whitespace
    product_name = re.sub(r'\s+', ' ', product_name).strip()

    return product_name.strip()

def process_gff3(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) == 9:
                attributes = fields[8]
                
                # Check for "product=" in the attributes column
                if 'product=' in attributes:
                    product_match = re.search(r'product=([^;]+)', attributes)
                    if product_match:
                        product_name = product_match.group(1)
                        
                        # Fix the product name
                        fixed_product_name = fix_product_name(product_name)

                        # Replace the old product name with the fixed one
                        attributes = attributes.replace(product_name, fixed_product_name)
                        fields[8] = attributes

                # Remove "fragment" from product descriptions or comments
                attributes = re.sub(r'fragment', '', attributes, flags=re.IGNORECASE)

                # in attributes, replace 'Homologous to' with 'Similar to'
                attributes = re.sub(r'Homologous to', 'Similar to', attributes)
                # replace Homology with an empty string
                attributes = re.sub(r'Homology', '', attributes)
                # fix product names that start with a
                attributes = re.sub(r'product=a ', 'product=', attributes)
                # replace Homologues of by Similar to
                attributes = re.sub(r'Homologues of', 'Similar to', attributes)
                # replace product=Uncharacterised protein family by product=Hypothetical protein
                attributes = re.sub(r'product=Uncharacterised protein family', 'product=Hypothetical protein', attributes)
                # if a product name ends in rRNA, extend by -modifying protein
                attributes = re.sub(r'rRNA;', 'rRNA-modifying protein;', attributes)
                # delete the phrase "and related protein" from attributes
                attributes = re.sub(r' and related protein', '', attributes)
                # replace "Conserved hypothetical protein" with "Conserved protein"
                attributes = re.sub(r'Conserved hypothetical protein', 'Conserved protein', attributes)
                # replace " and related enzymes" with empty string
                attributes = re.sub(r' and related enzymes', '', attributes)
                # replace "Shiftless antiviral inhibitor of ribosomal frameshifting protein homolog;note=EggNOG:UPF0515 protein C19orf66 homolog" by "Similar to shiftless antiviral inhibitor of ribosomal frameshifting protein"
                attributes = re.sub(r'Shiftless antiviral inhibitor of ribosomal frameshifting protein homolog;note=EggNOG:UPF0515 protein C19orf66 homolog', 'Similar to shiftless antiviral inhibitor of ribosomal frameshifting protein', attributes)
                fields[8] = attributes
                # Write the fixed line
                outfile.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Fix suspect product names in a GFF3 file.")
    parser.add_argument('-i', '--input', required=True, help='Input GFF3 file')
    parser.add_argument('-o', '--output', required=True, help='Output GFF3 file')

    args = parser.parse_args()

    process_gff3(args.input, args.output)

if __name__ == "__main__":
    main()
