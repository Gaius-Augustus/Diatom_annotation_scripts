#!/usr/bin/env python3

import re
import argparse

def fix_product_name(product_name):
    # Fix uncharacterized protein
    if product_name.startswith("uncharacterized protein"):
        product_name = "putative protein"
    # if a product name ends in putative, remove that
    product_name = re.sub(r'putative$', '', product_name, flags=re.IGNORECASE)
    
    # if a product name ends in uncharacterised, remove that
    product_name = re.sub(r'uncharacterised$', '', product_name, flags=re.IGNORECASE)
    # remove trailing whitespace
    product_name = product_name.strip()

    # Fix hypothetical protein only if followed by additional text, keeping the rest of the text
    if product_name.startswith("hypothetical protein "):
        product_name = product_name.replace("hypothetical protein", "putative protein", 1)
    
    # Fix proteins protein
    product_name = re.sub(r'proteins protein', 'protein', product_name)

    # Fix plural forms (simple heuristic)
    product_name = re.sub(r'proteins', 'protein', product_name)

    # Remove suspicious phrases
    product_name = re.sub(r'truncat(ed|ion)?', '', product_name)
    product_name = re.sub(r'Truncated ', r'Nonfunctional ', product_name)
    product_name = re.sub(r'Fragment', r'', product_name, flags=re.IGNORECASE)
    product_name = re.sub(r'partial', r'', product_name, flags=re.IGNORECASE)

    # Remove strings containing three or more consecutive digits (may be part of a larger string)
    product_name = re.sub(r'\S*\d{3,}\S*', r'', product_name)

    product_name = product_name.replace('_', ' ')  # Replace underscores with spaces

    # Change homologs/paralogs to '-like protein'
    if not(re.search(r'-like', product_name)):
        product_name = re.sub(r'\b(Homolog|paralog)\b', r'-like protein', product_name)
        product_name = re.sub(r'\shomolog', r'-like protein', product_name)
    else:
        product_name = re.sub(r'\b(Homolog|paralog)\b', r'', product_name)
        product_name = re.sub(r'\shomolog', r'', product_name)
    # also remove like from chloroplastic-like and mitochondrial-like
    product_name = re.sub(r'chloroplastic-like', r'chloroplastic', product_name)
    product_name = re.sub(r'mitochondrial-like', r'mitochondrial', product_name)
    # remove -like proteins from protein -like protein
    product_name = re.sub(r'protein -like protein', r'protein', product_name)

    # if a product name contains --, make it one -
    product_name = re.sub(r'--', '-', product_name)

    # someone seriously submitted a homeolog, replace by nothing
    product_name = re.sub(r'homeolog', '', product_name)

    # make forms singular
    product_name = re.sub(r'forms', 'form', product_name)

    # if a product name contains two times -like, remove the last -like
    product_name = re.sub(r'-like(.*)-like', r'\1-like', product_name)

    # if a product name contains -like putative, remove putative
    product_name = re.sub(r'-like putative', r'-like', product_name)

    # make yippee-like uniformal, if there is no dash, introduce it
    product_name = re.sub(r'yippee like', r'yippee-like', product_name)

    # if a product name ends in like \d+, remove the digits
    product_name = re.sub(r'like \d+$', r'like', product_name)
    # if a product name ends in containing \d+ remove digits
    product_name = re.sub(r'containing \d+$', r'containing', product_name)

    # Shorten long product names to 100 characters
    if len(product_name) > 100:
        # Find the last space before the 100th character
        last_space = product_name[:100].rfind(' ')
        # replace product_name with the first 100 characters up to the last space
        product_name = product_name[:last_space].rstrip()

    # Replace 'gene' with 'protein'
    product_name = re.sub(r'\bgene\b', r'protein', product_name)

    # Fix product names starting with "and"
    if re.match(r'\s?and\s', product_name):
        product_name = product_name[4:].strip()

    # Fix product names starting with "a"
    if re.match(r'\s?a\s', product_name):
        product_name = product_name[2:].strip()

    # Remove "FOG" from product names
    product_name = product_name.replace("FOG", "")

    # If a product name contains gp, replace the entire product name with "Hypothetical protein"
    if re.search(r'gp', product_name):
        product_name = "hypothetical protein"

    # Replace "possibly", "residue", and "unnamed" with empty strings or replacements
    product_name = re.sub(r'\bpossibly\b', r'', product_name, flags=re.IGNORECASE)
    product_name = re.sub(r'\bresidue\b', r'', product_name, flags=re.IGNORECASE)
    product_name = re.sub(r'\bunnamed\b', r'putative protein', product_name, flags=re.IGNORECASE)

    # if product ends in family, replace with family member
    product_name = re.sub(r'family$', r'family member', product_name, flags=re.IGNORECASE)

    # fix plurals
    product_name = re.sub(r'ases', r'ase', product_name)
    product_name = re.sub(r'proteins', r'protein', product_name)
    product_name = re.sub(r'units', r'unit', product_name)
    product_name = re.sub(r'ac-diamides', r'ac-diamide', product_name)
    product_name = re.sub(r'chaperonins', r'chaperonin', product_name)
    product_name = re.sub(r'complexes', r'complex', product_name)
    product_name = re.sub(r'condensins', r'condensin', product_name)
    product_name = re.sub(r'copines', r'copine', product_name)
    product_name = re.sub(r'biotransformers', r'biotransformer', product_name)
    product_name = re.sub(r'enzymes', r'enzyme', product_name)
    product_name = re.sub(r'cyclins', r'cyclin', product_name)
    product_name = re.sub(r'cylophilins', r'cylophilin', product_name)
    product_name = re.sub(r'channels', r'channel', product_name)
    product_name = re.sub(r'mutants', r'mutant', product_name)
    product_name = re.sub(r'receptors', r'receptor', product_name)
    product_name = re.sub(r'retroposons', r'retroposon', product_name)
    product_name = re.sub(r'factors', r'factor', product_name)
    product_name = re.sub(r'members', r'member', product_name)
    product_name = re.sub(r'antigens', r'antigen', product_name)

    # if something comes afte -like, move the -like to the end of the string
    product_name = re.sub(r'(.*)-like(.*)', r'\1\2-like', product_name)
    
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
