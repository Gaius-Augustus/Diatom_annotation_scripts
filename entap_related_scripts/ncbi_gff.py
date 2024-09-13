import subprocess
import csv
import re
import json
import sys

# Define inputs
annotation_file = sys.argv[1]
gff_file = sys.argv[2]
pfam_file = sys.argv[3]

# Define outputs
output_file = sys.argv[4]
accessions_file = 'accessions.txt'
refseq_json = 'refseq.json'
refseq_file = 'refseq.tsv'

# Define cleaning function for product descriptions
def clean_product_description(description):
    # Convert specific descriptions to 'hypothetical protein'
    if any(term in description.lower() for term in ["uncharacterized", "unknown", "low quality protein","predicted protein", "pseudo", "clone"]):
        return "hypothetical protein"
    
    # Remove trailing dashes
    description = re.sub(r'-$', '', description)
    
    # Remove occurrences of '. '
    description = re.sub(r'\. ', ' ', description)
    
    # Remove unbalanced brackets  
    description = re.sub(r'[^\w\s-]', '', description)
    description = re.sub(r'^-+', '', description).strip()
    
    # Convert descriptions that are all numbers to 'hypothetical protein'
    if description.isdigit():
        return "hypothetical protein"
    
    return description

# Step 1: Extract unique accession IDs and save to accessions.txt
accessions = set()
try:
    with open(annotation_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            accession = row.get('Subject Sequence', '').strip()
            if accession and accession != 'Subject Sequence':
                accessions.add(accession)
except FileNotFoundError:
    print(f"Error: The file {annotation_file} was not found.")
    exit(1)

with open(accessions_file, 'w') as file:
    for accession in sorted(accessions):
        file.write(f"{accession}\n")

# Step 2: Process accessions in chunks
command = f"datasets summary gene accession --report gene --inputfile {accessions_file} > {refseq_json}"
subprocess.run(command, shell=True)

with open(refseq_json, 'r') as file:
    data = json.load(file)

# Step 3: Write gene_id, query, and description to refseq.tsv
with open(refseq_file, 'w') as outfile:
    gene_id_to_info = {}
    for report in data.get("reports", []):
        gene_id = report["gene"]["gene_id"]
        description = report["gene"]["description"]
        queries = report.get("query", [])
        if gene_id not in gene_id_to_info:
            gene_id_to_info[gene_id] = {'description': description, 'queries': set()}
        for query in queries:
            gene_id_to_info[gene_id]['queries'].add(query)
    
    for gene_id, info in gene_id_to_info.items():
        description = info['description']
        for query in info['queries']:
            outfile.write(f"{query}\t{gene_id}\t{description}\n")

# Step 4: Load refseq.tsv into a dictionary for lookup
refseq_dict = {}
try:
    with open(refseq_file, 'r') as refseq:
        for line in refseq:
            cols = line.strip().split('\t')
            if len(cols) == 3:
                refseq_dict[cols[0]] = {'gene_id': cols[1], 'description': cols[2]}
except FileNotFoundError:
    print(f"Error: The file {refseq_file} was not found.")
    exit(1)

# Step 5: Load PFAM mappings into a dictionary
pfam_dict = {}
try:
    with open(pfam_file, 'r') as pfam:
        for line in pfam:
            cols = line.strip().split('\t')
            if len(cols) >= 6:
                pfam_dict[cols[5]] = cols[4]
except FileNotFoundError:
    print(f"Error: The file {pfam_file} was not found.")
    exit(1)

# Step 6: Process the GFF file
def get_dbxref_prefix(subject_sequence):
    if re.match(r'^(AC_|NC_|NG_|NT_|NW_|NZ_|NM_|NR_|XM_|XR_|AP_|NP_|YP_|XP_|WP_)', subject_sequence):
        return 'GeneID'
    elif subject_sequence.startswith('sp|') and '|' in subject_sequence:
        return 'UniProtKB/Swiss-Prot'
    elif subject_sequence.startswith('tr|') and '|' in subject_sequence:
        return 'UniProtKB/TrEMBL'
    return None

def extract_uniprot_id(subject_sequence):
    if '|' in subject_sequence:
        parts = subject_sequence.split('|')
        if len(parts) > 1:
            return parts[1]
    return subject_sequence

annotations, subject_sequences, go_terms, contaminants, eggnog, pfam = {}, {}, {}, {}, {}, {}

try:
    with open(annotation_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        entry_count = 0
        for row in reader:
            gene_id = row['Query Sequence'].strip()
            subject_sequence = row['Subject Sequence'].strip()
            eggnog_domain = row['EggNOG Protein Domains'].strip()

            if eggnog_domain.startswith('PFAM'):
                pfam_ids = eggnog_domain.split('(')[1].split(')')[0].split(', ')
                pfam[gene_id] = pfam_ids
            else:
                pfam_ids = eggnog_domain.split(', ')
                pfam[gene_id] = pfam_ids 
            
            contam = row.get('Contaminant', '').strip()
            eggnog_description = row.get('EggNOG Description', '').strip()
            eggnog[gene_id] = eggnog_description
            annotation = ' '.join(row['Description'].split()[1:]).strip()
            # Find the position of 'OS=' in the string (UniProt)
            os_index = annotation.find('OS=')
            # If 'OS=' is found, truncate the string to that point
            if os_index != -1:
                annotation = annotation[:os_index].strip()

            annotations[gene_id] = annotation
            subject_sequences[gene_id] = subject_sequence
            contaminants[gene_id] = contam

            go_term_values = [row.get('UniProt GO Biological', '').strip(),
                              row.get('UniProt GO Cellular', '').strip(),
                              row.get('UniProt GO Molecular', '').strip()]
            go_term_values = [value for value in go_term_values if value != 'NaN']
            go_terms[gene_id] = ','.join(go_term_values).rstrip(',')

            entry_count += 1

except FileNotFoundError:
    print(f"Error: The file {annotation_file} was not found.")
    exit(1)
except KeyError as e:
    print(f"Error: Missing expected header in {annotation_file}: {e}")
    exit(1)
except Exception as e:
    print(f"An error occurred while reading {annotation_file}: {e}")
    exit(1)

# Step 7: Count the number of transcripts and store start and end positions
transcript_count = 0
transcript_to_gene_map = {}
gene_to_transcript_map = {}
transcript_positions = {}  # Dictionary to store transcript positions

try:
    with open(gff_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) > 2 and fields[2] in ('transcript', 'mRNA'):
                transcript_count += 1
                feature_type = fields[2]
                attributes = fields[8]

                id_match = re.search(r'ID=([^;]+)', attributes)
                parent_match = re.search(r'Parent=([^;]+)', attributes)

                if id_match and parent_match:
                    transcript_id = id_match.group(1)
                    gene_id = parent_match.group(1)
                    transcript_to_gene_map[transcript_id] = gene_id
                    if gene_id not in gene_to_transcript_map:
                        gene_to_transcript_map[gene_id] = []
                    gene_to_transcript_map[gene_id].append(transcript_id)

                    # Store the start and end positions of the transcript
                    start_pos = int(fields[3])
                    end_pos = int(fields[4])
                    transcript_positions[transcript_id] = (start_pos, end_pos)
                
except FileNotFoundError:
    print(f"Error: The file {gff_file} was not found.")
    exit(1)
except Exception as e:
    print(f"An error occurred while reading {gff_file}: {e}")
    exit(1)

# Step 8: Process the GFF file and add new attributes
try:
    with open(gff_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            feature_type = fields[2]
            attributes = fields[8]

            # Extract ID and Parent from attributes
            id_match = re.search(r'ID=([^;]+)', attributes)
            parent_match = re.search(r'Parent=([^;]+)', attributes)

            if feature_type in ['mRNA', 'transcript']:
                if id_match:
                    gene_id = id_match.group(1)
                    fields[2] = 'mRNA' # Change transcript to mRNA because NCBI seems to like that better
                    new_attributes = re.sub(r'ID=[^;]+', f'ID={gene_id}', attributes)
                    fields[8] = new_attributes
                else:
                    outfile.write(line)
                    continue

            elif feature_type in ['exon', 'CDS']:
                if parent_match:
                    parent_id = parent_match.group(1)
                    new_attributes = re.sub(r'Parent=[^;]+', f'Parent={parent_id}', attributes)
                    fields[8] = new_attributes
                else:
                    outfile.write(line)
                    continue

            elif feature_type == 'gene':
                if id_match:
                    gene_id = id_match.group(1)
                    
                    # Retrieve associated transcript_ids
                    transcripts = gene_to_transcript_map.get(gene_id, [])
                    
                    # Find the longest transcript
                    if transcripts:
                        longest_transcript = max(transcripts, key=lambda x: transcript_positions[x][1] - transcript_positions[x][0])
                        gene_id = longest_transcript
                else:
                    outfile.write(line)
                    continue

            dbxrefs = []
            if gene_id in subject_sequences:
                subject_sequence = subject_sequences[gene_id]
                dbxref_prefix = get_dbxref_prefix(subject_sequence)

                if dbxref_prefix == 'GeneID':
                    gene_id_value = refseq_dict.get(subject_sequence, {}).get('gene_id')
                    if gene_id_value:
                        dbxrefs.append(f"{dbxref_prefix}:{gene_id_value}")
                elif dbxref_prefix in ['UniProtKB/Swiss-Prot', 'UniProtKB/TrEMBL']:
                    uniprot_id = extract_uniprot_id(subject_sequence)
                    if uniprot_id:
                        dbxrefs.append(f"{dbxref_prefix}:{uniprot_id}")
            
            if fields[2] != 'gene':
                if gene_id in pfam:
                    pfam_ids = pfam.get(gene_id, [])
                    pfam_entries = [f"PFAM:{pfam_dict.get(pfam_id, '')}" for pfam_id in pfam_ids if pfam_dict.get(pfam_id)]
                    if pfam_entries:
                        dbxrefs.append(','.join(pfam_entries))

                if gene_id in go_terms:
                    go_value = go_terms[gene_id]
                    if go_value:
                        go_value = re.sub(r',,+', ',', go_value)
                        dbxrefs.append(go_value)

            if dbxrefs:
                attributes += f";Dbxref={','.join(dbxrefs)}"

            # Determine product description based on dbxref_prefix
            if feature_type in ['mRNA', 'transcript','CDS']:
                annotation_description = 'hypothetical protein'
                if gene_id in subject_sequences:
                    subject_sequence = subject_sequences[gene_id]
                    dbxref_prefix = get_dbxref_prefix(subject_sequence)

                    if dbxref_prefix == 'GeneID':
                        annotation_description = refseq_dict.get(subject_sequence, {}).get('description', 'hypothetical protein')
                    else:
                        annotation_description = annotations.get(gene_id, 'hypothetical protein')

                description = clean_product_description(annotation_description)
                if description:
                    attributes += f";product={description}"
                else:
                    attributes += f";product=hypothetical protein"

                if gene_id in eggnog:
                    eggnog_description = eggnog[gene_id]
                    if eggnog_description != 'NaN':
                        attributes += f";note=EggNOG:{eggnog_description}"

            if feature_type == 'gene':
                attributes += ";gene_biotype=protein_coding"

            # if gene_id in contaminants and contaminants[gene_id] == 'Yes':
            #     attributes += ";note=Contaminant"

            fields[8] = attributes
            outfile.write('\t'.join(fields) + '\n')

except FileNotFoundError:
    print(f"Error: The file {gff_file} or {output_file} was not found.")
    exit(1)
except Exception as e:
    print(f"An error occurred while processing files: {e}")
    exit(1)

# Step 9: Output counts and check for consistency
print(f"Number of entries in {annotation_file}: {entry_count}")
print(f"Number of transcripts and mRNA in {gff_file}: {transcript_count}")

if entry_count == transcript_count:
    print("The number of entries in the annotation file matches the number of transcripts in the GFF file.")
else:
    print("Warning: The number of entries in the annotation file does not match the number of transcripts in the GFF file.")
print(f"Updated GFF file written to {output_file}")

