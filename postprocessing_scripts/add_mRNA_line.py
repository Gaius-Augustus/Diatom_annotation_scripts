#!/usr/bin/env python3

import os
import argparse

# Initialize a global counter for introns
intron_counter = 1

def add_mrna_lines(input_gff, output_gff):
    global intron_counter  # Ensure we can access the global counter
    with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
        gene_info = {}
        transcript_info = {}
        transcript_features = []
        exon_features = {}
        cds_features = {}

        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue  # Skip malformed lines

            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            start, end = int(start), int(end)
            attributes_dict = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes.split(';') if '=' in attr}

            if feature_type == "gene":
                # Before writing a new gene, flush any collected transcript features
                if transcript_features:
                    flush_transcript_features(outfile, transcript_info, transcript_features, exon_features, cds_features)
                    transcript_features = []
                    transcript_info = {}
                    exon_features = {}
                    cds_features = {}

                # Write the gene line immediately
                outfile.write(line)

            elif feature_type == "exon":
                parent_transcript = attributes_dict.get('Parent')
                if parent_transcript:
                    if parent_transcript not in exon_features:
                        exon_features[parent_transcript] = []
                    exon_features[parent_transcript].append((start, end, line))

            elif feature_type == "CDS":
                parent_transcript = attributes_dict.get('Parent')
                if parent_transcript:
                    if parent_transcript not in cds_features:
                        cds_features[parent_transcript] = []
                    cds_features[parent_transcript].append((start, end, line))  # Store the CDS feature line

            elif feature_type in ["start_codon", "stop_codon", "intron"]:
                transcript_features.append(line)

        # Flush the remaining transcript features for the last gene
        if transcript_features:
            flush_transcript_features(outfile, transcript_info, transcript_features, exon_features, cds_features)

def flush_transcript_features(outfile, transcript_info, transcript_features, exon_features, cds_features):
    # Write the mRNA feature for each transcript and adjust exon, CDS, and intron features
    for transcript_id, exons in exon_features.items():
        # Calculate mRNA boundaries from exons or CDS features
        mrna_start = min(exon[0] for exon in exons)
        mrna_end = max(exon[1] for exon in exons)
        
        if transcript_id in cds_features:
            mrna_start = min(mrna_start, min(cds[0] for cds in cds_features[transcript_id]))
            mrna_end = max(mrna_end, max(cds[1] for cds in cds_features[transcript_id]))

        # Write the mRNA feature before writing exons/CDS/introns
        write_mrna(outfile, transcript_id, exons[0][2], mrna_start, mrna_end)

        # Now adjust exons based on CDS and write exons/introns
        for exon_start, exon_end, exon_line in exons:
            cds_regions = [(cds_start, cds_end) for cds_start, cds_end, _ in cds_features.get(transcript_id, [])]
            adjusted_exon = adjust_exon_to_cds(exon_start, exon_end, cds_regions)

            if adjusted_exon:
                # Write the corrected exon
                write_feature_with_new_coords(outfile, exon_line, adjusted_exon[0], adjusted_exon[1])

                # If there's an intron part, write it
                if adjusted_exon[2]:  # intron_start and intron_end exist
                    create_and_write_intron(outfile, exon_line, adjusted_exon[2], adjusted_exon[3])

        # Write the CDS features for the transcript
        if transcript_id in cds_features:
            for cds_start, cds_end, cds_line in cds_features[transcript_id]:
                outfile.write(cds_line)

    # Write all other features (start_codon, stop_codon, etc.)
    for feature in transcript_features:
        outfile.write(feature)

def write_mrna(outfile, transcript_id, example_line, start, end):
    """
    Write the mRNA feature to the output file, using the example line for the structure.
    """
    fields = example_line.strip().split('\t')
    fields[2] = 'mRNA'  # Set the feature type to mRNA
    fields[3] = str(start)
    fields[4] = str(end)
    attributes = f"ID={transcript_id};Parent={transcript_id.split('.')[0]}"  # Adjust for correct ID and Parent
    fields[8] = attributes
    outfile.write("\t".join(fields) + "\n")

def adjust_exon_to_cds(exon_start, exon_end, cds_regions):
    """
    Adjust exon coordinates to match the CDS coordinates. 
    Return the adjusted exon coordinates and any leftover regions as introns.
    """
    for cds_start, cds_end in cds_regions:
        # Exon completely overlaps with CDS
        if exon_start >= cds_start and exon_end <= cds_end:
            return (exon_start, exon_end, None, None)  # No intron part

        # Exon partially overlaps with CDS on the left side
        elif exon_start < cds_start and exon_end <= cds_end and exon_end > cds_start:
            return (cds_start, exon_end, exon_start, cds_start - 1)  # Intron on the left

        # Exon partially overlaps with CDS on the right side
        elif exon_start >= cds_start and exon_start < cds_end and exon_end > cds_end:
            return (exon_start, cds_end, cds_end + 1, exon_end)  # Intron on the right

        # Exon surrounds the CDS completely (creates two introns)
        elif exon_start < cds_start and exon_end > cds_end:
            return (cds_start, cds_end, exon_start, cds_start - 1)  # Intron left, and exon_end -> intron

    # If no match with CDS, return exon as intron
    return (None, None, exon_start, exon_end)

def write_feature_with_new_coords(outfile, original_line, new_start, new_end):
    """
    Write the feature with adjusted coordinates.
    """
    fields = original_line.strip().split('\t')
    fields[3] = str(new_start)
    fields[4] = str(new_end)
    outfile.write("\t".join(fields) + "\n")

def create_and_write_intron(outfile, original_line, intron_start, intron_end):
    """
    Create a new intron feature with a unique ID and write it to the file.
    """
    global intron_counter
    fields = original_line.strip().split('\t')
    
    # Modify the ID attribute to replace agat-exon with fixed-intron and append a counter
    attributes = fields[8].replace("agat-exon", "fixed-intron")
    new_id = f"{attributes};ID=fixed-intron-{intron_counter}"
    fields[8] = new_id
    fields[2] = "intron"
    fields[3] = str(intron_start)
    fields[4] = str(intron_end)
    
    # outfile.write("\t".join(fields) + "\n") # realized the intron was in the file, already
    
    # Increment the intron counter for uniqueness
    intron_counter += 1

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Add missing mRNA lines and fix exon boundaries in a GFF3 file, this script was written to fix an AGAT-created gff3 file; AGAT messed up partial genes.')
    parser.add_argument('-i', '--input', required=True, help='Input GFF3 file.')
    parser.add_argument('-o', '--output', required=True, help='Output GFF3 file.')

    # Parse arguments
    args = parser.parse_args()

    # Check if input is a file
    if os.path.isfile(args.input):
        # Process the specified file
        add_mrna_lines(args.input, args.output)
        print(f"Processed: {args.input} -> {args.output}")
    else:
        print("Invalid input path. Please provide a valid GFF3 file.")

if __name__ == "__main__":
    main()
