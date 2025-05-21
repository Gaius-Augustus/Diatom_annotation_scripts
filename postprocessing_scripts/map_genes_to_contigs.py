#!/usr/bin/env python3
import os
import argparse
import re

def parse_suppl_table(filename):
    """
    Parses the input tsv file with assembly contamination according to EnTAP results (for example, suppl_table6.tsv)
    and returns a list of records of the form:
      { "species": species_name, "contigs": [contig1, contig2, ...] }
    The file is expected to have a header with columns 'Species' and 'Contaminated sequences'.
    """
    species_rows = []
    with open(filename, 'r', encoding='utf-8') as f:
        header_line = f.readline().strip()
        headers = header_line.split('\t')
        try:
            species_index = headers.index("Species")
            contigs_index = headers.index("Contaminated sequences")
        except ValueError:
            raise ValueError("In the header of this file should be 'Species' and 'Contaminated sequences'.")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) <= contigs_index:
                continue
            species_name = parts[species_index].strip()
            contigs_str = parts[contigs_index].strip()
            contigs = [c.strip() for c in contigs_str.split(',') if c.strip()]
            species_rows.append({
                "species": species_name,
                "contigs": contigs
            })
    return species_rows

def parse_gtf_file(gtf_file):
    """
    Parses a GTF file and returns a dictionary mapping: contig -> set of transcripts.
    Only strings with type ‘transcript’ are considered. If there is  transcript_id key in attributes,
    it is used, otherwise the first token is used.
    """
    mapping = {}
    transcript_re = re.compile(r'transcript_id "([^"]+)"')
    with open(gtf_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            contig = parts[0]
            attributes = parts[8].strip()
            if feature == "transcript":
                transcript = None
                if "transcript_id" in attributes:
                    match = transcript_re.search(attributes)
                    if match:
                        transcript = match.group(1)
                else:
                    transcript = attributes.split()[0]
                if transcript:
                    mapping.setdefault(contig, set()).add(transcript)
    return mapping

def main():
    parser = argparse.ArgumentParser(
        description="A script to match contigs from the TSV file with transcripts from GTF files for each species."
    )
    parser.add_argument(
        "-tsv", "--tsv_table", required=True,
        help="Path to the TSV file."
    )
    parser.add_argument(
        "-r", "--root_dir", default=".",
        help="Root directory with the folder with species (it is the current working directory by default)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path to the output TSV-file."
    )
    args = parser.parse_args()

    species_rows = parse_suppl_table(args.suppl_table)
    output_lines = []
    # Header for the output file
    output_lines.append("species\tcontig\tgene")
    
    for row in species_rows:
        species_name = row["species"]
        contigs = row["contigs"]
        # Build a path to the output directory
        species_folder = os.path.join(args.root_dir, species_name)
        if not os.path.isdir(species_folder):
            alt_species_folder = os.path.join(args.root_dir, species_name.replace(" ", "_"))
            if os.path.isdir(alt_species_folder):
                species_folder = alt_species_folder
            else:
                print(f"Warning: Folder for view ‘{species_name}’ not found. Skip.")
                continue
        # Build a path to the GTF insight braker #TODO: reaplace braker folder for the future runs
        gtf_path = os.path.join(species_folder, "braker", "filtered_braker.gtf")
        if not os.path.isfile(gtf_path):
            print(f"Warning: GTF file for type ‘{species_name}’ was not found on path {gtf_path}. Skip.")
            continue
        # Parse GTF for the current species
        mapping = parse_gtf_file(gtf_path)
        # For each contig from the TSV look for transcripts
        for contig in contigs:
            transcripts = mapping.get(contig, set())
            transcripts_str = ",".join(sorted(transcripts)) if transcripts else ""
            output_lines.append(f"{species_name}\t{contig}\t{transcripts_str}")

    with open(args.output, 'w', encoding='utf-8') as out:
        out.write("\n".join(output_lines))
    print(f"Reuslt is written into the {args.output}")

if __name__ == "__main__":
    main()

