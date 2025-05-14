# Processing of functionally decorated gff3 file

Contact: katharina.hoff@uni-greifswald.de

## Filtering single-exon genes

On the basis of functional annotation with EnTAP, results of an initial OrthoFinder run (with only the newly annotated protein sets), and a DIAMOND search against NCBI NR, the braker.gtf file was filtered as follows:

We extract the assigned homologs for each species from the `OrthoFinder/Orthogroups/Orthogroups.tsv` file:

```
./write_list_files.py
```

This generates a subfolder `output_files` with one file per species.

```
# Make intersection
overlapStat.pl single_exon_genes_without_hit.txt orthofinder_assigned_homolog.lst --outfiles
# the file orthofinder_assigned_homolog.lst is located in the Orthogroups folder of OrthoFinder, one file per species
    
# Rename the output file and clean up
mv combset.10.lst to_be_removed_before_entap.lst
rm combset*
    
# Get genes not annotated by entap
tail -n +2 entap_outfiles/final_results/unannotated.tsv | cut -f1 > no_entap.lst
    
# Make intersection
overlapStat.pl to_be_removed_before_entap.lst no_entap.lst --outfiles
    
# Rename the output file and clean up
mv combset.11.lst to_be_removed.lst
rm combset*
    
# Remove the bad IDs from the actual braker.gtf file
grep -v -f to_be_removed.lst ../braker/braker.gtf | gffread -T -o braker_filtered.gtf - &> /dev/null
    
# Generate annotation fasta from the filtered GTF file
getAnnoFastaFromJoingenes.py -g ../genome/genome.fa -f braker_filtered.gtf -o braker_filtered_raw -s True # produces stops.lst

# remove genes that have in-frame stop codons
grep -v -f stops.lst braker_filtered.gtf | gffread -T -o braker_filtered2.gtf - &> /dev/null
mv braker_filtered2.gtf braker_filtered.gtf
```

## Filtering gff3 file with functional decorations from EnTAP (generated from raw braker.gtf output)

```
grep ">" braker_filtered.aa" | perl -pe 's/>//' > good_tx.lst
filter_genes_from_uconn_gff3.py -g decorated.gff3 -l good_tx.lst -o braker.gff3
grep -v -P "\tAGAT\t" braker.gff3 > braker_noagat.gff3
add_mRNA_line.py -i decorated.gff3 -o mRNA.gff3
fix_product_names_ncbi.py -i mRNA.gff3 -o fixed_names.gff3
fix_Dbxref_attributes_in_genes.py -i fixed_names.gff3 -o fixed_dbxref.gff3
```
## Tagging contaminations in the GFF3 file

Contamination information was extracted from the EnTAP output as follows:

```
# List of species
species_list="Asterionella_formosa Craspedostauros_australis Cyclotella_cryptica Epithemia_pelagica Mediolabrus_comicus Pseudo-nitzschia_multiseries Skeletonema_tropicum Thalassiosira_exigua Thalassiosira_pacifica Asterionellopsis_glacialis Cyclostephanos_invisitatus Cylindrotheca_fusiformis Fistulifera_pelliculosa Nitzschia_palea Pseudo-nitzschia_pungens Stephanocyclus_meneghinianus Thalassiosira_gravida Thalassiosira_profunda Bacterosira_constricta Cyclostephanos_tholiformis Detonula_confervacea Fistulifera_solaris Nitzschia_putrida Skeletonema_costatum Stephanodiscus_minutulus Thalassiosira_livingstoniorum Thalassiosira_sundarbana Chaetoceros_muellerii Cyclotella_atomus Discostella_pseudostelligera Fragilaria_radians Porosira_glacialis Skeletonema_marinoi Stephanodiscus_triporus Thalassiosira_mediterranea Conticribra_guillardii Cyclotella_baltica Discostella_stelligera Fragilariopsis_cylindrus Psammoneis_japonica Skeletonema_menzelii Thalassiosira_allenii Thalassiosira_oceanica Conticribra_weissflogii Cyclotella_choctawhatcheeana Discostella_stelligeroides Licmophora_abbreviata Pseudo-nitzschia_delicatissima Skeletonema_potamos Thalassiosira_delicatula Thalassiosira_ordinaria"

# Loop through each species
for S in $species_list; do
    echo "$S"
    cd /nas-hs/projs/diatom-dl/braker-snake/data/entap/$S/braker/entap_outfiles/final_results
O=contam_contigs
mkdir ${O}
tail -n +2 annotated_contam.tsv | cut -f1 | sort -u > ${O}/contam_tx.lst
tail -n +2 annotated_without_contam.tsv | cut -f1 | sort -u > ${O}/no_contam_tx.lst
# retrieve the contigs of the transcripts including a count, i.e. how often have we seen this contig in what category?

# Extract contigs for contaminated transcripts
grep -F -f ${O}/contam_tx.lst /nas-hs/projs/diatom-dl/braker-snake/data/species/${S}/braker/braker.gtf | \
awk '{print $1"\t"$12}' | \
sed 's/"//g; s/;//g' | sort -u > ${O}/contam_contigs_full.tsv

# Extract contigs for non-contaminated transcripts
grep -F -f ${O}/no_contam_tx.lst /nas-hs/projs/diatom-dl/braker-snake/data/species/${S}/braker/braker.gtf | \
awk '{print $1"\t"$12}' | \
sed 's/"//g; s/;//g' | sort -u > ${O}/no_contam_contigs_full.tsv

# Count occurrences of contaminated contigs
cut -f1 ${O}/contam_contigs_full.tsv | sort | uniq -c | awk '{print $2"\t"$1}' > ${O}/contam_contigs_counts.txt

# Count occurrences of non-contaminated contigs
cut -f1 ${O}/no_contam_contigs_full.tsv | sort | uniq -c | awk '{print $2"\t"$1}' > ${O}/no_contam_contigs_counts.txt

# Define input files
contam_file="${O}/contam_contigs_counts.txt"
no_contam_file="${O}/no_contam_contigs_counts.txt"
output_file="${O}/contig_contamination_percentage.txt"

# Prepare a temporary file with all contigs and their counts
join -a1 -a2 -e 0 -o 0,1.2,2.2 -t $'\t' \
    <(sort "${contam_file}") \
    <(sort "${no_contam_file}") > "${O}/joined_contig_counts.txt"

# Process the joined file and compute contamination percentages
awk -F '\t' '{
    contam = $2;      # Contaminated count
    no_contam = $3;   # Non-contaminated count
    contig = $1;      # Contig name

    if (no_contam == 0) {
        # Contig exists only in contaminated
        print contig "\t100%" > "'${output_file}'"
    } else {
        # Calculate contamination percentage
        perc = (contam / (contam + no_contam)) * 100;
        print contig "\t" perc "%" > "'${output_file}'"
    }
}' "${O}/joined_contig_counts.txt"

# Clean up temporary files (optional)
rm "${O}/joined_contig_counts.txt"

echo "Results written to ${output_file}"

# Define input files
contam_file="${O}/contam_contigs_counts.txt"
no_contam_file="${O}/no_contam_contigs_counts.txt"
output_file="${O}/contig_contamination_percentage.txt"
filtered_output_file="${O}/high_contamination_contigs.txt"

# Prepare a temporary file with all contigs and their counts
join -a1 -a2 -e 0 -o 0,1.2,2.2 -t $'\t' \
    <(sort "${contam_file}") \
    <(sort "${no_contam_file}") > "${O}/joined_contig_counts.txt"

# Process the joined file and compute contamination percentages
awk -F '\t' '{
    contam = $2;      # Contaminated count
    no_contam = $3;   # Non-contaminated count
    contig = $1;      # Contig name

    if (no_contam == 0) {
        # Contig exists only in contaminated
        perc = 100;
    } else {
        # Calculate contamination percentage
        perc = (contam / (contam + no_contam)) * 100;
    }

    # Output all results with percentages
    print contig "\t" perc "%" > "'${output_file}'"

    # Filter for high contamination (greater than 75%)
    if (perc > 75) {
        print contig > "'${filtered_output_file}'"
    }
}' "${O}/joined_contig_counts.txt"

# Clean up temporary files (optional)
rm "${O}/joined_contig_counts.txt"

echo "Filtered results (contigs with >75% contamination) written to ${filtered_output_file}"
done
```

The resulting files were used to decorate the GFF3 files as follows:

```
for file in NCBI_gffs_filtered/*_ncbi.no_agat.mRNA.fixednames.dbxref.gff
do
    # Strip the suffix to get the base name (e.g., "Asterionella_formosa")
    base_name="${file%_ncbi.no_agat.mRNA.fixednames.dbxref.gff}"
    base_name="$(basename "$base_name")"

    # Path to the high-contamination contigs list
    contam_file="entap/${base_name}/braker/entap_outfiles/final_results/contam_contigs/high_contamination_contigs.txt"

    # Construct an output file name (you can change to suit your needs)
    output_gff="NCBI_gffs_filtered/${base_name}_ncbi.no_agat.mRNA.fixednames.dbxref.withContaminationTag.gff"

    # Check if the contamination file exists
    if [[ -f "$contam_file" ]]; then
        echo "Processing $file with contamination list $contam_file ..."
        awk -F'\t' -v OFS='\t' -v cf="$contam_file" '
            BEGIN {
                # If the file exists, load high-contamination contig names into an array
                while ((getline line < cf) > 0) {
                    contam[line] = 1
                }
                close(cf)
            }
            {
                # If line doesn't look like GFF feature lines (e.g. comment or FASTA), just print
                if ($1 ~ /^#/ || NF < 9) {
                    print; next
                }

                # 1) Capitalize any existing note= to Note=
                gsub(/note=/, "Note=", $9);

                # 2) If this line is on a contaminated contig and is a CDS feature,
                #    append or create the contamination note in the attributes.
                if (($1 in contam) && ($3 == "CDS")) {
                    if ($9 ~ /Note=/) {
                        # Append to existing Note
                        sub(/(Note=[^;"]*)/, "&, The genomic contig of this gene structure is likely a contamination", $9)
                    } else {
                        # Create a new Note attribute
                        $9 = $9 ";Note=The genomic contig of this gene structure is likely a contamination"
                    }
                }

                print
            }
        ' "$file" > "$output_gff"
    else
        # If no contamination file is found, we still need to capitalize 'note=' -> 'Note='
        echo "No contamination file found for $base_name; still capitalizing note fields ..."
        # Run awk to only capitalize note= to Note=, no contamination logic
        awk -F'\t' -v OFS='\t' '
            {
                if ($1 ~ /^#/ || NF < 9) {
                    print; next
                }
                # Capitalize any existing note= to Note=
                gsub(/note=/, "Note=", $9)
                print
            }
        ' "$file" > "$output_gff"
    fi
done
```
## Horizontal Gene Transfer Analysis

In addition to utilizing the data already gathered by EnTAP in the annotation phase, every diatom species was blasted (with DIAMOND version 2.1.8, same version as the primary annotation) against donor and recipient databases with the following general DIAMOND command:

```
diamond-2.1.8 blastp -o path/to/output/file -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle --very-sensitive -p 8 --max-target-seqs 3 --evalue 0.000010 --subject-cover 50.000000 --query-cover 50.000000 -q path/to/fasta/file -d path/to/diamond/database
```

Note the parameters of evalue (0.000010), subject coverage (50), query coverage (50), sensitivity (very-sensitive).

The donor databases:
* NCBI Refseq Bacteria Protein 228
* NCBI Refseq Plant Protein 216

The recipient database:
* NCBI Refseq Orchrophyta with the Diatoms removed
    * This was generated through NCBI query `txid2696291[organism:exp] AND refseq[filter] NOT txid2836[organism:exp]` on January 29th 2025 resulting in 34,013 proteins at the time.

The results of running DIAMOND against the donor and recipient databases were then analyzed using the following criteria to determine potential HGT candidates. The longest isoform was used to represent the primary gene for HGT analysis.

* If gene aligned against the bacteria donor database, and not plant donor
* If gene did not align against the recipient database

If all of the above were true, this provided the list of potential candidates. The list was further refined using the following criteria with locational data from the GFF. 

* Remove HGT candidate if it does not have two flanks/neighboring genes
* Remove HGT candidate if either of the flanks were flagged as a contaminant
* Contaminant defined as the primary annotation taxonomy belonging to bacteria or fungi lineage
* Remove HGT candidate if either flanks are an HGT candidate
* Remove HGT candidate if either flanks have an alignment against a donor database

Once all of the additional filtering was complete, we were left with the final HGTs.

## Retrieving longest isoform for final OrthoFinder analysis (from braker.aa)

```
get_longest_isoform_from_braker_aa.py -i braker_filtered.aa -o braker_longest.aa
```

The contents of species-specific fixed_dbxref.gff3 files are provided at Zenodo with doi 10.5281/zenodo.13745090

## Clean faa files (remove contaminated genes and HGT genes)
### Step 1: Create file with contaminated contigs

```
# make gtf files from gff3 if necessary with
# gffread input.gff3  -T  -o output.gtf #conda install -c bioconda gffread

export=~/OrthoFinder/species_braker_gtf_only # I made a separete directory for files, but it is not necessary; 
cd $WORKDIR 
# here should be a set of folders with braker_filtered.gtf files, for example: Asterionella_formosa/braker/braker_filtered.gtf

# a. map genes to contigs

# suppl_table6.tsv: https://docs.google.com/spreadsheets/d/1s2khPazgCiE_juLIOb_cU-yog9vZF4eNTCTblEtn3qg/edit?gid=0#gid=0
python3 ./map_genes_to_contigs.py -s suppl_table6.tsv -r ./ -o ./mapped_genes.tsv

# mapped_genes.tsv:
#species contig  gene
#Asterionella formosa    NKIB01000138.1  g191.t1
#Asterionella formosa    NKIB01000374.1  g549.t1
#Asterionella formosa    NKIB01000390.1  g570.t1

# b. prepare tsv file with HGT genes
# suppl_table7.tsv: https://docs.google.com/spreadsheets/d/1uHG6mxcYsl9Ar1GymwgRBxbjiODnbkWmPeyZ7sD6bFA/edit?gid=1437366730#gid=1437366730
cat suppl_table7.tsv | cut -f1,8 | head --lines=-2 > tmp.tsv
./add_HGT_column.sh tmp.tsv suppl_table7.tmp.tsv

# suppl_table7.tmp.tsv
#Species HGT     Final Transcript IDs in filtered GFF3
#Asterionella formosa    HGT     g8223.t1,g5938.t1,g1909.t1,g2541.t1,g2541.t2,g11180.t1,g8789.t1
#Asterionellopsis glacialis      HGT     g8351.t1,g4875.t1,g3876.t1,g7383.t1,g14445.t1,g5941.t1,..
#Bacterosira constricta  HGT     
#Chaetoceros muellerii   HGT     g4117.t1,g3341.t1,g2374.t1,g523.t1,g10377.t1,g10377.t2,..
#Conticribra guillardii  HGT     g4925.t1,g2580.t1

# c. Concatenate file with HGT genes and file with genes from contaminated contigs
tail -n +2 suppl_table7.tmp.tsv >> mapped_genes.tsv
# sort -o mapped_genes.tsv mapped_genes.tsv
# rm suppl_table7.tmp.tsv

# mapped_genes.tsv
# species contig  gene
# Asterionella formosa    HGT     g8223.t1,g5938.t1,g1909.t1,g2541.t1,g2541.t2,g11180.t1,g8789.t1
# Asterionella formosa    NKIB01000138.1  g191.t1
# Asterionella formosa    NKIB01000374.1  g549.t1
```

### Step 2: Remove contaminants and HTG genes
```
ANNOTATIONS_DIR="/nas-hs/projs/diatom-dl/braker-snake/data/Bacillariophyta_annotations"  	# input folder with .gff3.gz or .gff3 files
LONGEST_ISOFORMS_DIR="longest_isoforms"                  	# temporary folder for GTF files with longest isoforms

CLEANED_LONGEST_ISOFORMS_DIR="cleaned_longest_isoforms" 	# folder with GTF files without contaminated genes
FAA_OUTPUT_DIR="cleaned_faa_files"               					# output folder for faa files (will be used as input for busco)
GENOME_FASTA_DIR="/nas-hs/projs/diatom-dl/braker-snake/data/species/"   # Folder with genome FASTA files

# mapped genes file (with header)
MAPPED_GENES_FILE="/home/natalia/busco_input/mapped_genes.tsv"
# temporary file for gene IDs to exclude
EXCLUDE_LIST="/tmp/exclude_genes.lst"

# tools
LONGEST_ISOFORM_SCRIPT="/home/natalia/TSEBRA/bin/get_longest_isoform.py"
GET_ANNO_FASTA_SCRIPT="/home/natalia/Augustus/scripts/getAnnoFastaFromJoingenes.py"  # GFF-to-FAA script
PYTHON="python"

# create output directories
mkdir -p "$LONGEST_ISOFORMS_DIR"
mkdir -p "$FAA_OUTPUT_DIR"

# create exclusion list from mapped_genes.tsv
grep -v '^#' "$MAPPED_GENES_FILE" | cut -f3 | tr ',' '\n' | sed '/^$/d' | sort -u > "$EXCLUDE_LIST"
echo "Exclusion list created: $EXCLUDE_LIST"

# process each unpacked gff3 file
echo "Processing GFF files to remove contaminants..."
for gff_file in "$ANNOTATIONS_DIR"/*.gff3; do
    species_name=$(basename "$gff_file" | sed 's/.gff3//')
    genome_file="$GENOME_FASTA_DIR/${species_name}/genome/genome.fa"

    # check if the genome file exists
    if [[ ! -f "$genome_file" ]]; then
        echo "ERROR: Genome file not found for species $species_name ($genome_file)"
        continue
    fi

    # filter the gff3 file to remove genes present in mapped_genes.tsv
    filtered_gff_file="${ANNOTATIONS_DIR}/${species_name}_filtered.gff3" 
    grep -Fv -f "$EXCLUDE_LIST" "$gff_file" > "$filtered_gff_file"
    echo "Filtered GFF3 for $species_name -> $filtered_gff_file"

    # Convert the filtered GFF3 to GTF 
    gtf_file="${ANNOTATIONS_DIR}/${species_name}.gtf"
    awk 'BEGIN {OFS="\t"} $1 !~ /^#/ {
        # Split the attributes field to extract gene_id and transcript_id
        split($9, attributes, ";");
        gene_id = "";
        transcript_id = "";
        for (i in attributes) {
            if (attributes[i] ~ /ID=/) {
                gene_id = substr(attributes[i], index(attributes[i], "=") + 1);
            } else if (attributes[i] ~ /Parent=/) {
                transcript_id = substr(attributes[i], index(attributes[i], "=") + 1);
            }
        }
        print $1, "source", $3, $4, $5, $6, $7, $8, "gene_id \"" gene_id "\"; transcript_id \"" transcript_id "\";"
    }' "$filtered_gff_file" > "$gtf_file"
    echo "Converted filtered gff3 to gtf for $species_name -> $gtf_file"

    # b. get the longest isoforms from the gtf file (can be useful here)
    #output_gtf="$LONGEST_ISOFORMS_DIR/${species_name}_longest_isoforms.gtf"
    #$PYTHON "$LONGEST_ISOFORM_SCRIPT" -g "$gtf_file" -o "$output_gtf"
    #echo "Longest isoforms extracted for $species_name -> $output_gtf"

    # c. convert gtf to faa using getAnnoFastaFromJoingenes.py (also possible to convert gff3 to faa)
    faa_base_output="$FAA_OUTPUT_DIR/${species_name}"  # Base name for output files
    $PYTHON "$GET_ANNO_FASTA_SCRIPT" -g "$genome_file" -o "$faa_base_output" --gtf "$output_gtf"
    echo "Converted to faa: ${faa_base_output}.aa"

    # optionally, remove the temporary filtered GFF file
    # rm "$filtered_gff_file"
done

# clean up temporary exclusion list
rm "$EXCLUDE_LIST"

echo "Processing complete. Faa files are in $FAA_OUTPUT_DIR."
# $FAA_OUTPUT_DIR directory can be used as an input for the next step.
```

## OrthoFinder analysis

The bash scripts and command to perform OrthoFinder analysis are described in [orthofinder.md](orthofinder.md).

## Ploidy estimation with smudgeplot  
Ploidy was estimated for diatoms with paired-end DNA sequences available on NCBI (35 total):
```
nextflow run smudgeplot.nf -profile singularity --sra "sra.txt"
```
**Scripts/Files:**
* [smudgeplot.nf](smudgeplot.nf) - Nextflow ploidy workflow (FETCH_SRA **->** KRAKEN **->** FASTK **->** SMUDGEPLOT)
* [nextflow.config](nextflow.config) - Nextflow config file; declares singularity profile and parameters
* sra.txt - List of diatom SRA accessions to be pulled from NCBI
```
$ head -5 sra.txt
SRR26112839
SRR14100021
SRR18733581
SRR18733586
SRR18733505
```
