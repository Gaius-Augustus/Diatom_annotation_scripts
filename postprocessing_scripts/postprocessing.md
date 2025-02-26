# Processing of functionally decorated gff3 file

Contact: katharina.hoff@uni-greifswald.de

## Filtering single-exon genes

On the basis of functional annotation with EnTAP, results of an initial OrthoFinder run (with only the newly annotated protein sets), and a DIAMOND search against NCBI NR, the braker.gtf file was filtered as follows:

```
# Make intersection
overlapStat.pl single_exon_genes_without_hit.txt orthofinder_assigned_homolog.lst --outfiles
    
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

## OrthoFinder analysis

The bash scripts and command to perform OrthoFinder analysis are described in [orthofinder.md](orthofinder.md).
