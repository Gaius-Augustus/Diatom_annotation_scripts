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

## Retrieving longest isoform for final OrthoFinder analysis (from braker.aa)

```
get_longest_isoform_from_braker_aa.py -i braker_filtered.aa -o braker_longest.aa
```

The contents of species-specific fixed_dbxref.gff3 files are provided at Zenodo with doi 10.5281/zenodo.13745090

## OrthoFinder analysis

The bash scripts and command to perform OrthoFinder analysis are described in (orthofinder.md)[orthofinder.md].
