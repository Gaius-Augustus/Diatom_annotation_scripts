#!/bin/bash

output_file="aggregated_stats.tsv"

# Initialize the output file with headers
echo -e "Species\tTotal_Input_Sequences\tUnique_with_Alignment\tFlagged_Sim_Search_Contaminant\tNot_Flagged_Sim_Search_Contaminant\tUnique_without_Alignment\tUnique_with_Family_Assignment\tUnique_without_Family_Assignment\tUnique_with_GO_Term\tUnique_with_Pathway_Assignment\tFlagged_EggNOG_Contaminant\tNot_Flagged_EggNOG_Contaminant\tAnnotated_Similarity_Search_Only\tAnnotated_Gene_Family_Only\tAnnotated_Gene_Family_and/or_Similarity_Search\tAnnotated_Contaminant_Gene_Family_and/or_Similarity_Search\tAnnotated_Without_Contaminant_Gene_Family_and/or_Similarity_Search\tUnique_Unannotated\tTotal_Runtime_Minutes" > "$output_file"

# Loop through directories
for dir in $(ls -dr shared_with_uconn/*/); do
    dir=${dir%*/}      # remove the trailing "/"
    species="${dir##*/}"
    
    # Find the most recent log file
    out_dir=$(ls -t "${dir}/braker/entap_outfiles/"log*.txt | head -n 1)
    
    if [[ -f "$out_dir" ]]; then
        total_input=$(grep -m 1 "Total Input Sequences:" "$out_dir" | awk -F': ' '{print $2}')
        alignments=$(grep -m 1 "Total unique sequences with an alignment:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        flagged_sim_cont=$(grep -m 1 "Total alignments flagged as a Similarity Search contaminant:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        not_flagged_sim_cont=$(grep -m 1 "Total alignments NOT flagged as a Similarity Search contaminant:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        no_align=$(grep -m 1 "Total unique sequences without an alignment:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        fam_assigned=$(grep -m 1 "Total unique sequences with family assignment:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        fam_not_assigned=$(grep -m 1 "Total unique sequences without family assignment:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        go_term=$(grep -m 1 "Total unique sequences with at least one GO term:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        pathway_assigned=$(grep -m 1 "Total unique sequences with at least one pathway (KEGG) assignment:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        flagged_eggnog_cont=$(grep -m 1 "Total unique sequences flagged as an EggNOG contaminant:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        not_flagged_eggnog_cont=$(grep -m 1 "Total unique sequences NOT flagged as an EggNOG contaminant:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
      #  retained=$(grep -m 1 "Total retained sequences (after filtering and/or frame selection):" "$out_dir" | awk -F': ' '{print $2}')
        annotated_sim_search=$(grep -m 1 "Total unique sequences annotated (similarity search alignments only):" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        annotated_gene_fam=$(grep -m 1 "Total unique sequences annotated (gene family assignment only):" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        annotated_both=$(grep -m 1 "Total unique sequences annotated (gene family and/or similarity search):" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        annotated_cont=$(grep -m 1 "Total annotated sequences flagged as a contaminant from either Similarity Search or EggNOG:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        annotated_no_cont=$(grep -m 1 "Total annotated sequences NOT flagged as a contaminant:" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        unannotated=$(grep -m 1 "Total unique sequences unannotated (gene family and/or similarity search):" "$out_dir" | awk -F': ' '{print $2}' | awk '{print $1}')
        runtime=$(grep -m 1 "Total runtime (minutes):" "$out_dir" | awk -F': ' '{print $2}')

        # Append the results to the output file
        echo -e "$species\t$total_input\t$alignments\t$flagged_sim_cont\t$not_flagged_sim_cont\t$no_align\t$fam_assigned\t$fam_not_assigned\t$go_term\t$pathway_assigned\t$flagged_eggnog_cont\t$not_flagged_eggnog_cont\t$annotated_sim_search\t$annotated_gene_fam\t$annotated_both\t$annotated_cont\t$annotated_no_cont\t$unannotated\t$runtime" >> "$output_file"
    else
        echo "No .txt log file found in $dir"
    fi
done

sed -i 's/_/ /g' "$output_file"
