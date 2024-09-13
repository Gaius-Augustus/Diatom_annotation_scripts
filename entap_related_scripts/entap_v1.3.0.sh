#!/bin/bash
#SBATCH --job-name=entap_ecar_run
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=70G
#SBATCH --qos=general
#SBATCH -o ecar_%j.out
#SBATCH -e ecar_%j.err
#SBATCH --partition=general

unset PYTHONPATH
unset PYTHONHOME

module load eggnog-mapper/2.1.7
module load python/3.8.1
module load sqlite/3.31.1
module load perl/5.30.1
module load interproscan/5.25-64.0
module load TransDecoder/5.3.0
module load rsem/1.3.0
module load singularity/3.9.2

unset PYTHONPATH
unset PYTHONHOME

#normal run
for dir in shared_with_uconn/*/     # list directories in the form "/tmp/dirname/"
do
    dir=${dir%*/}      # remove the trailing "/"
    species="${dir##*/}"
    fasta_file="${dir}/braker/braker.aa"
    out_dir="${dir}/braker/entap_outfiles"
    entap_cmd="/core/labs/Wegrzyn/EnTAP/EnTAP_v1.3.0/EnTAP/EnTAP --runP --entap-ini /home/FCAM/ahart/EnTAP/paper_collab/entap_v1.3.0_config.ini --run-ini /home/FCAM/ahart/EnTAP/paper_collab/entap_v1.3.0_params --out-dir ${out_dir} --input ${fasta_file} --taxon ${species}"
    echo "${entap_cmd}"
    ${entap_cmd}
done
