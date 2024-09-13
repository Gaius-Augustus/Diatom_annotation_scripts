#!/bin/bash
#SBATCH --job-name=orp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G

module load singularity
module load nextflow

nextflow run agat.nf -profile singularity
