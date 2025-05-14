## The OrthoFinder program was run on 63 samples, including 9 previously annotated species (_Chaetoceros tenuissimus, Cylindrotheca closterium, Fragilaria crotonensis, Mayamaea pseudoterrestris, Nitzschia inconspicua, Phaeodactylum tricornutum, Pseudo-nitzschia multistriata, Seminavis robusta, Thalassiosira pseudonana_) and 5 species used as an ougroup (_Bremia lactucae, Phytophthora cinnamomi, Phytophthora infestans, Phytophthora ramorum, Phytophthora sojae_). 

As input data, you need to provide a folder with .faa files, for example Bacterosira_constricta.faa

Please note that all slurm script launch parameters should be adapted to your resources.

## Slurm script for the execution of OrthoFinder
```
#!/bin/bash
#SBATCH --job-name=OrthoFinder                 # job name
#SBATCH --output=orthofinder_%j.log    		   # std output and error log (%j: job ID)
#SBATCH --ntasks=1                             # run on a single CPU
#SBATCH --cpus-per-task=64                     # number of CPU cores per task
#SBATCH --mem=96000          
#SBATCH --time=72:00:00                        # time limit hrs:min:sec
#SBATCH --partition=snowball

source ~/.bashrc

# Define directories
#SPECIES=diatoms
BASE_DIR=/home/nenashen66/OrthoFinder/
CLEANED_DATA=/home/natalia/busco_input/cleaned_faa_files/
INPUT_DIR=${BASE_DIR}/input_diatoms/cleaned_faa_files/
OUTPUT_DIR=${BASE_DIR}/output_diatoms_cleaned

mv $CLEANED_DATA $INPUT_DIR # input data for the OrhoFinder run should be inside the $BASE_DIR

ulimit -n 10000
echo "Limit was changed:"
ulimit -Sn

/home/nenashen66/anaconda3/bin/python ~/OrthoFinder/orthofinder.py -f $INPUT_DIR -t 64 -a 16 -M msa -o $OUTPUT_DIR

# execution of the OrthoFinder from the step of generation of the orthogroups:
# also possible to use another method to generate a tree: -T fasttree
#/home/nenashen66/anaconda3/bin/python ~/OrthoFinder/orthofinder.py -fg $BASE_DIR -t 64 -a 16 -M msa -T fasttree
```

script submission: sbatch slurm_run.sh
The result of the program will be located in the directory: ${BASE_DIR}/output_diatoms_cleaned/Results_Apr07/

Statistical data with results are presented in the directory Comparative_Genomics_Statistic in the file  Statistics_PerSpecies.tsv
