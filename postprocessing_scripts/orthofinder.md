# The OrthoFinder program was run on 59 samples in two steps: an initial run on 40 samples (selected in alphabetical order) and on the remaining 19. 

In the second step, we use an additional parameter that allows us to add samples to an already completed analysis.
As input data, you need to provide a folder with .faa files, for example Bacterosira_constricta.faa

Please note that all slurm script launch parameters should be adapted to your resources.

## 1. Slurm script for the first run:
```
#!/bin/bash
#SBATCH --job-name=OrthoFinderInitial
#SBATCH --output=orthofinder_initial_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96000          
#SBATCH --time=48:00:00
#SBATCH --partition=snowball
# Define directories
BASE_DIR=~/diatoms
FIRST_SUBSET_DIR=${BASE_DIR}/input_set1
INPUTPUT_DIR=${BASE_DIR}/OF_initial
# Run OrthoFinder for the first subset
OrthoFinder/orthofinder.py -f $FIRST_SUBSET_DIR -t 48 -o $OUTPUT_DIR
```
script submission: sbatch slurm_run1.sh
The result of the program will be located in the directory: ${BASE_DIR}/OF_initial/Results_Sep12/

## 2. Slurm script for the second run:
```
#!/bin/bash
#SBATCH --job-name=OrthoFinderAddSubset
#SBATCH --output=orthofinder_add_%A_%a.log
#SBATCH --ntasks=1                        
#SBATCH --partition=snowball
BASE_DIR=~/diatoms
PREVIOUS_OUTPUT_DIR=${BASE_DIR}/OF_initial/Results_Sep12/
SUBSET_DIR=$1
# before the run of orthofinder, check and change the limit for the number of opened files, if it is necessary
ulimit -n 10000 
OrthoFinder/orthofinder.py -b $PREVIOUS_OUTPUT_DIR/WorkingDirectory/ -f $SUBSET_DIR -t 48
```
script submission: sbatch slurm_run2.sh ~/diatoms/input_set2

The result of the program will be located in the directory: ${BASE_DIR}/OF_initial/Results_Sep12/WorkingDirectory/OrthoFinder/Results_Sep12/
Statistical data with results are presented in the directory Comparative_Genomics_Statistic in the file  Statistics_PerSpecies.tsv
