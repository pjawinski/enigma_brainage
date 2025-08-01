#!/bin/bash
#SBATCH --job-name=wen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/phesant_wen.out
#SBATCH --error=logs/phesant_wen.err

# Load environment
cd /lustre/psychologie/jawinskp/projects/enigma_brainage
source ~/.bashrc                        # Make sure conda is available
conda activate envs/phesant

# Define variables
covsFile="results/phewas/sample/r2025.vars.discovery"
covs="sex,age,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
phenoFile="results/phewas/sample/ukb678162.phesant.2025.wba.mergedInstances.csv"
phesantDir="/lustre/psychologie/jawinskp/software/PHESANT"
nparts=25
partList="$(echo {1..25})"
standardise="TRUE"
genetic="FALSE"
sensitivity="FALSE"

# Run PHESANT for each trait
for trait in wen; do
    targetDir="results/phewas/mergedInstances/${trait}/"
    rm -rf "${targetDir}"
    mkdir -p "${targetDir}"
    traitFile="results/phewas/sample/pgs.sbayesrc.${trait}.score"
    code/phesant.sh "SCORE1_SUM" "${traitFile}" "${covsFile}" "${covs}" \
        "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${partList}" \
        "${standardise}" "${genetic}" "${sensitivity}" \
        2>&1 | tee "${targetDir}/phesant.log"
done
