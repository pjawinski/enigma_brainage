#!/bin/bash

# =====================
# === Calculate PGS ===
# =====================

# get arguments
weightFile="${1}" # weightFile="results/gap_gm/discovery/sbayesrc/sbayesrc.snpRes.gz"
fileHandler="${2}" # fileHandler='data/genetics/chr${i}/imp_mri_qc_EUR/chr${i}_mri_qc'
fileType="${3}" # fileType='pgen' | 'bgen'
idCol="${4}" # idCol="Name"
a1Col="${5}" # a1Col="A1"
effectCol="${6}" # effectCol="A1Effect"
maf="${7}" # maf=0.01
threads=${8} # threads=100
outFile="${9}" # outFile="data/genetics/prs/imp_mri_qc_EUR/gap_gm.sbayesrc"

# echo settings
echo $'\n'"--- Calculate PGS | Settings ---"
echo "weightFile: ${weightFile}"
echo "fileHandler: ${fileHandler}"
echo "fileType: ${fileType}"
echo "idCol: ${idCol}"
echo "a1Col: ${a1Col}"
echo "effectCol: ${effectCol}"
echo "maf: ${maf}"
echo "threads: ${threads}"
echo "outFile: ${outFile}"$'\n'

# set targetdir and make folder
targetDir="$(dirname "${outFile}")"
mkdir -p "${targetDir}"
mkdir -p "${outFile}.chr"
targetDir="$(readlink -f "${targetDir}")"

# creating weight file (filter variants with effect weight 0)
echo "Creating exposure file."
awk -v cols="${idCol},${a1Col},${effectCol}" 'BEGIN { ncols=split(cols,colnames,",") }
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colnums[i]=j } } } }
    $colnums[3]!=0 { print $colnums[1], $colnums[2], $colnums[3]}' <(zcat "${weightFile}") > "${outFile}".weights

# Calculate the number of parallel jobs (minimum 1 to avoid zero)
parallel_jobs=$(( threads / 10 ))
parallel_jobs=$(( parallel_jobs > 0 ? parallel_jobs : 1 ))

# Calculate threads per task
threads_per_task=$(( threads / parallel_jobs ))
threads_per_task=$(( threads_per_task > 0 ? threads_per_task : 1 ))

# calculate polygenic score
echo "Calculating scores for individual chromosomes"
> "${outfile}".jobs
if [ "$fileType" == "bgen" ]; then
    for i in {1..22}; do
        echo plink2 --bgen "$(eval echo "${fileHandler}.bgen")" --sample "$(eval echo "${fileHandler}.sample")" --maf "${maf}" --score "${outFile}".weights 1 2 3 header cols=+scoresums --threads "${threads_per_task}" --out "${outFile}.chr/chr${i}" >> "${outfile}".jobs
    done
elif [ "$fileType" == "pgen" ]; then
    for i in {1..22}; do
        echo plink2 --pfile "$(eval echo "${fileHandler}")" --maf "${maf}" --score "${outFile}".weights 1 2 3 header cols=+scoresums --threads "${threads_per_task}" --out "${outFile}.chr/chr${i}" >> "${outfile}".jobs
    done
fi

# Run the jobs in parallel with a maximum of 10 concurrent tasks
parallel --jobs ${parallel_jobs} --line-buffer < "${outfile}".jobs

# aggregate score
awk 'NR==1 { print "FID", "IID", "CHR_COUNT", "NMISS_ALLELE_CT", "SCORE1_SUM"; next }
     NR==FNR { id[$1]=$1; chrcount[$1]=1; nmiss[$1]=$3; score[$1]=$6; next }
     FNR==1 { next }
     { chrcount[$1]++; nmiss[$1]=nmiss[$1]+$3; score[$1]=score[$1]+$6 }
     END { for (i in score) { print id[i], id[i], chrcount[i], nmiss[i], score[i] } }' OFS='\t' "${outFile}.chr/"*.sscore \
     > ${outFile}.score

# clean up
echo "Cleaning up."
rm -f "${outFile}".weights
rm -f "${outfile}".jobs
tar -cvzf "${outFile}.chr.tar.gz" --directory="${targetDir}" "$(basename ${outFile}.chr)" --remove-files
chmod 770 "${outFile}"*
echo "--- Completed: Calculate PGS --- "


