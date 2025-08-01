#!/bin/bash

# =============================================
# === run phenome-wide association analysis ===
# =============================================

# get arguments
trait="${1}" # trait="gap_gm"
traitFile="${2}" # traitFile="data/${trait}/${trait}.txt"
covsFile="${3}" # covsFile="data/${trait}/covs.txt"
covs="${4}" # covs="sex,age,age2,ac1,ac2,TIV"
targetDir="${5}" # targetDir="results/${trait}/phesant"
phenoFile="$(readlink -f "$6")" # phenoFile="data/basket/20200409_2007685/data/ukb41573_imaging_wba_phesant.csv"
phesantDir="${7}" # phesantDir="/fast/software/PHESANT"
nparts=${8} # nparts=1
partList="${9}" # partList="1 6 11 16"
standardise="${10}" # standardise="FALSE"
genetic="${11}" # genetic="FALSE"
sensitivity="${12}" # sensitivity="FALSE"

# echo settings
echo $'\n'"--- PHESANT analysis Settings ---"
echo "trait: ${trait}"
echo "traitFile: ${traitFile}"
echo "covsFile: ${covsFile}"
echo "covs: ${covs}"
echo "targetDir: ${targetDir}"
echo "phenoFile: ${phenoFile}"
echo "phesantDir: ${phesantDir}"
echo "nparts: ${nparts}"
echo "partList: ${partList}"
echo "standardise: ${standardise}"
echo "genetic: ${genetic}"
echo "sensitivity: ${sensitivity}"$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# get exposure file
if [[ ! -f "${targetDir}/phesant.exposure.csv" ]]; then
	echo "Creating exposure file."
	awk -F '\t' -v trait="${trait}" 'BEGIN { print "eid,"trait; cols[1]="IID"; cols[2]=trait } 
		NR==1 { for(i=1;i<=length(cols);++i) { for(j=i;j<=NF;++j) { if($j==cols[i]) { colnums[i]=j } } } }
		NR>1 { output=$colnums[1]","$colnums[2]; print output}
		' OFS=',' "${traitFile}" > "${targetDir}"/phesant.exposure.csv
fi

# get confounder file
if [[ ! -f "${targetDir}/phesant.confounder.csv" ]]; then
	echo "Creating confounder file."
	header="eid,${covs}" 
	awk -F '\t' -v header="${header}" -v covs="${covs}" 'BEGIN { print header; ncovs=split(covs,covnames,",") }
		NR==1 { for(i=1;i<=ncovs;++i) { for(j=2;j<=NF;++j) { if($j==covnames[i]) { colnums[i]=j } } } }
		NR>1 { output=$1; for(i=1;i<=ncovs;++i) { output=output","$colnums[i] }; print output}' "${covsFile}" > "${targetDir}"/phesant.confounder.csv
fi

# run analysis
echo "Running phenome-wide scan."$'\n'
mkdir -p "${targetDir}"/phesant.output
cd "${phesantDir}"/WAS

if [[ ${sensitivity} == "TRUE" ]]; then
	for i in ${partList}; do (
	Rscript phenomeScan.r \
	--phenofile="${phenoFile}" \
	--traitofinterestfile="${targetDir}"/phesant.exposure.csv \
	--variablelistfile="${phesantDir}"/variable-info/outcome-info.tsv \
	--datacodingfile="${phesantDir}"/variable-info/data-coding-ordinal-info.txt \
	--traitofinterest="${trait}" \
	--resDir="${targetDir}"/phesant.output/ \
	--userId="eid" \
	--partIdx="${i}" \
	--numParts="${nparts}" \
	--standardise="${standardise}" \
	--genetic="${genetic}" \
	--sensitivity
	) &
	done
else 
	for i in ${partList}; do (
	Rscript phenomeScan.r \
	--phenofile="${phenoFile}" \
	--traitofinterestfile="${targetDir}"/phesant.exposure.csv \
	--variablelistfile="${phesantDir}"/variable-info/outcome-info.tsv \
	--datacodingfile="${phesantDir}"/variable-info/data-coding-ordinal-info.txt \
	--traitofinterest="${trait}" \
	--resDir="${targetDir}"/phesant.output/ \
	--userId="eid" \
	--partIdx="${i}" \
	--numParts="${nparts}" \
	--confounderfile="${targetDir}"/phesant.confounder.csv \
	--standardise="${standardise}" \
	) &
	done
fi
wait
echo "Finished."

# post phenome scan results processing
numPartsProcessed=$(wc -w <<< "${partList}")
if [[ "${numPartsProcessed}" -eq "${nparts}" ]]; then
	cd "${phesantDir}"/resultsProcessing/
	Rscript mainCombineResults.r \
	--resDir="${targetDir}"/phesant.output/ \
	--variablelistfile="${phesantDir}"/variable-info/outcome-info.tsv \
	--numParts="${nparts}"
else
	echo "Skipping result combination: only ${numPartsProcessed} out of ${nparts} parts were run."
fi

# create visualization
# mkdir -p "${targetDir}/results/web"
# cd ${phesantDir}/PHESANT-viz/bin/
# java -cp .:../jar/json-simple-1.1\ 2.jar ResultsToJSON "${targetDir}/results/results-combined.txt" "../node-positions.csv" "${targetDir}/results/web/java-json.json"

# set permissions
chmod -R 770 "${targetDir}"

