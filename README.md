[![Codacy Badge](https://app.codacy.com/project/badge/Grade/5b75cbd059774762b4ac2a184e37386d)](https://app.codacy.com/gh/pjawinski/enigma_brainage/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![ENIGMA | Brain Age](https://img.shields.io/badge/ENIGMA-Brain%20Age-brightgreen)](https://enigma.ini.usc.edu/ongoing/enigma-brain-age/)
[![Bluesky](https://img.shields.io/badge/Bluesky-pjawinski.bsky.social-blue?logo=bluesky)](https://bsky.app/profile/pjawinski.bsky.social)

# Genetic architecture of brain age gap
This repository contains the analysis scripts to reproduce the UK Biobank-related results from our article **"Shared genetic architecture of brain age gap across 30 cohorts worldwide"**. The individual-level data used in this project were obtained from the [UK Biobank](https://www.ukbiobank.ac.uk/) under application number 42032. Details on UK Biobank sample selection and genetic data preprocessing follow the workflow described in our [previous article](https://doi.org/10.1038/s43587-025-00962-7), with the corresponding code available [here](https://github.com/pjawinski/ukb_brainage). This work contributes to the overarching multi-cohort ENIGMA genetic analyses available in [this related repository](https://github.com/VilteBaltra/genetic-architecture-of-brain-age-gap).

## Folder structure
[code/](code/) - preparation scripts, helper functions, and analysis workflows  <br>
[envs/](envs/) - conda `.yml` files to recreate the software environments  <br>
[results/](results/) - output files (no individual-level results are shared due to UK Biobank data privacy policies)<br>
[run.sh](run.sh) - main analysis script for brain age prediction, GWAS of brain age, polygenic score (PGS) analyses, and PheWAS <br>

## Software Environment
Analyses were run on **Debian GNU/Linux 11 (bullseye)** with  **kernel version 5.10.0-23-amd64**. The [code/prepare/](code/prepare/) directory contains scripts to facilitate the installation of the necessary bioinformatic tools for reproducing our analyses. For managing conda environments, we recommend using [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), which offers faster dependency resolution and package installation compared to `conda`.

**Required Tools:**
Below is a list of the primary tools utilized in our analysis, along with their respective versions and roles:

- **[R](https://www.r-project.org/)** `4.1.1 and 4.4.2` | Statistical computing and plotting, included in conda environments
- **[PHOTON-AI](https://photon-ai.com/enigma_brainage)** `ENIGMA brainage` singularity container | Brain age estimation from MRI data
- **[BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/downloads/)** `v2.4.1` | Linear mixed-model GWAS for brain age gap
- **[PHESANT](https://github.com/MRCIEU/PHESANT)** `v1.1` | Phenome-wide association analysis (PheWAS) in UK Biobank
- **[PLINK 1.9](https://www.cog-genomics.org/plink/)** `v1.90b6.8` 64-bit | Genomic preprocessing
- **[PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)** `v2.00a2LM` 64-bit Intel | Genomic preprocessing
- **[GCTB](https://cnsgenomics.com/software/gctb/)** `v2.5.2` | Polygenic score weighting using SBayesRC

## Cloning the Repository
Navigate to your preferred local directory and clone this repository via the following commands:
```
git clone https://github.com/pjawinski/enigma_brainage
cd enigma_brainage
```

- For genetics file preparation, follow the scripts of our previous brain age project [code/prepare.genetics.sh](https://github.com/pjawinski/ukb_brainage/blob/main/code/prepare.genetics.sh)
- For a step-by-step overview of the analysis workflow, refer to the main analysis script: [run.sh](run.sh)
> **Note:** Direct downloads of genetic data are deprecated. All data access and preprocessing must now be performed using the [UK Biobank Research Analysis Platform (RAP)](https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform).

## Contact
Philippe Jawinski | Humboldt-Universität zu Berlin, Germany | philippe.jawinski[at]hu-berlin.de <br>
Vilte Baltramonaityte | University of Bath, United Kingdom | vb506[at]bath.ac.uk
