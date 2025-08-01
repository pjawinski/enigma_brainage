[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Bluesky](https://img.shields.io/badge/Bluesky-pjawinski.bsky.social-blue?logo=bluesky)](https://bsky.app/profile/pjawinski.bsky.social)

# ENIGMA brainageFactor analysis
This repository contains the analysis scripts used to reproduce the UK Biobank-related results from our article "Cross-trait genetic architecture across six brain age models". The individual-level data used in this project were obtained from the [UK Biobank](https://www.ukbiobank.ac.uk/) under application number 42032. Details on sample selection and genetic data preprocessing were based on the workflow described in our previous UKB project, "Genome-wide analysis of brain age gap identifies 59 associated loci and unveils relationships with mental and physical health", available at https://github.com/pjawinski/ukb_brainage.

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

## Contact
Philippe Jawinski | Humboldt-Universität zu Berlin | philippe.jawinski[at]hu-berlin.de <br>
Sebastian Markett | Humboldt-Universität zu Berlin | sebastian.markett[at]hu-berlin.de
