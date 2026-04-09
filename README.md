# Bioinformatic Processing of 16S rRNA Amplicon and Long-Read Metagenomics Datasets from Freshwater Larval Gut Bacterial Communities

## Overview
This repository contains the code used to process and analyse 16S rRNA amplicon and long-read metagenomics datasets derived from freshwater larval gut bacterial communities.

## Repository Contents

| File | Description |
|------|-------------|
| `larvae_exp_script_bash.Rmd` | R Markdown file containing bash code chunks for preprocessing and pipeline steps, executed via SLURM on an HPC cluster |
| `16S_analyses_script_BLIK_III.R` | R script to reproduce the full analysis of the 16S rRNA amplicon dataset |

## Requirements
- R (≥ 4.3.1)
- SLURM-based HPC cluster (for bash pipeline)
- DADA2, phyloseq, etc.

## Usage
1. Clone the repository:
```bash
   git clone https://github.com/your-username/your-repo-name.git
```
2. Run the bash pipeline by submitting the chunks in `larvae_exp_script_bash.Rmd` via SLURM.
3. Reproduce the 16S rRNA amplicon analyses by running `16S_analyses_script_BLIK_III.R` in R.
