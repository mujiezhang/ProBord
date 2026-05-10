<div align="center">
  <img src="https://github.com/user-attachments/assets/1e6a744b-1286-49a9-89b0-fe815ea30a35" alt="ProBord" width="450" />
</div>

# ProBord: **Pro**virus **Bord**er Delimiter ✨
ProBord (**Pro**virus **Bord**er Delimiter) is a bioinformatics tool that predicts the precise borders of proviruses by identifying attL/R sites.

## Table of contents
<!-- TOC -->
- [ProBord: Provirus Border Delimiter](#probord-provirus-border-delimiter-)
- [Introduction](#introduction)
  - [Provirus integration process](#-provirus-integration-process)
  - [Workflow of ProBord](#-workflow-of-probord)
- [Instructions](#instructions)
  - [Dependencies](#dependencies)
  - [**Installation**](#installation)
  - [**Database preparation**](#database-preparation)
  - [**How to run**](#how-to-run)
  - [**Output files**](#output-files)
- [Citation](#citation)
- [Contact](#-contact)

<!-- /TOC -->

---

# Introduction
## 🧬 Provirus integration process

A provirus usually refers to a virus integrated into a prokaryotic chromosome as a stable genetic element. Before integration, the phage attP site and the host attB site—share core sequence—undergo site-specific recombination catalyzed by integrase (Int), producing attL and attR sites flanking the prophage in the host genome. During excision, attL and attR recombine in reverse, mediated by integrase and excisionase (Xis), restoring attP on the free phage DNA and attB on the host chromosome. Unless otherwise stated, attB, attP, and attL/R refer to their core sequences.

<img width="1441" height="355" alt="integration" src="https://github.com/user-attachments/assets/6526af29-8f12-4e69-ac0a-8e521f33c703" />

## 💡 Workflow of ProBord

- Step1: Preprocessing viral region
  - The inputs to ProBord are proviruses predicted using external tools such as geNomad, and host contamination is removed from the predicted viral regions using CheckV, and the trimmed provirus sequence is extended by 5 kb on both sides to generate a “host–attL–provirus–attR–host” mixed sequence (mix-seq) that captures the potential integration neighborhood;
- Step2: Identifying candidate att cores (CACs) using length-dependent strategies
  - For short CACs (5–11 bp), mix-seq is aligned against prokaryotic reference genomes, and att-hot regions are located based on cumulative base coverage, within which short CACs are scanned.
  - For long CACs (≥12 bp), BLASTn is used to align the two 25 kb terminal regions of mix-seq to identify matching terminal CACs.
- Step3: Candidate attB comparing and scoring
  - The left and right CACs are extended by 100 bp into the flanking host regions and assembled into a candidate attB, which is then aligned against prokaryotic reference genomes and scored. The highest-scoring candidate attB is used to trace back the corresponding attL/attR positions, thereby inferring the precise integration boundary of the provirus.

<img width="1782" height="1104" alt="probord-wokflow" src="https://github.com/user-attachments/assets/ce23bd73-c60e-403d-a052-47f174082ef0" />


# Instructions

## Dependencies
- ProBord is a Python script that relies on:
  
```
blastn
python3
biopython
checkv=1.0.3
ncbi-genome-download
```

## Installation

- Install miniconda and add channels (**If already installed, please skip**)
```
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels bioconda
conda config --add channels conda-forge
```
- Install ProBord
```
conda create -n probord probord
conda activate probord
```

## Database preparation
- Prepare the CheckV database (if needed; otherwise skip):  `checkv download_database ./ `
- Prepare blastn database for attB detection (**required**):
  - (**✅ recommended**) If your provirus originates from a specific bacterial/archaeal genus, you only need to download bacterial/archaeal genomes and create a blastn database for that genus using the script `prepare_blastn_db.sh`. For example, for the genus "Mannheimia": `bash prepare_blastn_db.sh Mannheimia bacteria`.
  - If you have numerous proviruses from diverse genera, or if you don't know your provirus host classification, you can download the NCBI nt database ( https://ftp.ncbi.nlm.nih.gov/blast/db/ ) or all bacterial/archaeal genomes from NCBI RefSeq (bacteria: https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/ , archaea: https://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/ ), create a blastn database, and then run probord. (This approach consumes substantial storage space and memory, and will significantly increase probord's runtime.)
    
    **Note**: We are currently developing algorithms to compress DNA sequences while preserving potential attB sites, aiming to reduce runtime memory consumption.

## How to run
- ▶️ Command line options: `probord -h`:
```
usage: probord [-h] [-v] {run,prepare_db} ...

ProBord: Provirus Border Delimiter v1.0

positional arguments:
  {run,prepare_db}  Available commands
    run             Run the ProBord pipeline for provirus border delimitation
    prepare_db      Download NCBI genomes and build BLAST database for a genus

options:
  -h, --help        show this help message and exit
  -v, --version     show program's version number and exit
```

We provide two test datasets:

🚩 `Mannheimia phage vB_MhM_3927AP2` and its host contig: `NZ_CP017531.1.fna`: This transposable phage features exceptionally short attL/R sites (5 bp).

🚩 `Haemophilus phage HP2` and its host contig: `LR134490.1.fna`: This phage contains long attL/R sites (182 bp).

These datasets respectively represent: `Short-att phages (5–11 bp att sites)` and `Long-att phages (≥12 bp att sites)`
- run an example
```
# prepare blastn db for genera Mannheimia
probord prepare_db Mannheimia bacteria

# run ProBord with default parameters
probord run -hf test/data/phage_vB_MhM_3927AP2/NZ_CP017531.1.fna  -vf test/data/phage_vB_MhM_3927AP2/phage_vB_MhM_3927AP2_location.tsv -wd phage_vB_MhM_3927AP2_prediction -cv checkv-db-v1.5/ -db Mannheimia/Mannheimia
```

## Output files
In this example, the results of ProBord's analysis will be written to the `phage_vB_MhM_3927AP2_prediction` directory, which will look like this:
```
phage_vB_MhM_3927AP2_prediction/
├── attB_blast
│   └── attB_mix_outfmt.txt
├── att_prediction.tsv
└── probord.log
```
1. `attB_mix_outfmt.txt`: blast result of all candidate attB
2. `att_prediction.tsv`: the prediction result
3. `probord.log`: log file

A detailed overview of `att_prediction.tsv`:
| original_name | host_contig | original_start | original_end | temp_name | attL_start | attL_end | attR_start | attR_end | att_length | attL_sequence | attR_sequence | attB_score |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Mannheimia phage vB_MhM_3927AP2 | NZ_CP017531.1 | 829886 | 863606 | Mannheimia--phage--vB_MhM_3927AP2__NZ_CP017531.1__829886__863606_829886-829891-863601-863606-100-100-6 | 829886 | 829891 | 863601 | 863606 | 6 | AATACT | AATACT | 100.0 |



# Citation
......

# 📬 Contact
```
# Mujie Zhang
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University
# Email: zhangmujie@sjtu.edu.cn
```
