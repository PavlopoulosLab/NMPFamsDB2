# NMPFamsDB2
We identified 608,258 previously unknown protein families (≥100 members) and 6.5 million families (≥25 members) from 40.3B sequences spanning metagenomes, metatranscriptomes, and reference genomes. Structural and genomic analyses reveal novel functions, greatly expanding the known protein universe.

# Protein Family Space Analysis Scripts

A collection of scripts and Unix command examples developed for the study **"Quadrupling the Protein Family Space with Global Metagenomics."** This repository contains helper scripts and command line workflows for filtering, searching, clustering, sequence alignment, trimming, and redundancy removal of protein sequence datasets.

## Repository Contents

### 1. Filtering_Searching_Clustering

A collection of Unix/Linux commands for:

* Filtering protein FASTA and tabular files from IMG/M
* Searching protein families against Pfam and reference genome databases
* Protein sequence clustering using MMSeqs

### 2. Scripts

#### `parser.py`
Splits large FASTA files into smaller chunks for efficient parallel processing.

#### `mafft_parallel.py`
Runs MAFFT multiple sequence alignments in parallel across multiple FASTA files.

#### `trimmer.py`
Trims multiple sequence alignments by removing poorly aligned regions and unnecessary columns.

#### `redundancy_removal.py`
Removes identical protein sequences from FASTA files.

## Requirements

* Python 3.8 or later
* Unix/Linux environment
* MAFFT
* MMSeqs2
* HHfilter
* tantan

## Citation

If you use this repository in your research, please cite:

> *"Quadrupling the Protein Family Space with Global Metagenomics."*

## License

See the `LICENSE` file for licensing terms.
