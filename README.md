# MitoCatch

An automated, reference-independent bioinformatics pipeline for the recovery, assembly, and annotation of complete mitochondrial genomes from NGS data. Developed in Snakemake to ensure efficiency, scalability, and reproducibility in evolutionary and conservation studies.

## Prerequisites
* **Conda** or **Mamba** installed.
* Linux-based environment (or WSL on Windows).

## Installation and Setup

### 1. Clone the repository
```bash
git clone [https://github.com/evochange/MitoCatch.git](https://github.com/evochange/MitoCatch.git)
cd MitoCatch
```

### 2. Environment Setup
You can set up the required environments using either `conda` or `mamba`.

#### Option A: Using Conda

```bash
# Main Pipeline Environment
conda create -n mitocatch_env -c conda-forge -c bioconda -c defaults snakemake=9.23.0 biopython trimmomatic sra-tools spades magic-wormhole blast -y

# QC and Pre-processing Environments
conda create -n mitocatch_qc -c bioconda -c conda-forge -c defaults fastqc=0.12.1 multiqc=0.9.1a0 -y
conda create -n mitocatch_trimmomatic -c bioconda -c conda-forge -c defaults trimmomatic=0.39 -y
conda create -n mitocatch_seqtk -c bioconda -c conda-forge -c defaults seqtk=1.5 -y

# Assembly Environments
conda create -n mitocatch_getorganelle -c conda-forge -c bioconda python=3.9 getorganelle=1.7.7.1 blast=2.5.0 spades=4.2.0 bowtie2=2.5.5 -y
conda create -n mitocatch_trinity -c bioconda -c conda-forge -c defaults trinity=2.1.1 samtools=1.21 htslib=1.21 jellyfish=2.2.10 -y

# Annotation and Validation Environment
conda create -n mitocatch_annotation -c bioconda -c conda-forge -c defaults blast=2.17.0 biopython=1.85 entrez-direct=24.0 -y
```

#### Option B: Using Mamba (Recommended for speed)
```bash
conda install -n base -c conda-forge mamba

# Run the following commands substituting 'conda' with 'mamba'
mamba create -n mitocatch_env -c conda-forge -c bioconda snakemake=9.23.0 biopython trimmomatic sra-tools spades magic-wormhole blast -y
mamba create -n mitocatch_qc -c bioconda -c conda-forge -c defaults fastqc=0.12.1 multiqc=0.9.1a0 -y
mamba create -n mitocatch_trimmomatic -c bioconda -c conda-forge -c defaults trimmomatic=0.39 -y
mamba create -n mitocatch_seqtk -c bioconda -c conda-forge -c defaults seqtk=1.5 -y
mamba create -n mitocatch_getorganelle -c conda-forge -c bioconda python=3.9 getorganelle=1.7.7.1 blast=2.5.0 spades=4.2.0 bowtie2=2.5.5 -y
mamba create -n mitocatch_trinity -c conda-forge -c bioconda trinity=2.1.1 samtools=1.21 htslib=1.21 jellyfish=2.2.10 -y
mamba create -n mitocatch_annotation -c conda-forge -c bioconda blast=2.17.0 biopython=1.85 entrez-direct=24.0 -y
```

### 3. Initialize GetOrganelle Database
Required only upon first installation.

```bash
conda activate mitocatch_getorganelle
get_organelle_config.py -a all
conda deactivate
```

## How to Run the Pipeline

### 1. Activate the main environment:
```bash
conda activate mitocatch_env
```

### 2. Execute the pipeline:
```bash
bash mitocatch.sh --r1 path/to/read1.fastq.gz --r2 path/to/read2.fastq.gz --modo getorganelle --organelle_type animal_mt
```

## Usage Notes:

* --r1 / --r2: Replace path/to/read1.fastq.gz with the actual path and filename of your input FASTQ files.
* --modo: Supports getorganelle (recommended) or trinity.
* --organelle_type: Required only for getorganelle mode. Use the appropriate nomenclature for your sample (e.g., animal_mt, plant_mt, fungus_mt).

## Test Dataset
To verify that the pipeline is correctly installed and functioning, a test dataset is available for download from the NCBI SRA database:

* **Data Source**: [MyxoHares - Unravelling the genomic landscape of myxomatosis susceptibility in Iberian hares](https://www.ncbi.nlm.nih.gov/sra/ERX15431688[accn])
* **Accession Number**: ERR16040998

### How to run the test
#### 1. Prepare the directory:
Ensure you are inside your project's root folder and navigate to the data directory:
```bash
mkdir -p data/raw #If you don't have it already
cd data/raw
```

#### 2. Download the data:
Activate your `mitocatch_env` and download the FASTQ files using `sra-tools`:
```bash
conda activate mitocatch_env
fastq-dump --split-files --gzip ERR16040998
```

#### 3. Execute the pipeline:
Return to the project root and run the pipeline. You can choose between `getorganelle` or `trinity` modes:
* Using GetOrganelle mode:
```bash
cd ../..
bash mitocatch.sh --r1 data/raw/ERR16040998_1.fastq.gz --r2 data/raw/ERR16040998_2.fastq.gz --modo getorganelle --organelle_type animal_mt
```
* Using Trinity mode:
```bash
cd ../..
bash mitocatch.sh --r1 data/raw/ERR16040998_1.fastq.gz --r2 data/raw/ERR16040998_2.fastq.gz --modo trinity 
```

## Video Tutorial
To facilitate understanding of the workflow and demonstrate the practical execution of the pipeline, a video tutorial is available:

* **Watch online:** [MitoCatch - How to use (Tutorial)](https://youtu.be/RVB9Mt4O4sg?si=TTJXA-aC7H9JL7JL)
* **Local download:** The video file (`MitoCatch - How to use.mp4`) is also available directly within the repository for offline consultation.

