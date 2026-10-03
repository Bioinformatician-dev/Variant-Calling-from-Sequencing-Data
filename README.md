# 🧬 Variant Calling from Sequencing Data

A bioinformatics workflow for identifying **genetic variants from next-generation sequencing (NGS) data**.

The project demonstrates the computational steps involved in transforming sequencing reads into a set of candidate genomic variants, providing a foundation for downstream genomic and variant interpretation analyses.

---

## 🔬 Overview

Variant calling is a fundamental task in genomics used to identify differences between an individual's sequenced genome and a reference genome.

These differences may include:

* Single-nucleotide variants (SNVs)
* Single-nucleotide polymorphisms (SNPs)
* Insertions and deletions (INDELs)
* Other small genomic variants

A typical variant-calling workflow is:

```text
        Raw Sequencing Reads
                 │
                 ▼
        Quality Assessment
                 │
                 ▼
           Read Trimming
                 │
                 ▼
       Reference Alignment
                 │
                 ▼
       BAM Processing/Sorting
                 │
                 ▼
          Variant Calling
                 │
                 ▼
        Variant Filtering
                 │
                 ▼
        Variant Annotation
                 │
                 ▼
          Final VCF File
```

---

## ✨ Key Features

The project is designed around the major stages of an NGS variant-calling workflow:

* 🧬 Sequencing-read processing
* 🔬 Quality control
* 🧭 Reference-genome alignment
* 📦 BAM/SAM processing
* 🧬 Variant identification
* 🔎 Variant filtering
* 📝 VCF generation
* 🧪 Downstream variant annotation

---

## 🧬 What Is Variant Calling?

Variant calling is the process of identifying genomic differences from sequencing data.

For example:

```text
Reference:
ATGCGTACGATCGATCG

Sample:
ATGCGTTCGATCGATCG
       ↑
      Variant
```

The observed difference can be represented in a standardized variant format such as VCF.

---

## 🧪 Workflow

### 1. Raw Sequencing Data

The analysis begins with sequencing reads, commonly provided as FASTQ files.

Example:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

For paired-end sequencing:

```text
R1 → Forward reads
R2 → Reverse reads
```

---

### 2. Quality Control

Raw reads should first be assessed for sequencing quality.

Typical quality-control metrics include:

* Per-base sequence quality
* Per-sequence quality
* Adapter contamination
* Sequence duplication
* GC-content
* Read-length distribution

A commonly used tool is **FastQC**.

---

### 3. Read Trimming

Low-quality bases and sequencing adapters can be removed before alignment.

Possible tools include:

* Fastp
* Trimmomatic
* Cutadapt

Example:

```text
Raw FASTQ
    │
    ▼
Adapter / Quality Filtering
    │
    ▼
Clean FASTQ
```

---

### 4. Reference Genome Alignment

Clean reads are aligned against a reference genome.

For short-read DNA sequencing, **BWA-MEM** is a commonly used aligner.

Conceptually:

```text
Clean Reads
     │
     ▼
Reference Genome
     │
     ▼
Aligned Reads
     │
     ▼
SAM/BAM
```

---

### 5. BAM Processing

The resulting alignment files are processed and sorted.

Typical operations include:

```text
SAM
 │
 ▼
BAM
 │
 ▼
Sort
 │
 ▼
Index
```

Tools such as **SAMtools** can be used for these operations.

---

### 6. Variant Calling

Variant callers examine aligned reads and identify positions where the sample differs from the reference.

Possible variant-calling frameworks include:

* GATK
* bcftools
* FreeBayes

The output is commonly stored as a **VCF (Variant Call Format)** file.

---

## 📄 VCF Output

A simplified VCF record may look like:

```text
#CHROM POS ID REF ALT QUAL FILTER INFO
chr1   12345 .  A   G   99   PASS   DP=40
```

Where:

| Field  | Meaning                        |
| ------ | ------------------------------ |
| CHROM  | Chromosome/contig              |
| POS    | Genomic position               |
| ID     | Variant identifier             |
| REF    | Reference allele               |
| ALT    | Alternative allele             |
| QUAL   | Variant quality                |
| FILTER | Filtering status               |
| INFO   | Additional variant information |

---

## 🔎 Variant Filtering

Raw variant calls may contain low-confidence variants.

Filtering can consider metrics such as:

* Read depth
* Variant quality
* Allele balance
* Mapping quality
* Strand bias
* Genotype quality

A simplified workflow is:

```text
Raw Variants
     │
     ▼
Quality Filters
     │
     ▼
High-Confidence Variants
```

> Filtering thresholds should be defined according to the sequencing platform, experimental design, variant caller, and biological application rather than using arbitrary universal cutoffs.

---

## 🧬 Variant Annotation

After high-confidence variants are identified, they can be annotated to determine their potential biological context.

Annotation may include:

* Gene affected
* Transcript
* Coding consequence
* Amino-acid change
* Known database identifiers
* Population frequency
* Clinical significance, where applicable

Possible annotation resources/tools include:

* ANNOVAR
* SnpEff
* VEP
* ClinVar

---

## 🛠️ Technologies & Tools

| Tool                | Purpose                          |
| ------------------- | -------------------------------- |
| FastQC              | Read quality assessment          |
| Fastp / Trimmomatic | Read preprocessing               |
| BWA-MEM             | Reference alignment              |
| SAMtools            | BAM/SAM processing               |
| GATK                | Variant discovery and processing |
| BCFtools            | Variant manipulation             |
| SnpEff / VEP        | Variant annotation               |
| Python              | Workflow automation              |
| Linux               | Computational environment        |

> The exact tools used by the current implementation should be documented here as the pipeline develops.

---

## 📂 Repository Structure

```text
Variant-Calling-from-Sequencing-Data/
│
├── README.md
├── main.py
│
├── data/
│   ├── raw/
│   └── processed/
│
├── reference/
│
├── results/
│   ├── bam/
│   ├── vcf/
│   └── annotation/
│
└── logs/
```

---

## 🚀 Installation

### Clone the repository

```bash
git clone https://github.com/Bioinformatician-dev/Variant-Calling-from-Sequencing-Data.git
cd Variant-Calling-from-Sequencing-Data
```

### Create a Conda environment

```bash
conda create -n variant_calling python=3.10
conda activate variant_calling
```

Bioinformatics dependencies can then be installed through Conda/Bioconda.

Example:

```bash
conda install -c bioconda fastqc fastp bwa samtools bcftools
```

Additional tools such as GATK and annotation software can be installed according to the workflow configuration.

---

## ▶️ Usage

Run the main Python workflow:

```bash
python main.py
```

A complete command-line implementation can eventually support parameters such as:

```bash
python main.py \
    --reads sample_R1.fastq.gz sample_R2.fastq.gz \
    --reference reference.fasta \
    --output results/
```

---

## 📊 Expected Outputs

A complete workflow can generate:

```text
results/
│
├── qc/
│   ├── fastqc/
│   └── multiqc/
│
├── trimmed/
│
├── alignment/
│   ├── sample.sorted.bam
│   └── sample.sorted.bam.bai
│
├── variants/
│   ├── raw.vcf.gz
│   └── filtered.vcf.gz
│
└── annotation/
    └── annotated_variants.tsv
```

---

## 🧪 Quality-Control Strategy

Quality control should be performed at multiple stages:

```text
Raw Reads
   │
   ├── FastQC
   │
   ▼
Trimmed Reads
   │
   ├── FastQC
   │
   ▼
Aligned Reads
   │
   ├── Alignment statistics
   ├── Mapping quality
   └── Coverage
   │
   ▼
Variant Calls
   │
   ├── Depth
   ├── Quality
   ├── Allele balance
   └── Filtering
```

This helps identify technical problems before interpreting biological results.

---

## 🎯 Applications

Variant-calling workflows are widely used in:

* 🧬 Whole-genome sequencing (WGS)
* 🧬 Whole-exome sequencing (WES)
* 🧪 Cancer genomics
* 🧬 Population genomics
* 🦠 Microbial genomics
* 🧬 Rare-disease research
* 🔬 Comparative genomics
* 🧪 Functional genomics

---

## 🧠 Learning Objectives

This project demonstrates how raw sequencing data can be transformed into interpretable genomic variant information.

Key concepts include:

* FASTQ processing
* Sequence quality control
* Reference alignment
* SAM/BAM processing
* Variant discovery
* VCF files
* Variant filtering
* Variant annotation
* Reproducible genomic workflows

---

## 🔮 Future Improvements

### Pipeline Development

* [ ] Automated FastQC
* [ ] MultiQC reporting
* [ ] Automated read trimming
* [ ] BWA-MEM integration
* [ ] SAMtools processing
* [ ] GATK variant calling
* [ ] BCFtools alternative workflow
* [ ] Variant filtering
* [ ] Variant annotation

### Reproducibility

* [ ] Conda `environment.yml`
* [ ] `requirements.txt`
* [ ] YAML configuration
* [ ] Docker/Singularity support
* [ ] Workflow automation with Snakemake
* [ ] Nextflow implementation
* [ ] Automated testing
* [ ] Version tracking

### Reporting

* [ ] HTML pipeline report
* [ ] QC summary
* [ ] Alignment statistics
* [ ] Variant statistics
* [ ] Interactive plots
* [ ] Variant annotation summary

---

## 🚀 Proposed Production Workflow

```text
                 ┌─────────────────────┐
                 │   FASTQ R1 / R2     │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │      FastQC         │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │   Read Trimming     │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │   BWA-MEM Alignment │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │   SAMtools / BAM    │
                 │   Sort + Index      │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │   Variant Calling   │
                 │     GATK / bcftools │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │ Variant Filtering   │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │ Variant Annotation  │
                 └──────────┬──────────┘
                            │
                            ▼
                 ┌─────────────────────┐
                 │ Final VCF + Report  │
                 └─────────────────────┘
```

---

## 📈 Reproducibility Principles

A research-grade version of this project should record:

* Reference genome version
* Reference genome index
* Sequencing platform
* Tool versions
* Parameters
* Filtering thresholds
* Annotation database versions
* Input sample identifiers

This makes the analysis easier to reproduce and audit.

---

## 👩‍💻 Author

**Salma Hafeez**

Bioinformatics • Computational Biology • Genomics

GitHub: **Bioinformatician-dev**

---

## ⭐ Contributing

Contributions and suggestions are welcome.

Create a feature branch:

```bash
git checkout -b feature/new-variant-module
```

Make your changes, test the workflow, and submit a pull request.

