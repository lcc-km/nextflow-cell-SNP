# nextflow-cell-SNP

A Nextflow-based pipeline for cell SNP (Single Nucleotide Polymorphism) detection and annotation from whole-exome sequencing (WES) data.

## Overview

`nextflow-cell-SNP` is a robust, scalable bioinformatics pipeline designed for germline variant calling from paired-end sequencing data. The pipeline follows GATK (Genome Analysis Toolkit) best practices and includes quality control, alignment, deduplication, base quality score recalibration, variant calling, variant quality score recalibration, and functional annotation.

## Features

- **Quality Control**: FASTP for adapter trimming and quality filtering
- **Alignment**: DRAGMAP (Dragen Read Mapper) for fast and accurate read alignment
- **Deduplication**: Sambamba MarkDuplicates for removing PCR duplicates
- **Quality Metrics**: SamTools stats and Mosdepth for coverage analysis
- **Base Recalibration**: GATK BaseRecalibrator + ApplyBQSR for base quality score recalibration
- **Variant Calling**: GATK HaplotypeCaller for germline SNP and INDEL detection
- **Variant Recalibration**: GATK VQSR (Variant Quality Score Recalibration) for both SNPs and INDELs
- **Variant Annotation**: SnpEff and SnpSift for functional annotation with multiple databases (dbSNP, frequency databases, ClinVar)
- **Container Support**: Docker and Singularity profiles for reproducible execution
- **Modular Design**: Built with nf-core style modules for maintainability

## Pipeline Workflow

```
FASTQ Input
    ↓
[FASTP] - Quality trimming & adapter removal
    ↓
[DRAGMAP] - Read alignment to reference genome
    ↓
[SAMBAMBA MARKDUP] - Mark/remove duplicates
    ↓
[SAMTOOLS STATS + MOSDEPTH] - Quality metrics & coverage analysis
    ↓
[GATK BQSR] - Base Quality Score Recalibration
    ↓
[GATK HAPLOTYPECALLER] - Germline variant calling
    ↓
[GATK VQSR] - Variant Quality Score Recalibration (SNP + INDEL)
    ↓
[SNPEFF + SNPSIFT] - Functional variant annotation
    ↓
Annotated VCF Output
```

## Requirements

- **Nextflow**: >= 23.10.0
- **Container Engine**: Docker or Singularity/Apptainer (recommended)
- **Reference Genome**: hg38 (GRCh38) with associated index files
- **Database Resources**:
  - Known variant sites (dbSNP 146, 1000G phase1, OMNI, HapMap, Mills indels)
  - Annotation databases (dbSNP, frequency database, ClinVar)
  - SnpEff database (GRCh38.99)
  - DRAGMAP hash table
  - Probe BED file (for targeted/WES sequencing)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/lcc-km/nextflow-cell-SNP.git
cd nextflow-cell-SNP
```

2. Install Nextflow (if not already installed):
```bash
curl -s https://get.nextflow.io | bash
```

3. Ensure Docker or Singularity is available on your system.

## Usage

### Basic Usage

```bash
nextflow run cell_snp.nf \
    --input samplesheet.csv \
    --outdir /path/to/output \
    -profile docker
```

### Profiles

- `local`: Run processes locally without containers
- `docker`: Run with Docker containers (recommended)
- `singularity`: Run with Singularity/Apptainer containers

### Example with Singularity

```bash
nextflow run cell_snp.nf \
    --input samplesheet.csv \
    --outdir ./results \
    -profile singularity
```

## Input

### Sample Sheet Format

The pipeline accepts a CSV sample sheet with the following columns:

| Column    | Description                                | Required |
|-----------|--------------------------------------------|----------|
| patient   | Patient identifier                         | Yes      |
| sample    | Sample identifier                          | Yes      |
| lane      | Sequencing lane identifier                 | Yes      |
| fastq_1   | Path to R1 (forward) FASTQ file            | Yes      |
| fastq_2   | Path to R2 (reverse) FASTQ file            | No (for single-end) |

### Example samplesheet.csv

```csv
patient,sample,lane,fastq_1,fastq_2
patient1,sample1,LANE001,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
patient2,sample2,LANE001,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

## Configuration

### Main Parameters

Edit `nextflow.config` to configure the pipeline parameters:

#### Input/Output
- `input`: Path to input sample sheet CSV
- `outdir`: Output directory path

#### Genome Reference
- `genome`: Reference genome build (default: 'hg38')
- `genomes_base`: Base path to reference genome files
- `fasta`: Reference genome FASTA file
- `fasta_fai`: FASTA index file
- `fasta_dict`: Sequence dictionary file

#### Known Variant Sites (for BQSR/VQSR)
- `dbsnp_146`: dbSNP 146 VCF
- `phase1_1000G`: 1000 Genomes phase1 high-confidence SNPs
- `omni_1000G`: 1000 Genomes OMNI SNPs
- `hapmap`: HapMap SNPs
- `indels_1000G`: Mills and 1000G gold standard INDELs

#### Annotation Databases
- `dbSNP`: dbSNP VCF for annotation
- `freq`: Population frequency database VCF
- `clinvar`: ClinVar VCF for clinical significance
- `snpeff_db`: SnpEff database version (default: GRCh38.99)
- `snpeff_cache`: Path to SnpEff cache directory

#### Alignment
- `dragmap_hash`: Path to DRAGMAP hash table directory
- `probe_bed`: Path to probe/target BED file (for WES)

### Resource Configuration

The pipeline includes three process labels with different resource allocations:

| Label           | CPUs | Memory  | Max Forks |
|-----------------|------|---------|-----------|
| `process_high`  | 24   | 60 GB   | 4         |
| `process_medium`| 8    | 16 GB   | 10        |
| `process_low`   | 2    | 4 GB    | 20        |

Global executor settings (local mode):
- Total CPUs: 104
- Total Memory: 370 GB
- Queue size: 20

## Output

The pipeline generates the following output files per sample:

### Alignment & QC
- Trimmed FASTQ files (FASTP output)
- Aligned BAM/CRAM files
- Deduplicated BAM files
- BAM index files (.bai)
- SamTools statistics
- Mosdepth coverage reports

### Variant Calling
- Raw VCF from HaplotypeCaller
- VQSR recalibrated SNP VCF
- VQSR recalibrated INDEL VCF
- Merged final VCF (SNPs + INDELs)

### Annotation
- SnpEff annotated VCF
- SnpSift annotated VCF (with dbSNP, frequency, ClinVar annotations)
- Filtered VCF (PASS variants with DP >= 10 and QUAL >= 20)

## Directory Structure

```
nextflow-cell-SNP/
├── README.md
├── cell_snp.nf              # Main pipeline workflow
├── simple_pipeline.nf       # Simplified example pipeline
├── nextflow.config          # Configuration file
├── samplesheet.csv          # Example sample sheet
├── modules/
│   └── nf-core/             # nf-core style modules
│       ├── fastp/           # FASTP module
│       ├── dragmap/         # DRAGMAP alignment module
│       ├── gatk4/           # GATK4 modules (BQSR, HC, VQSR, etc.)
│       ├── samtools/        # SamTools modules
│       ├── sambamba/        # Sambamba modules
│       ├── mosdepth/        # Mosdepth module
│       ├── snpeff/          # SnpEff annotation module
│       └── snpsift/         # SnpSift annotation module
└── sub_flow/
    └── GATK_GERMLINE_VARIANT_CALLING  # Sub-workflow for GATK germline calling
```

## Pipeline Structure

### Main Workflow (`cell_snp.nf`)

The main pipeline orchestrates the full analysis:

1. **SAMPLE_PARSING**: Parses the input sample sheet and creates sample channels
2. **FASTP_SIMPLE**: Quality control and adapter trimming
3. **DRAGMAP_SIMPLE**: Read alignment to reference genome
4. **SAMBAMBA_MARKDUP**: Mark and remove duplicate reads
5. **Quality Control**: SamTools stats and Mosdepth coverage analysis
6. **GATK BQSR**: Base quality score recalibration
7. **GATK HaplotypeCaller**: Germline variant calling
8. **VQSR**: Variant quality score recalibration (SNP and INDEL separately)
9. **GATK4_MERGEVCFS**: Merge SNP and INDEL VCFs
10. **SNPSIFT_FILTER**: Filter variants based on quality criteria
11. **SNPSIFT_ANNMEM**: Annotate variants with multiple databases

### Sub-workflow (`GATK_GERMLINE_VARIANT_CALLING`)

A reusable sub-workflow that encapsulates the GATK germline variant calling steps (BQSR → HaplotypeCaller → VQSR → Merge).

## Troubleshooting

### Common Issues

1. **Missing reference files**: Ensure all reference genome files and databases exist at the paths specified in the config.

2. **Docker permissions**: If using Docker, ensure your user has proper Docker permissions or use `sudo`.

3. **Memory issues**: For large datasets, you may need to increase memory allocation in the config.

4. **SnpEff database download**: The pipeline can automatically download the SnpEff database if `download_cache = true`.

### Debugging

Run with Nextflow's debug options:
```bash
nextflow run cell_snp.nf --input samplesheet.csv -profile docker -resume -log nextflow.log
```

Use `-resume` to restart from the last successfully completed step.

## License

This project is provided as-is for research purposes.

## Authors

- lucc

## Acknowledgments

- [nf-core](https://nf-co.re/) for module design patterns
- [GATK](https://gatk.broadinstitute.org/) Best Practices
- [Nextflow](https://www.nextflow.io/) for workflow management
