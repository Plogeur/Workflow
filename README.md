
# Variant Calling Pipeline for Illumina Sequencing Data
This pipeline is designed to perform variant calling from Illumina sequencing data using a Snakemake workflow. It includes steps for reads quality check, trimming, reference indexing, alignment, variant calling, and filtering.

## Usage
```bash
snakemake --cores <nb_core_max>
```

## Dependencies
- FastQC
- Fastp
- BWA
- Samtools
- Picard
- Fgbio
- GATK
- Qualimap

## Pipeline Steps
Reads Quality Check: Perform quality check on input reads using FastQC.

<p align="center">
  <img src="https://github.com/Plogeur/Worflow/blob/main/img/worflow.jpg" alt="Nom_de_l'image">
</p>

1.Reads Trimming: Trim adapter sequences and low-quality bases from reads using Fastp.

2.Indexing Reference: Index the reference genome if necessary.

3.Alignment: Align trimmed reads to the reference genome using BWA.

4.Statistics: Generate alignment statistics using Samtools and Picard tools.

5.MarkDuplicates: Mark duplicate reads using Picard tools.

6.Clip Overlaps: Clip overlapping reads using Fgbio.

7.Variant Calling: Call variants using HaplotypeCaller from GATK.

8.Genotype gVCF: Genotype the produced GVCF file using GATK.

9.Index Feature File: Index the produced VCF file.

10.Filter SNPs: Apply filtering criteria for SNPs using GATK VariantFiltration.

11.Filter InDels: Apply filtering criteria for InDels using GATK VariantFiltration.

12.Merge Variants: Merge filtered SNPs and InDels into a single VCF file using GATK MergeVcfs.

13.Compress Files: Compress VCF and gVCF files using bgzip.

## Note
Ensure all dependencies are installed and accessible in the system's PATH.
Adjust parameters and filtering criteria as per your specific requirements.
