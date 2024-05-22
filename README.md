![](https://img.shields.io/badge/language-R-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/PhD-thesis) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/PhD-thesis)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/PhD-thesis) ![GitHub](https://img.shields.io/github/license/viktormiok/PhD-thesis)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/PhD-thesis) 

Targeted expression profiling identifies mechanisms of potential breast carcinogens
========================

- [Scope](#scope)
- [Summary](#summary)
- [License](#license)
- [Assay](#assay)
- [Protocol](#protocol)
  
# Scope

Differential gene expression and pathway analysis of MCF-7 cell lines exposed to multiple chemicals belonging to specific chemical classes at different doses under estrogen-starved and treated conditions employing breast cancer panel targeting 500 genes using TempO-Seq platform and a targeted gene expression panel, Breast Carcinogen Screen (BCScreen).

The study aims to identify mechanisms of potential breast carcinogens using targeted expression profiling, differential expression, and pathway analysis.

# Summary

Over 200 chemicals have been associated with altered mammary gland development and cancer in rodent studies. Yet, the molecular mechanisms by which this diverse set of chemicals affects the breast remain poorly understood. We examined gene expression profiles of eight known and suspected mammary toxicants using the TempO-Seq platform and a targeted gene expression panel, Breast Carcinogen Screen (BCScreen). The BCScreen panel has 500 genes with at least 30 genes strongly associated with each of the 14 biologic processes relevant to breast cancer and mammary gland development. 

We compared gene expression profiles in MCF-7 cells treated for 24hrs with four concentrations (1 nM to 10 uM) of the following chemicals: three estrogenic chemicals (genistein, bisphenol A, butyl benzyl phthalate); the selective estrogen receptor modulator, tamoxifen (TAM); 1,4 benzoquinone (BQ), a genotoxic metabolite of the rodent mammary carcinogen benzene; perfluorooctanoic acid (PFOA), a persistent pollutant that alters mammary gland development; and two closely related perfluorinated chemicals that have not been evaluated for mammary toxicity/effects (PFHXA, PFNA). We identified expression patterns activated by mammary toxicants with known ER agonist activity (genistein, BPA, BBP). PFOA and BQ and TAM, mammary toxicants with antagonists or no known direct ER activity suppressed some of these genes. 

The known weak estrogens genistein, BPA, and BBP had similar gene expression profiles and were enriched for genes associated with cell cycle and DNA repair, particularly at higher doses. Notably, PFOA, Tam, and BQ had greater gene expression changes at the lowest (1nM) vs. higher doses and were enriched for genes associated with the cell cycle. For context, serum PFOA levels in the general US population are generally five to ten times higher than the 1 nM dose where we observed these changes (i.e., > 0.4 micrograms/L). PFOA profiles were enriched for several transcription factors including E2F1 and E2F4, which regulate cell cycle genes and important processes in mammary development and tumor progression. Utilizing the BCScreen gene panel, we have identified mechanistic targets relevant to breast cancer perturbed by mammary toxicants.

# Assay
Sequencing libraries for targeted panels were generated as described briefly and depicted in the figure below. In TempO-Seq, each Detector Oligo (DO) consists of a sequence complementary to an mRNA target plus a universal (i.e. same for every targeted gene) primer binding site. They anneal in immediate juxtaposition to each other on the targeted RNA template such that they can be ligated together. Ligated detector oligos are PCR-amplified using a primer set (single-plex PCR reaction, with a single primer pair for each sample) that introduces both the adaptors required for sequencing and a sample-specific barcode. The barcode sequences flank the target sequence and are inserted appropriately into the standard Illumina adaptors to permit standard dual-index sequencing of the barcodes and deconvolution of sample-specific reads from the sequencing data using the standard Illumina software. All the PCR-amplified and barcoded samples are pooled into a single library for sequencing. Sequencing reads are demultiplexed using the standard sequencing instrument software for each sample using the barcodes to give a FASTQ file for each.

# Protocol
TempO-Seq data are analyzed using the Tempo-SeqR software package using the statistical computing language R. The input for TempO-Seq data analysis is a folder of zipped FASTQ files. Each FASTQ file contains the reads and quality scores for one sample. Each FASTQ file is aligned using the Bowtie algorithm to a pseudo-transcriptome input by the user. The output for TempO-Seq data analysis is a table of counts with each column representing a sample and each row representing a gene.

# License

__`TempOseq_CarcinoChemMCF-7`__ is distributed under the GPL-3.0 License. Please read the license before using __`TempOseq_CarcinoChemMCF-7`__, which is distributed in the `LICENSE` file.
