# Variant Sequence Metrics Parser

A lightweight Python utility designed to parse standard VCF files against a reference genome (FASTA) and transcriptome (GTF). This script maps variants to their corresponding transcripts and extracts basic structural and sequence-level metrics, focusing on the translated footprints of frameshift variants.

## Overview
When processing large variant call formats (VCFs), it is often necessary to extract basic thermodynamic and sequence-level properties of the resulting translated peptides. This tool automates the genomic-to-mRNA mapping process and calculates simple biophysical metrics for the divergent sequences.

## Dependencies
This script requires **Python 3.7+** and the following libraries:
* `pandas`
* `biopython`
* `pyfaidx`

You can install the dependencies via pip:
```bash
pip install pandas biopython pyfaidx
