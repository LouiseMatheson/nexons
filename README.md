# Nexons
A pipeline for quantitating transcript level abundances from Nanopore RNA sequencing data

*Note: Nexons has been extensively rewritten and we would recommend using the [new improved, more user-friendly version](https://github.com/s-andrews/nexons). This is being actively developed and will soon allow identification of novel isoform and extraction of read start/end sites, as implemented in this version.*

The nexons.py code is a wrapper around the [chexons](https://github.com/s-andrews/CHBT) program which maps cDNA sequences to genomic segments.  In this pipeline we use chexons to perform a local re-alignment of nanopore sequencing reads which have been genomically mapped with a general purpose mapper to extract specific splice junction coordinates.  We then collate and quantitate these both within and between samples.

The Nexons_postprocessing.R script takes the output from nexons.py, and performs further refinement of the output, to remove or flag potential artefacts, perform additional annotation and collate common read start/end sites. 

## nexons.py

### Requirements

- Python 3
- [chexons](https://github.com/s-andrews/CHBT)
- Samtools

### Basic usage

```
nexons.py -g SIRV6 -o nexons_output/SIRV6.txt --no_chr_prefix SIRV6.gtf SIRV6.fasta test.bam
```

The nexons.py script takes BAM files as input (for example, mapped with Minimap2), together with a FASTA and GTF file and a gene name or ID. Reads mapped to the gene of interest are extracted, realigned, and splice variants identified and quantified, with collation of variants that differ by less than a specified threshold. 

Splice junction correction is also performed based on annotated junctions, and, optionally, user-provided splice junctions derived from short-read data, in [this format]. 

<img width="1056" height="866" alt="nexons_docs" src="https://github.com/user-attachments/assets/c7a9ffa7-e021-4d03-a332-eef62010c787" />
