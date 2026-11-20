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

Since the nexons.py script is run on one gene at a time, it would typically be run as a loop over a gene list, with gene names used as the output file names, and saving into an output directory. It should be run on all BAM files for an experiment together, such that collation and naming of novel isoforms is performed across all samples together. 

<img width="1056" height="866" alt="nexons_docs" src="https://github.com/user-attachments/assets/c7a9ffa7-e021-4d03-a332-eef62010c787" />

Splice junction correction is also performed based on annotated junctions, and, optionally, user-provided splice junctions derived from short-read data, in [this format](https://github.com/LouiseMatheson/nexons/blob/master/test/example_splice_junctions.txt). 

### Optional arguments

```
usage: nexons.py [-h] [--outfile OUTFILE] [--both_out] [--gtf_out] [--flexibility FLEXIBILITY] [--mincount MINCOUNT]
                 [--mincoverage MINCOVERAGE] [--mapthreshold MAPTHRESHOLD] [--minexons MINEXONS] [--gene GENE]
                 [--splice_sites SPLICE_SITES] [--direction DIRECTION] [--no_chr_prefix] [--verbose] [--verbose_proportions]
                 [--report_all] [--quiet] [--suppress_warnings] [--version]
                 gtf fasta bam [bam ...]

positional arguments:
  gtf                   A GTF file containing the genes you want to analyse
  fasta                 A multi-fasta file containing the genome sequence
  bam                   One or more BAM files to quantitate

options:
  -h, --help            show this help message and exit
  --outfile OUTFILE, -o OUTFILE
                        The file to write the output count table to
  --both_out, -b        Write out to gtf and custom nexons format
  --gtf_out, -t         Write to gtf file instead of custom nexons format
  --flexibility FLEXIBILITY, -f FLEXIBILITY
                        Maximum difference between exon boundaries for splice patterns to be merged (applies to individual coordinates: as
                        long as differences for all sites are within the flexibility, the splice patterns will be merged)
  --mincount MINCOUNT, -m MINCOUNT
                        Minimum number of observations for any given variant to report
  --mincoverage MINCOVERAGE, -c MINCOVERAGE
                        Minimum proportion of the gene which must be covered by a variant
  --mapthreshold MAPTHRESHOLD, -mt MAPTHRESHOLD
                        Minimum proportion of the sequence that must be mapped to a transcript
  --minexons MINEXONS, -e MINEXONS
                        Minimum number of exons which must be present to count a transcript
  --gene GENE, -g GENE  Name or ID of a single gene to quantitate
  --splice_sites SPLICE_SITES
                        Text file containing known splice donors and acceptors
  --direction DIRECTION, -d DIRECTION
                        Directionality of the library (none, same, opposing)
  --no_chr_prefix, -nc  Don't add a 'chr' prefix in temp bed file (required when the fasta file used for mapping does not use this prefix)
  --verbose, -v         Produce more verbose output
  --verbose_proportions, -vp
                        Print out all proportions of sequences mapped
  --report_all          Report all transcripts, including seeded ones with no counts. Only works with --gtf_out or --both_out. Writes to a
                        separate file named match_info.txt
  --quiet               Suppress all messages
  --suppress_warnings   Suppress warnings (eg about lack of names or ids)
  --version             Print version and exit

```

### Output

The default output from nexons.py is a tab-delimited text file. Optionally, a GTF file containing the identified transcript isoforms can be output in addition or instead. 

The tab-delimited text output contains one row per sample (BAM file) per transcript isoform, with columns as follows:

| Column  | Contents |
|---------|----------|
| File    | Name of input BAM file |
| GeneID  | Gene ID. For genes where multiple gene IDs in the GTF file are assigned the same gene name, all gene IDs will be listed, separated by ':' |
| GeneName  | Gene name (if no name in GTF, gene ID) |
| Chr    | Chromosome |
| Strand | Strand |
| SplicePattern | Coordinates of splice pattern. Read start/end coordinates are omitted; splice junctions are separate by ':'; exon start/end coordinates by '-'. Coordinates are provided from 5'-3' with respect to the transcript. |
| TranscriptID | Transcript ID assigned to this splice pattern. Either derived from the GTF file if the splice pattern is annotated, or 'Variant' followed by a number |
| Count  | Number of reads in the BAM file assigned to this splice pattern, following collation of similar splice patterns |
| Starts | Comma separated list of coordinates of the first mapped position in the transcript, ie the transcription start site |
| Ends | Comma separated list of coordinates of the last mapped position in the transcript, ie the transcription end site. Read start and end lists are provided in the same order such that coordinates at the same position in the list originate from the same read. |     
| MergedIsoforms
