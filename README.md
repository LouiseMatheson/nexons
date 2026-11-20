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

<img width="528" height="433" alt="nexons_docs" src="https://github.com/user-attachments/assets/c7a9ffa7-e021-4d03-a332-eef62010c787" />

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
| TranscriptID | Transcript ID assigned to this splice pattern. Either derived from the GTF file if the splice pattern is annotated, or 'Variant' followed by a number. If more than one annotated isoform shares an identical splice pattern, all will be listed, separated by ':' |
| Count  | Number of reads in the BAM file assigned to this splice pattern, following collation of similar splice patterns |
| Starts | Comma separated list of coordinates of the first mapped position in the transcript, ie the transcription start site |
| Ends | Comma separated list of coordinates of the last mapped position in the transcript, ie the transcription end site. Read start and end lists are provided in the same order such that coordinates at the same position in the list originate from the same read. |     
| MergedIsoforms | If any annotated isoforms do not share the exact splice pattern, but were collated with this isoform due to their splice junctions differing by less than the allowed flexibility, their transcript ID(s) will be listed here |


------------------


## Nexons postprocessing.R

### Basic usage

```
Rscript --vanilla Nexons_postprocessing.R [options] SIRV6.gtf nexons_output/
```
This script takes a directory containing the tab-delimited text file outputs from nexons.py (nexons_output in the command above) as input, in addition to a GTF file (which should be the same file as provided to nexons.py), and performs several steps to clean up, flag potential artefacts and further annotate the data. It also collates the read start and end sites, identifying and quantifying common (high density) sites. 

<img width="592" height="575" alt="nexonspost_docs" src="https://github.com/user-attachments/assets/7e55dc2a-10a7-4996-aab2-62239bc0b55b" />

It can be run across multiple genes (ie nexons output files) together, but since each gene is analysed separately it supports parallelisation. Performance is best, however, when run across multiple batches of nexons output files. An example wrapper for generating and running across batches using the Slurm workload manager is provided (nexons_postprocessing.sh - see details below), but may need to be modified for different environments.

### Optional arguments

| Argument (default) | Description |
---------------------------
| -o OUTFILE (nexons_output_refined.txt) | The file to write the output count table to  |
| -q ENDS_OUTFILE (OUTFILE + '_completeEnds' tag) | The file to write the collated/quantitated read start and ends to  |
| -f FLEXIBILITY (5) | Maximum difference between exon boundaries for splice patterns to be merged (applies to individual coordinates: as long as differences for all sites are within the flexibility, the splice patterns will be merged). Should be the same as was used for nexons.py |
| -m MINEXONS (2) | Minimum number of exons which must be present to count a transcript. Should be the same as for nexons.py -e |
| -i MININTRONWIDTH (50) | Introns below this width will be removed from splice patterns; splice patterns that would match an annotated isoform if an intron of this size were included will be flagged |
| -x MISSEDEXONWIDTH (30) | Splice patterns that match an annotated isoform other than excluding exons below this width will be flagged |
| -g OVERLAPGAP (500) | Genes within this distance will be considered potentially overlapping, and splice patterns checked to see whether they are likely to originate from any overlapping gene |
| -r ENDREADS (5) | Start/end sites must have at least this number of reads with a 5 base window to be designated a common/high density site (across all samples and isoforms) |
| -e ENDFLEX (10) | Flexibility used when defining common start/end sites. High density sites closer than this distance from a site with higher density will be collated with that site. |
| -d ENDFLEXQUANT (30) | Flexibility used when quantitating common start/end sites. Read start/ends are assigned to their closest high density site, if there is any site within this distance (otherwise labelled 'unclear'). This is also the flexibility used when determining the proportion of reads start/ending at sites corresponding to annotated isoforms, even if another, closer, high density unannotated site was found. |
| -p ENDPROPORTION (0.1) | Start/end combinations above this proportion will be included in OUTFILE (all high density sites are included in ENDS_OUTFILE regardless of proportion) |
| -c CORES (1) | Number of cores to parallelise over |
| --no_header | If included, OUTFILE and ENDS_OUTFILE will be written out without column headers. Used by wrapper for batch processing to facilitate concatenation. |
| -z FINISHED_LOG | Filename for the log file to show that processing has completed (used by wrapper for batch processing) |
| --log_finish | Always TRUE if -z is provided. If --log_finish is included but -z is not, default filename 'FINISHED.log' will be used. |


### Output

Identical to nexons.py txt output but with additional columns

completeEnds output...


## nexons_postprocessing.sh

### Basic usage

```
nohup nexons_postprocessing.sh Nexons_postprocessing.R SIRV6.gtf nexons_output/ &
```

This provides a wrapper for generating and running across batches when using a Slurm workload manager, although it may need to be modified for different working environments. 

It first generates subdirectories within the nexons_output directory, each containing a subset of the nexons.py output txt files. A separate job to run Nexons_postprocessing.R is then submitted for each of these subdirectories, omitting column headers from the outputs and naming each output by its subset number, and a folder (logfiles_nexons_postprocessing) generated in which log files for each job are recorded. The column headers are written separately to a '...subset0...' file. A 'FINSHED_subset[x].log' file is also generated as each job completes; once this has been detected for all subsets, the output files are concatenated together with the column headers, the 'FINISHED...' files removed and the file structure within the nexons output directory returned to its original state. Note that if all of the jobs do not complete successfully, this clean-up will need to be done manually. 

Arguments are identical to Nexons_postprocessing.R except:
- -s can be provided to define batch size (default 100)
- -b can be provided to define memory allowance for each job (default 50G)
- --no_header is always used to facilitate concatenation of output files
- -z (and therefore --log_finish) is defined by the wrapper such that log files showing completion of each batch are generated and recognised

