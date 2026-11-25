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
|-------------------|--------|
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

The main output file (OUTFILE above) contains the same columns as the txt output from nexons.py, with additional columns as follows:

| Column  | Contents |
|---------|----------|
| shortIntron_nIsoform | Number of splice patterns that have been collated into this isoform following removal of introns shorter than the minimum width allowed |
| shortIntron_count | Total count coming from splice patterns that contained short introns (this will also have been added to the overall Count column, thus allows the proportion of reads coming from the collated, modified splice patterns to be assessed) |
| shortIntron_nUnique | Number of distinct short introns found within the collated splice patterns |
| shortIntron_minWidth | Minimum width of the short introns found and removed |
| shortIntron_maxWidth | Maximum width of the short introns found and removed |
| MissingExonIntron_TranscriptID | Transcript ID of any annotated isoform(s) that match the splice pattern, except for the omission of a short intron or exon (below the defined thresholds), separated by '_' if > 1 |
| MissingExons | Coordinates of missing exons: exons for each transcript ID above are listed in the same order, separated by '_'; if more than one exon is missing, these are separated by ';'. 'NA' indicates no exons were missing (only intron(s)). | 
| MissingExonWidths | Width(s) of the missing exons, ordered and separated as for 'MissingExons' |
| MissingFirstLastExon | TRUE/FALSE for whether at least one of the missing exons is at the start or end of the transcript; separated by '_' for each of the missing exon/intron transcript IDs |
| MissingIntrons | Coordinates of missing introns, ordered and separated as for 'MissingExons' |
| MissingIntronWidths | Width(s) of missing introns, ordered and separated as for 'MissingExons' |
| Overlapping | Whether the splice pattern is likely to originate from a different gene, which overlaps or is within the defined distance (see above) from this gene, based on matching splice junctions and relative overlaps between the read and genes. 'likely_overlapping' indicates the splice pattern is a better match to the overlapping gene; 'possible_overlapping' indicates that there are shared splice junctions with the overlapping gene but more shared with the original gene, or no clear evidence from splice junctions or relative overlap; 'none' indicates that the original gene is a better match, or no genes within the defined distance (and on the same strand) were found. |
| OverlappingGene | The overlapping/nearby gene(s) from which reads were found to be likely or possibly originating from; otherwise 'none' |
| SpliceMatch | 'all', 'partial' or 'none' indicates whether all/some of the splice donors/acceptors within the splice pattern are annotated. If no splice junctions are present in the splice pattern, 'single_exon'. |
| ReadLocation | Annotation for the transcript feature(s) overlapped by the read. 'complete_tx' indicates the read starts within the first exon, and ends within the last exon, of an annotated transcript (this does not necessarily indicate a matching splice pattern); 'complete_tx_singleExon' indicates that it is fully within an annotated single exon isoform (although the read may contain introns). 'exon', 'intron' or 'outside' indicate that the read is fully within a single exon or intron, or outside of the annotated gene coordinates. 'exon_intron_exon' or 'intron_exon_intron' indicate that the read starts and ends in an exon or intron, respectively, spanning at least one intron or exon in between. 'exon_intron', 'exon_outside' and 'exon_intron_outside' indicate that these features are spanned by the read, but do not necessarily show the order. |
| n_junctions | Number of splice junctions within the isoform |
| TruncationOrigin | If the splice pattern is completely contained within another splice pattern for this gene, ie is a potential truncation product, a comma-separated list of the transcript IDs for those longer splice patterns; otherwise 'none' |
| TruncationStartEndScore | Percentage of reads that terminate before, or within 10 nt of, the splice donor/acceptor that would be used if the read represents a truncation of a longer isoform. Comma separated list for each unique combination of donor/acceptor (where only one end is truncated, the range extends to 0 or Inf at the other end). High score indicates high likelihood that the read represents a truncation, eg 100 indicates that 100% of reads terminate before reaching the next end of the exon. | 
| TruncationScoreMax | Maximum truncation score across all potential truncation origins | 
| TruncationScoreMin | Minimum truncation score across all potential truncation origins |
| TruncationAnnoScore | Comma separated list of truncation scores for any potential truncation origin that is an annotated transcript isoform |
| TruncationAnnoScoreMax | Maximum truncation score across annotated truncation origins |
| TruncationAnnoScoreMin | Minimum truncation score across annotated truncation origins |
| TxEndsWithinTrunc | TRUE/FALSE for whether there is an annotated transcription start/end site before reaching at least one of the splice donor/acceptors used by potential truncation origins. If TRUE, a high truncation score may be less likely to indicate that the transcript is truly truncated (especially if DefinedStart/End is also high) |
| TxEndsWithinAnnoTrunc | TRUE/FALSE for whether there is an annotated transcription start/end site before reaching at least one of the splice donor/acceptors used by annotated potential truncation origins |
| TruncationType | 'threeprime', 'fiveprime' or 'both', to indicate which end is missing. If different truncation origins are truncated at different ends, all possible values will be listed (although not matched) |
| TruncationTypeAnno | As above but only for annotated truncation origins |
| StartEndProportions | For combinations of read start/end coordinates that are above the defined ENDPROPORTION, the proportion of reads that start/end at these positions. List of format 'start-end:proportion' separated by '_' for different start/end combinations |
| DefinedStart | Proportion of reads for which the start is within ENDFLEXQUANT of one of the identified common start sites (ie, is not labelled 'unclear'), and does not lie within 10 nt of a splice acceptor used by a potential 5' truncation origin. A high score indicates that read starts accumulate at specific positions, but are not generated just by loss of the subsequent exon, and are thus more likely to represent true start sites rather than truncation. |
| DefinedEnd | As above but for read ends |
| TranscriptProportions | For annotated transcript isoforms, proportion of reads whether the start and end sites are within ENDFLEXQUANT of the annotated start and end sites. Where more than one isoform shares the same splice pattern but differs in start/end coordinates, all are listed, separated by '_' |


For the summarised transcript start/ends (ENDS_OUTFILE), the output contains the following columns:

| Column  | Contents |
|---------|----------|
| GeneName  | Gene name (if no name in GTF, gene ID) |
| TranscriptID | Transcript ID assigned to this splice pattern. Either derived from the GTF file if the splice pattern is annotated, or 'Variant' followed by a number. If more than one annotated isoform shares an identical splice pattern, all will be listed, separated by ':' |
| Start | Mapping coordinate of read start (or 'unclear' for reads > ENDFLEXQUANT from one of the identified common starts) |
| End | Mapping coordinate of read end (or 'unclear' for reads > ENDFLEXQUANT from one of the identified common starts) |
| SplicePattern | Coordinates of splice pattern. Read start/end coordinates are omitted; splice junctions are separate by ':'; exon start/end coordinates by '-'. Coordinates are provided from 5'-3' with respect to the transcript. |
| File    | Name of input BAM file |
| Count | Number of reads utilising this read start/end combination |
| TranscriptEndMatch | Transcript ID of annotated transcript(s) matching this read start/end and splice pattern |
| Proportion | Proportion of reads in this sample, for this transcript isoform, that use this read start/end combination |



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

