#!/usr/bin/env Rscript
renv::load("/bi/group/turner/Louise/Ozge_data/nexons_postprocessing")

library(progressr)
library(rtracklayer)
library(tidyverse)
library(furrr)

options(progressr.enable = T)

handlers(list(handler_progress(enable = T)))

#############################################################
# Parameters
#############################################################

print(Sys.time())
print("Reading parameters...")

Args = commandArgs(trailingOnly=TRUE)

if(length(Args) > 1) {
  gtf_file <- Args[length(Args)-1]
} else { gtf_file <- "../../Annotations/gencode.vM28.annotation.gtf.gz" }
print(paste("GTF file:", gtf_file))

if(length(Args) > 1) {
  nexons_directory <- Args[length(Args)]
} else { nexons_directory <- "per_gene_outs/" }
print(paste("Nexons output directory:", nexons_directory))


if("-o" %in% Args) {
  output_file <- Args[which(Args == "-o") + 1]
} else { output_file <- "nexons_output_refined.txt" }
print(paste("Output file:", output_file))

if("-q" %in% Args) {
  ends_outfile <- Args[which(Args == "-q") + 1]
} else { ends_outfile <- "nexons_output_refined_completeEnds.txt" }
print(paste("Output file (complete ends):", ends_outfile))

if("-z" %in% Args) {
  FINISHED <- Args[which(Args == "-z") + 1]
} else { FINISHED <- "FINISHED_subset.log" }
print(paste("Completion logfile:", FINISHED))

if("-f" %in% Args) {
  flexibility <- as.numeric(Args[which(Args == "-f") + 1])
} else { flexibility <- 5 }
print(paste("Flexibility:", flexibility))

if("-m" %in% Args) {
  min_exons <- as.numeric(Args[which(Args == "-m") + 1])
} else { min_exons <- 2 }
print(paste("Minimum exons:", min_exons))

if("-i" %in% Args) {
  minIntronWidth <- as.numeric(Args[which(Args == "-i") + 1])
} else { minIntronWidth <- 50 } # introns below this will be removed from splice patterns (but splice patterns that would match an annotated isoform if this intron were included will be flagged)
print(paste("Minimum intron width:", minIntronWidth))


if("-x" %in% Args) {
  missedExonWidth <- as.numeric(Args[which(Args == "-x") + 1])
} else { missedExonWidth <- 30 } # will check for splice patterns that match an annotated isoform other than excluding exons below this width
print(paste("Missing exon width:", missedExonWidth))

if("-g" %in% Args) {
  overlapGap <- as.numeric(Args[which(Args == "-g") + 1])
} else { overlapGap <- 500 } # genes closer than this distance will be designated 'overlapping' and splice patterns checked in case they might originate from the overlapping gene
print(paste("Overlap gap:", overlapGap))


if("-e" %in% Args) {
  flexibility_ends <- as.numeric(Args[which(Args == "-e") + 1])
} else { flexibility_ends <- 10 } # ends < this flexibility apart will be merged
print(paste("Flexibility for identifying common ends:", flexibility_ends))

if("-d" %in% Args) {
  flexibility_ends2 <- as.numeric(Args[which(Args == "-d") + 1])
} else { flexibility_ends2 <- 30 } # ends assigned to closest collapsed end, if < this flexibility away. Same flexibility used for calculating proportion of reads matching annotated isoforms
print(paste("Flexibility for collating ends:", flexibility_ends2))

if("-p" %in% Args) {
  start_end_proportion <- as.numeric(Args[which(Args == "-p") + 1])
} else { start_end_proportion <- 0.1 } # start/end combinations above this proportion will be included in the collated (1 row per splice pattern) output
print(paste("Proportion above which start/end combinations are reported:", start_end_proportion))

if("-r" %in% Args) {
  start_end_minreads <- as.numeric(Args[which(Args == "-r") + 1])
} else { start_end_minreads <- 5 } # Starts/ends must have at least this number of reads with a 5 base window to be included (across all samples and isoforms)
print(paste("Number of reads required for designation of start/end coordinates:", start_end_minreads))

if("-c" %in% Args) {
  n_cores <- as.numeric(Args[which(Args == "-c") + 1])
} else { n_cores <- 1 } # Will process in parallel over this number of cores
print(paste("Number of cores:", n_cores))

try(future::plan(multisession, workers = (2*n_cores)))

options(future.globals.maxSize = 5000000000)



identifyPotentialTruncations <- function(parsed_splices, gtf, flexibility = 10){
  # identify any annotated isoforms that weren't detected, as there could be truncations of these
  undetected_isoforms <- identifyUndetectedTranscripts(gtf, parsed_splices) # note this assumes there is a single gene contained in parsed_splices - will only detect for the gene on the first row
  # pull out starts and ends of transcripts
  gtf %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::distinct(transcript_id, start, end) -> transcripts
  
  # extract splice patterns and variant numbers from parsed_splices.
  #
  # put them into a list where names are variant names/numbers
  parsed_splices %>%
    dplyr::mutate(n_junctions = stringr::str_count(SplicePattern, ":")) -> parsed_splices
  assertthat::validate_that(all(parsed_splices$n_junctions > 0), msg = "Cannot have splice patterns with fewer than 1 junctions, removing these...")
  #
  all_junctions <- parsed_splices %>%
    dplyr::filter(n_junctions > 0) %>%
    dplyr::distinct(TranscriptID, SplicePattern) %>%
    tibble::deframe() %>%
    purrr::map(.f=spliceStartsEnds)
  
  all_junctions_with_anno <- parsed_splices %>%
    dplyr::filter(n_junctions > 0) %>%
    dplyr::distinct(TranscriptID, SplicePattern) %>%
    dplyr::bind_rows(undetected_isoforms) %>%
    tibble::deframe() %>%
    purrr::map(.f=spliceStartsEnds)
  
  undetected_isoforms %>%
    dplyr::bind_rows(dplyr::select(parsed_splices, TranscriptID, SplicePattern)) %>%
    dplyr::distinct() %>%
    dplyr::filter(!grepl("^Variant", TranscriptID)) %>%
    dplyr::left_join(transcripts, by = c(TranscriptID = "transcript_id")) -> all_isoforms
  
  trunc_results <- purrr::map(.x=all_junctions, function(trunc) {
    purrr::map(.x=all_junctions_with_anno, .f=checkTruncation, potential_truncs=trunc, flexibility=flexibility) %>%
      dplyr::bind_rows()
  })
  
  individual_results <- purrr::map(trunc_results, ~ {
    .x %>% 
      tibble::add_column(TruncationOrigin = names(all_junctions_with_anno)) %>% 
      dplyr::filter(truncation == T) %>%
      dplyr::group_by(TruncationOrigin) %>%
      dplyr::mutate(TruncationStart = if_else(is.na(TruncationStart), {
        if(parsed_splices$Strand[1] == "+") { 0 } else { Inf }
      }, TruncationStart)) %>%
      dplyr::mutate(TruncationEnd = if_else(is.na(TruncationEnd), {
        if(parsed_splices$Strand[1] == "+") { Inf } else { 0 }
      }, TruncationEnd)) %>%
      dplyr::ungroup()
  })
  
  individual_results %>%
    dplyr::bind_rows(.id = "TranscriptID") -> tibble_results
  
  if(nrow(tibble_results) > 0) {
    tibble_results %>%
      dplyr::mutate(OriginComb = if_else(grepl("^Variant", TruncationOrigin), "Novel", TruncationOrigin)) %>%
      dplyr::group_by(TranscriptID, TruncationStart, TruncationEnd, OriginComb) %>%
      dplyr::summarise(TruncationOrigin = paste(TruncationOrigin, collapse = ","), .groups = "drop") %>%
      tidyr::nest(TruncationStartEnd = -TranscriptID) -> tibble_results
    parsed_splices %>%
      dplyr::select(TranscriptID, SplicePattern, Strand, Starts, Ends) %>%
      dplyr::group_by(TranscriptID, SplicePattern, Strand) %>%
      dplyr::summarise(
        Starts = paste(Starts, collapse = ","),
        Ends = paste(Ends, collapse = ","),
        .groups = "drop"
      ) %>%
      dplyr::left_join(transcripts, by = c(TranscriptID = "transcript_id")) %>%
      dplyr::mutate(start = dplyr::if_else(is.na(start), purrr::map2_dbl(SplicePattern, Strand, function(spl, strnd) {
        if(strnd == "+") {
          max(all_isoforms$start[grepl(paste0("^", sub(":.*", "", spl)), all_isoforms$SplicePattern)])
        } else {
          min(all_isoforms$start[grepl(paste0(sub(".*:", "", spl), "$"), all_isoforms$SplicePattern)])
        }
      }), start)) %>% 
      dplyr::mutate(end = dplyr::if_else(is.na(end), purrr::map2_dbl(SplicePattern, Strand, function(spl, strnd) {
        dplyr::if_else(strnd == "+", 
                min(all_isoforms$end[grepl(paste0(sub(".*:", "", spl), "$"), all_isoforms$SplicePattern)]),
                max(all_isoforms$end[grepl(paste0("^", sub(":.*", "", spl)), all_isoforms$SplicePattern)])
        )
      }), end)) %>%
      dplyr::mutate(
        start = dplyr::if_else(is.infinite(start), NA_real_, start),
        end = dplyr::if_else(is.infinite(end), NA_real_, end)
      ) %>%
      dplyr::right_join(tibble_results, by = "TranscriptID") %>% 
      dplyr::mutate(TruncationStartEnd = purrr::pmap(list(TruncationStartEnd, Starts, Ends, start, end, TranscriptID, Strand), function(data, starts, ends, annoStart, annoEnd, txid, strnd) {
        if(nrow(data) == 0) {
          tibble(
            TruncationStartEndScore = NA_character_,
            TruncationScoreMax = 0,
            TruncationScoreMin = 0,
            TruncationAnnoScore = NA_character_,
            TruncationAnnoScoreMax = 0,
            TruncationAnnoScoreMin = 0)
        } else {
          stringr::str_split(starts, ",") %>% unlist() %>% as.integer() -> read_starts
          stringr::str_split(ends, ",") %>% unlist() %>% as.integer() -> read_ends
          data %>%
            dplyr::group_by(TruncationStart, TruncationEnd, OriginComb) %>%
            dplyr::mutate(TruncType = case_when(
              TruncationStart %in% c(0, Inf) ~ "threeprime",
              TruncationEnd %in% c(0,Inf) ~ "fiveprime",
              TRUE ~ "both"
            )) %>%
            dplyr::mutate(score = signif(100*sum(purrr::map2_lgl(read_starts, read_ends, ~ dplyr::between(.x, min(TruncationStart, TruncationEnd)-10, max(TruncationStart, TruncationEnd)+10) & dplyr::between(.y, min(TruncationStart, TruncationEnd)-10, max(TruncationStart, TruncationEnd)+10)))/length(read_starts), 4)) %>%
            dplyr::mutate(TruncationStartEndScore = paste0(TruncationStart, "-", TruncationEnd, ":", score)) %>%
            dplyr::mutate(TruncationAnnoScore = paste0(OriginComb, ":", score)) %>%
            dplyr::mutate(TruncationOutsideTxEnds = dplyr::between(annoStart, min(TruncationStart, TruncationEnd)-10, max(TruncationStart, TruncationEnd)+10) & dplyr::between(annoEnd, min(TruncationStart, TruncationEnd)-10, max(TruncationStart, TruncationEnd)+10)) %>%
            dplyr::ungroup() %>% 
            dplyr::summarise(
              TruncationOrigin = paste(TruncationOrigin, collapse = ","),
              TruncationStartEndScore = paste(unique(TruncationStartEndScore), collapse = ","), 
              TruncationScoreMax = max(score), 
              TruncationScoreMin = min(score),
              TruncationAnnoScore = dplyr::if_else(sum(OriginComb == "Novel") < n(), paste(TruncationAnnoScore[OriginComb != "Novel"], collapse = ","), NA_character_), 
              TruncationAnnoScoreMax = dplyr::if_else(sum(OriginComb == "Novel") < n(), max(score[OriginComb != "Novel"]), 0), 
              TruncationAnnoScoreMin = dplyr::if_else(sum(OriginComb == "Novel") < n(), min(score[OriginComb != "Novel"]), 0),
              TxEndsWithinTrunc = sum(TruncationOutsideTxEnds, na.rm = T) > 0,
              TxEndsWithinAnnoTrunc = dplyr::if_else(sum(OriginComb == "Novel") < n(), sum(TruncationOutsideTxEnds[OriginComb != "Novel"], na.rm = T) > 0, F),
              TruncationType = paste(unique(TruncType), collapse = "_"),
              TruncationTypeAnno = dplyr::if_else(sum(OriginComb == "Novel") < n(), paste(unique(TruncType[OriginComb != "Novel"]), collapse = "_"), NA_character_),
              TruncStartsAnno = dplyr::if_else(sum(OriginComb == "Novel") < n(), paste(unique(TruncationStart[OriginComb != "Novel"]), collapse = "_"), NA_character_),
              TruncEndsAnno = dplyr::if_else(sum(OriginComb == "Novel") < n(), paste(unique(TruncationEnd[OriginComb != "Novel"]), collapse = "_"), NA_character_)
            )
        }
      })) %>% 
      tidyr::unnest(TruncationStartEnd) %>%
      dplyr::select(TranscriptID, contains("Trunc")) -> tibble_results
    
  } else {
    tibble_results <- tibble::tibble(TranscriptID = character(), TruncationOrigin = character(), TruncationStartEndScore = character(), TruncationScoreMax = numeric(), TruncationScoreMin = numeric(), TruncationAnnoScore = character(), TruncationAnnoScoreMax = numeric(), TruncationAnnoScoreMin = numeric(), TxEndsWithinTrunc = logical(), TxEndsWithinAnnoTrunc = logical(), TruncationType = character(), TruncationTypeAnno = character(), TruncStartsAnno = character(), TruncEndsAnno = character())
  }
  
  
  
  parsed_splices %>%
    dplyr::left_join(tibble_results, by = "TranscriptID") %>%
    dplyr::mutate(TruncationOrigin = if_else(is.na(TruncationOrigin), "none", TruncationOrigin))
  
}


identifyUndetectedTranscripts <- function(gtf, parsed_splices) {
  gtf %>%
    dplyr::filter(gene_name == parsed_splices$GeneName[1]) %>%
    dplyr::filter(!transcript_id %in% parsed_splices$TranscriptID) %>%
    dplyr::filter(type == "exon") %>%
    dplyr::arrange(exon_number) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::mutate(junction = dplyr::case_when(
      strand == "+" & exon_number == 1 ~ as.character(end),
      strand == "-" & exon_number == 1 ~ as.character(start),
      strand == "+" & exon_number == max(exon_number) ~ as.character(start),
      strand == "-" & exon_number == max(exon_number) ~ as.character(end),
      strand == "+" ~ paste(start, end, sep = "-"),
      strand == "-" ~ paste(end, start, sep = "-")
    )) %>%
    dplyr::summarise(SplicePattern = paste(junction, collapse = ":"), .groups = "drop") %>%
    dplyr::rename(TranscriptID = transcript_id) %>%
    dplyr::distinct(TranscriptID, SplicePattern)
}

#' spliceStartsEnds
#' Takes a splice string in the format "8381:8455-8585:10859-10991:11312" and converts it to a tibble of start and end points of splice junctions.
#' Used internally by \code{\link{identifyPotentialTruncations}} function.
#'
#' @param splice_string set of splices in the format 8381:8455-8585:10859-10991:11312"
#'
#' @return tibble of start and end points of splice junctions
#' @export
#'
#' @examples
#' splice <- "8381:8455-8585:10859-10991:11312"
#' spliceStartsEnds(splice)
spliceStartsEnds <- function(splice_string){
  splices <- stringr::str_split(splice_string, pattern = "-")[[1]]
  tibble::enframe(splices, name = "junction_no") %>%
    tidyr::separate(value, into = c("start", "end"), sep = ":", convert=TRUE)
}

#' Check junction matches
#' Check if 2 sets of splice junctions match.
#'
#' @param set1 tibble with column names junction_no, start, end
#' @param set2 tibble with column names junction_no, start, end
#' @param flexibility integer value - number of bases away to still be classified as a match. For exact matches set flexibility to 0.
#'
#' @return vector of TRUE and FALSE values
#' @export
#'
#' @examples
#' NA
checkJunctionMatches <- function(set1, set2, flexibility){
  set1$start <= set2$start+flexibility &
    set1$start >= set2$start-flexibility &
    set1$end <= set2$end+flexibility &
    set1$end >= set2$end-flexibility
}


#' checkTruncation
#'
#' Checks one set of splice junctions against another set, to see if one is a potential truncation of the other.
#' Does an initial check for the first splice junction location, then do an exact test (within flexibility parameter) for vectors of the same length.
#'
#' @param other_set tibble - one set of splice junctions
#' @param potential_truncs tibble - one set of splice junctions to test to see if they might be truncations
#' @param flexibility integer value - number of bases away to still be classified as a match. For exact matches set flexibility to 0.
#'
#' @return TRUE or FALSE value
#' @export
#'
#' @examples
#' NA
checkTruncation <- function(other_set, potential_truncs, flexibility) {
  
  if(nrow(potential_truncs) >= nrow(other_set)) return(tibble(truncation = FALSE, TruncationStart = NA_integer_, TruncationEnd = NA_integer_)) # if potential_trunc is longer then reject
  
  # check for first junction match
  junction1 <- checkJunctionMatches(dplyr::slice(potential_truncs, 1), other_set, flexibility=flexibility)
  if(all(junction1 == FALSE)) return(tibble(truncation = FALSE, TruncationStart = NA_integer_, TruncationEnd = NA_integer_)) # reject if no match
  
  # extract the first junction number that mapped
  junction1_no <- dplyr::filter(other_set, junction1) %>% dplyr::pull(junction_no)
  
  if(nrow(potential_truncs) > (1+nrow(other_set)-junction1_no)) return(tibble(truncation = FALSE, TruncationStart = NA_integer_, TruncationEnd = NA_integer_)) # if too few remaining junctions after match
  
  # subset consecutive junctions
  subset <- other_set %>%
    dplyr::slice(junction1_no:dplyr::n()) %>%
    dplyr::slice_head(n=nrow(potential_truncs))
  
  # test whole set to see if each junction matches
  x <- checkJunctionMatches(potential_truncs, subset, flexibility=flexibility)
  
  if(any(x == FALSE)) return(tibble(truncation = FALSE, TruncationStart = NA_integer_, TruncationEnd = NA_integer_))
  if(junction1_no == 1) { # capture the expected start and end of the exons at either end, if it is a truncation -> starts/ends should be within this
    TruncationStart <- NA_integer_
  } else {
    TruncationStart = other_set %>% dplyr::slice(junction1_no - 1) %>% pull(end)
  }
  if((junction1_no + nrow(potential_truncs)) > nrow(other_set)) {
    TruncationEnd <- NA_integer_
  } else {
    TruncationEnd = other_set %>% dplyr::slice(junction1_no + nrow(potential_truncs)) %>% pull(start)
  }
  return(tibble(truncation = TRUE, TruncationStart = TruncationStart, TruncationEnd = TruncationEnd)) # potential truncation
}


#############################################################
# Read in GTF file
#############################################################

print(Sys.time())
print("Reading in GTF file...")

import(gtf_file) %>%
  tidyr::as_tibble() -> gtf

if("exon_number" %in% colnames(gtf)) {
  gtf %>%
    dplyr::mutate(exon_number = as.integer(exon_number))-> gtf
} else {
  gtf %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::mutate(exon_number = if_else(strand == "+", row_number(), rev(row_number()))) %>% 
    dplyr::ungroup() -> gtf
}

if(! "transcript_name" %in% colnames(gtf)) {
  gtf %>% 
    dplyr::mutate(transcript_name = transcript_id) -> gtf
} 
if(! "gene_name" %in% colnames(gtf)) {
  gtf %>% 
    dplyr::mutate(gene_name = gene_id) -> gtf
} 

if(! "transcript_support_level" %in% colnames(gtf)) {
  gtf %>% 
    tibble::add_column(transcript_support_level = NA) -> gtf
} 

if(!"ccdsid" %in% colnames(gtf)) {
  if("ccds_id" %in% colnames(gtf)) {
    gtf %>%
      dplyr::rename(ccdsid = ccds_id) -> gtf
  } else {
    gtf %>%
      tibble::add_column(ccdsid = NA) -> gtf
  }
}

if(!"gene" %in% gtf$type) {
  gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::group_by(seqnames, gene_name, gene_id, strand) %>%
    dplyr::summarise(start = min(start), end = min(end), .groups = "drop") %>%
    tibble::add_column(type = "gene") %>%
    dplyr::mutate(width = 1+end-start) %>%
    dplyr::bind_rows(gtf) -> gtf
}

if(!"transcript" %in% gtf$type) {
  gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::group_by(seqnames, gene_name, gene_id, transcript_name, transcript_id, strand, ccdsid, transcript_support_level) %>%
    dplyr::summarise(start = min(start), end = min(end), .groups = "drop") %>%
    tibble::add_column(type = "transcript") %>%
    dplyr::mutate(width = 1+end-start) %>%
    dplyr::bind_rows(gtf) -> gtf
}


gtf %>%
  dplyr::select(seqnames, start, end, strand, type, gene_name, gene_id, transcript_id, exon_number, width, transcript_name, ccdsid, transcript_support_level) %>%
  dplyr::filter(type %in% c("gene", "transcript", "exon")) -> gtf


################################################################
# Identify nexons output files and genes within them
################################################################

print(Sys.time())
print("Identifying nexons output files and genes; generating short exon/intron files...")

list.files(nexons_directory, "txt$", full.names = T) -> all_files_to_process

if(length(all_files_to_process) == 0) { stop(paste0("No txt files found in nexons output directory (", nexons_directory, ")!")) }


all_files_to_process %>%
  furrr::future_map(~ {
    read_tsv(.x, show_col_types = F) %>% distinct(GeneName) %>% pull() -> gene
    if(length(gene) == 0) { NA_character_ } else { gene }
    }) -> all_genes_to_analyse

all_genes_to_analyse %>%
  purrr::map_lgl(is.na) %>%
  which() -> exclude_files

if(length(exclude_files) > 0) {
  print("Excluding files (no data or unable to read):")
  print(all_files_to_process[exclude_files])
  all_files_to_process <- all_files_to_process[-exclude_files]
}

all_genes_to_analyse %>%
  unlist() %>%
  na.omit() %>%
  unique() -> all_genes_to_analyse

if(length(all_genes_to_analyse) == 0) { stop(paste0("No readable txt files containing data found in nexons output directory (", nexons_directory, ")!")) }

#######################################################################
# Overlapping genes file
#######################################################################

gtf %>%
  dplyr::filter(type == "gene") %>%
  with(GRanges(seqnames, IRanges(start, end), strand, gene_name, gene_id)) -> genes_gr

findOverlaps(genes_gr, genes_gr, maxgap = overlapGap) %>% 
  tidyr::as_tibble() %>%
  dplyr::mutate(
    gene1 = genes_gr$gene_name[queryHits],
    gene2 = genes_gr$gene_name[subjectHits],
    gene2_id = genes_gr$gene_id[subjectHits]
  ) %>%
  dplyr::select(starts_with("gene")) %>%
  dplyr::filter(gene1 %in% all_genes_to_analyse) %>%
  dplyr::filter(gene1 != gene2) -> all_overlap_genes


######################################################################
# Identify combinations of short exons that could have been missed
######################################################################

gtf %>% 
  dplyr::filter(gene_name %in% all_genes_to_analyse & type == "exon") %>% 
  dplyr::group_by(transcript_id) %>%
  dplyr::filter(n() > min_exons) %>% # since if it only has min_exons, a splice pattern missing 1 or more exons won't have been detected
  dplyr::mutate(exon_coord = case_when(
    exon_number == 1 & strand  == "+" ~ as.character(end),
    exon_number == 1 ~ as.character(start),
    exon_number == max(exon_number) & strand == "+" ~ as.character(start),
    exon_number == max(exon_number) ~ as.character(end),
    strand == "+" ~ paste(start, end, sep = "-"),
    TRUE ~ paste(end, start, sep = "-")
  )) %>%
  dplyr::ungroup() %>%
  dplyr::filter(width <= missedExonWidth) %>%
  dplyr::arrange(transcript_name, width) %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::slice(1:5) %>% # if > 5 short exons, only include shortest 5 (or first 5 if the same length) - otherwise gets exponentially (n^n) more memory intensive due to possible exon combinations (few genes exceed this; have checked in our data none would be picked up as missing >5 anyway)
  ungroup() %>%
  tidyr::nest(exons = -transcript_id) %>%
  dplyr::mutate(exon_combinations = furrr::future_map(exons, ~ {
    .x %>%
      dplyr::pull(exon_coord) -> all_coord
    purrr::map(seq_along(1:nrow(.x)), function(y) {
      combn(all_coord, y, simplify = F)
    }) %>%
      unlist(recursive = F)
  })) -> short_exon_combinations


######################################################################
# Add in short introns from annotated isoforms that will be missing
######################################################################

gtf %>%
  dplyr::filter(gene_name %in% all_genes_to_analyse & type == "exon") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(exon_number) %>%
  dplyr::mutate(start_junction = if_else(strand == "+", end, start)) %>%
  dplyr::mutate(end_junction = if_else(strand == "+", dplyr::lead(start), dplyr::lead(end))) %>%
  dplyr::mutate(intron_width = abs(end_junction - start_junction)-1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(intron_width < minIntronWidth) %>%
  tidyr::unite(intron_coord, start_junction, end_junction, sep = ":") %>%
  dplyr::select(transcript_id, intron_coord, intron_width) %>%
  tidyr::nest(introns = c(intron_coord, intron_width)) %>%
  dplyr::full_join(short_exon_combinations) %>% 
  dplyr::select(-exons) -> short_exons_introns

#################################################################################
# Create the splice patterns missing each combination of each short exon/intron
#################################################################################

if(nrow(short_exons_introns) > 0) {
  
  # collapse any splice patterns that could match more than one annotated isoform (with different combination of missing exons)
  
  gtf %>%
    dplyr::filter(transcript_id %in% short_exons_introns$transcript_id & type == "exon") %>%
    dplyr::arrange(transcript_name, exon_number) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(exon_coord = case_when(
      exon_number == 1 & strand  == "+" ~ as.character(end),
      exon_number == 1 ~ as.character(start),
      exon_number == max(exon_number) & strand == "+" ~ as.character(start),
      exon_number == max(exon_number) ~ as.character(end),
      strand == "+" ~ paste(start, end, sep = "-"),
      TRUE ~ paste(end, start, sep = "-")
    )) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_name, transcript_id, exon_coord, width) %>%
    tidyr::nest(exons = c(exon_coord, width)) %>% 
    dplyr::left_join(short_exons_introns) %>%
    dplyr::mutate(alternative_splice_patterns = furrr::future_map2(exons, exon_combinations, ~ { # first deal with the exons that may have been missed. Then remove introns (otherwise will potentially disrupt exon coord)
      .y %>%
        append("") %>% # will retain whole splice pattern, for removal of just introns (can remove later for those without short introns)
        as.list() %>% # for those with only short introns
        tibble::as_tibble_col("missing_exons") %>%
        dplyr::mutate(missing_exon_widths = purrr::map(missing_exons, function(z) .x$width[match(z, .x$exon_coord)])) %>%
        dplyr::mutate(splice_pattern = purrr::map_chr(missing_exons, function(z) {
          .x %>%
            dplyr::filter(!exon_coord %in% z) %>%
            dplyr::pull(exon_coord) %>%
            paste(collapse = ":")
        })) %>%
        dplyr::mutate(first_exon = purrr::map_lgl(missing_exons, function(z) 1 %in% which(.x$exon_coord %in% z))) %>%
        dplyr::mutate(last_exon = purrr::map_lgl(missing_exons, function(z) nrow(.x) %in% which(.x$exon_coord %in% z))) %>%
        dplyr::mutate(splice_pattern = sub("-[0-9]+$", "", sub("^[0-9]+-", "", splice_pattern))) %>% # if 1 exon permitted, it's the end coordinate that will be retained so needs to be substituted in this order
        dplyr::mutate(missing_exons = purrr::map_chr(missing_exons, function(z) paste(z, collapse = ";"))) %>%
        dplyr::mutate(missing_exon_widths = purrr::map_chr(missing_exon_widths, function(z) paste(z, collapse = ";")))
    })) %>% 
    dplyr::mutate(alternative_splice_patterns = furrr::future_map2(alternative_splice_patterns, introns, ~ {
      if(length(.y) == 0) {
        .x %>%
          dplyr::filter(!is.na(missing_exons)) %>%
          dplyr::mutate(missing_introns = NA_character_) %>%
          dplyr::mutate(missing_intron_widths = NA_character_)
      } else {
        .x %>%
          dplyr::mutate(missing_introns = list(.y$intron_coord)) -> new_tibble
        for(z in .y$intron_coord) { 
          new_tibble %>% dplyr::mutate(missing_introns = dplyr::if_else(grepl(z, splice_pattern), missing_introns, purrr::map(missing_introns, function(a) {
            a[!a == z]
          }))) %>%
            dplyr::mutate(splice_pattern = sub("--", "-", sub(z, "", splice_pattern))) -> new_tibble
        }
        new_tibble %>%
          dplyr::mutate(splice_pattern = gsub("-$|^-", "", splice_pattern)) %>%
          dplyr::mutate(missing_intron_widths = purrr::map_chr(missing_introns, function(z) {
            if(length(z) == 0) { "" } else {
              .y %>%
                dplyr::filter(intron_coord %in% z) %>%
                dplyr::pull(intron_width) %>%
                paste(collapse = ";")
            }
          })) %>%
          dplyr::mutate(missing_introns = purrr::map_chr(missing_introns, function(z) paste(z, collapse = ";"))) -> new_tibble
        return(new_tibble)
      }
    })) %>%
    tidyr::unnest(alternative_splice_patterns) %>%
    dplyr::mutate(splice_pattern = dplyr::if_else(nchar(splice_pattern) > 0, splice_pattern, furrr::future_map_chr(exons, ~ {
      .x %>% 
        dplyr::slice(nrow(.x)) %>% 
        dplyr::mutate(exon_end = sub("^[0-9]+-", "", exon_coord)) %>%
        dplyr::pull(exon_end)
    }))) %>% # annotated end of the last exon (ie, how single exon expected to appear)
    dplyr::select(-exons, -introns, -exon_combinations) %>% 
    dplyr::mutate(includes_firstlast_exon = dplyr::if_any(ends_with("exon"))) %>%
    dplyr::mutate(dplyr::across(starts_with("missing"), ~ dplyr::if_else(grepl("[0-9]", .x), .x, NA_character_))) %>%
    dplyr::filter(!(is.na(missing_exons) & is.na(missing_introns))) %>%
    dplyr::group_by(gene_name, splice_pattern) %>% 
    dplyr::summarise(
      MissingExonIntron_TranscriptID = paste(transcript_id, collapse = "_"),
      MissingExons = paste(missing_exons, collapse = "_"), 
      MissingExonWidths = paste(missing_exon_widths, collapse = "_"),
      MissingFirstLastExon = paste(includes_firstlast_exon, collapse = "_"),
      MissingIntrons = paste(missing_introns, collapse = "_"), 
      MissingIntronWidths = paste(missing_intron_widths, collapse = "_"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(matches("^MissingEx|^MissingInt"), ~ dplyr::if_else(grepl("[0-9]", .x), .x, NA_character_))) -> missing_exon_intron_splicePatterns
  
}
  
  


with_progress({
  p <- progressor(steps = 6*length(all_files_to_process))
  furrr::future_map(all_files_to_process, function(txtFile) {
    p(message = paste0("Processing ", txtFile, "; ", Sys.time()), class = "sticky")
    txtFile %>%
      readr::read_tsv(col_types = cols(Chr = col_character(), Count = col_double(), Starts = col_character(), Ends = col_character(), SplicePattern = col_character())) %>%
      dplyr::rename_with(~ gsub(" ", "", .x)) -> nexons_output_all
    
    #############################################################
    # Remove introns < minIntronWidth
    #############################################################
    
    p(message = paste0("Short introns/exons: ", txtFile, "; ", Sys.time()), class = "sticky")
    
    nexons_output_all %>%
      dplyr::mutate(introns = stringr::str_split(SplicePattern, "-")) %>% 
      dplyr::mutate(introns = map(introns, ~ {
        tibble::as_tibble_col(.x, "coord") %>%
          tidyr::separate(coord, c("start", "end"), convert = T, remove = F) %>%
          dplyr::mutate(width = abs(end-start)-1) %>%
          dplyr::filter(width < minIntronWidth)
      })) %>%
      dplyr::mutate(nShortIntrons = purrr::map_dbl(introns, nrow)) %>%
      dplyr::mutate(ShortIntrons = purrr::map_chr(introns, ~ paste(.x$coord, collapse = "_"))) %>%
      dplyr::mutate(ShortIntronWidths = purrr::map_chr(introns, ~ paste(.x$width, collapse = "_"))) %>%
      dplyr::mutate(minShortIntron = purrr::map_dbl(introns, ~ min(.x$width))) %>%
      dplyr::mutate(SplicePattern_new = purrr::map2_chr(introns, SplicePattern, ~ {
        if(nrow(.x) == 0) {
          .y 
        } else {
          excluded_pattern <- .y
          for(i in 1:nrow(.x)) {
            excluded_pattern <- sub("--", "-", sub(.x$coord[i], "", excluded_pattern))
          }
          excluded_pattern <- gsub("-$|^-", "", excluded_pattern) # Removes trailing dash if first or last junction removed
          if(nchar(excluded_pattern) == 0) { "single_exon" } else { # don't have access to ends here; below replace this with the mode of the ends (closest approximation to how single exon appear, although this could be a combination of multiple)
            excluded_pattern
          }
        }
      })) %>% 
      group_by(GeneID, SplicePattern) %>%
      dplyr::mutate(SplicePattern_new = if_else(SplicePattern_new != "single_exon", SplicePattern_new, {
        stringr::str_split(Ends, ",") %>%
          unlist() %>%
          tibble::as_tibble_col("End") %>%
          dplyr::group_by(End) %>%
          dplyr::count() %>%
          dplyr::ungroup() %>%
          dplyr::arrange(desc(n)) %>%
          dplyr::slice(1) %>% dplyr::pull(End)
      })) %>%
      group_by(File, GeneID, SplicePattern) %>%
      dplyr::mutate(ShortIntrons = if_else(ShortIntrons == "", NA_character_, ShortIntrons)) %>%
      dplyr::mutate(ShortIntronWidths = if_else(ShortIntronWidths == "", NA_character_, ShortIntronWidths)) %>%
      dplyr::mutate(minShortIntron = if_else(is.infinite(minShortIntron), NA_real_, minShortIntron)) %>%
      dplyr::select(-introns) %>%
      dplyr::ungroup() -> nexons_output_refined
    
    
    
    ###############################################################################################
    # Reduce to merge splice patterns that are now identical - or that are within the flexibility
    ###############################################################################################
    
    nexons_output_refined %>% 
      dplyr::arrange(nShortIntrons, desc(Count)) -> nexons_output_ordered
    
    nexons_output_ordered %>%
      dplyr::filter(nShortIntrons == 0) -> nexons_output_refined # these can be accepted immediately
    
    nexons_output_ordered %>%
      dplyr::filter(nShortIntrons > 0) -> nexons_output_ordered
    
    if(nrow(nexons_output_ordered) > 0) {
      nexons_output_refined %>%
        dplyr::distinct(SplicePattern_new) %>%
        dplyr::pull() %>%
        split(paste0("j", stringr::str_count(., ":"))) -> all_splices 
      all_splices %>% 
        purrr::map(~ stringr::str_split(.x, "[-:]") %>% purrr::map(as.integer)) -> all_splice_coord
    
      for(i in 1:nrow(nexons_output_ordered)) {
        test_splice <- nexons_output_ordered$SplicePattern_new[i]
        test_splice %>%
          stringr::str_split("[-:]") %>% 
          purrr::map(as.integer) -> test_splice_coord
        paste0("j", stringr::str_count(test_splice, ":")) -> test_index
        
        if(test_index %in% names(all_splice_coord)) {
          purrr::map(all_splice_coord[[test_index]], ~ all(abs(test_splice_coord[[1]] - .x) < flexibility)) %>%
            unlist() %>%
            which() -> matching_splice
          if(length(matching_splice) > 0) {
            nexons_output_ordered$SplicePattern_new[i] <- all_splices[[test_index]][matching_splice][1]
          } else {
            all_splices[[test_index]] <- c(all_splices[[test_index]], test_splice)
            all_splice_coord[[test_index]] <- c(all_splice_coord[[test_index]], test_splice_coord)
          }
        } else {
          all_splices[[test_index]] <- test_splice
          all_splice_coord[[test_index]] <- test_splice_coord
        }
        
        
        
      }
      
      nexons_output_refined %>%
        dplyr::bind_rows(nexons_output_ordered) -> nexons_output_refined
    }
    nexons_output_refined %>%
      dplyr::group_by(GeneID, GeneName, Chr, Strand, SplicePattern_new) %>% # need to replace transcript ID without grouping by file
      dplyr::mutate(TranscriptID = TranscriptID[1]) %>% 
      dplyr::group_by(File, GeneID, GeneName, Chr, Strand, SplicePattern_new, TranscriptID) %>% 
      dplyr::summarise(
        Count_new = sum(Count),
        Starts = paste(Starts, collapse = ","),
        Ends = paste(Ends, collapse = ","),
        MergedIsoforms = MergedIsoforms[1],
        shortIntron_nIsoform = sum(!is.na(ShortIntrons)),
        shortIntron_count = sum(Count[!is.na(ShortIntrons)]),
        shortIntron_nUnique = {
          ShortIntrons %>%
            stringr::str_split("_") %>%
            unlist() %>%
            n_distinct(na.rm = T)
        },
        shortIntron_minWidth = {
          ShortIntronWidths %>%
            stringr::str_split("_") %>%
            unlist() %>%
            as.integer() %>%
            min(na.rm = T)
        },
        shortIntron_maxWidth = {
          ShortIntronWidths %>%
            stringr::str_split("_") %>%
            unlist() %>%
            as.integer() %>%
            max(na.rm = T)
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(shortIntron_minWidth = if_else(is.finite(shortIntron_minWidth), shortIntron_minWidth, NA_real_)) %>%
      dplyr::mutate(shortIntron_maxWidth = if_else(is.finite(shortIntron_maxWidth), shortIntron_maxWidth, NA_real_)) %>%
      dplyr::rename(SplicePattern = SplicePattern_new, Count = Count_new) -> nexons_output_refined
    
    
    
    
    
      
    #################################################################################
    # Join annotated short exon/intron table (generated at start) to nexons results
    #################################################################################
      
    if(nrow(short_exons_introns) > 0) {
      
      nexons_output_refined %>%
        dplyr::left_join(missing_exon_intron_splicePatterns, by = c(GeneName = "gene_name", SplicePattern = "splice_pattern")) -> nexons_output_refined
      
    } else {
      nexons_output_refined %>% 
        tibble::add_column(
          MissingExonIntron_TranscriptID = NA_character_,
          MissingExons = NA_character_,
          MissingExonWidths = NA_character_,
          MissingFirstLastExon = NA_character_,
          MissingIntrons = NA_character_,
          MissingIntronWidths = NA_character_
        ) -> nexons_output_refined
      print(Sys.time())
      print(paste("No annotated exons/introns below specified widths identified for", txtFile))
      
    }
    
    
    
    ###############################################################################
    # Identify genes that overlap (or are within overlapGap) on the same strand
    ###############################################################################
    
    p(message = paste0("Overlapping genes: ", txtFile, "; ", Sys.time()), class = "sticky")
    
    all_overlap_genes %>%
      dplyr::filter(gene1 %in% nexons_output_all$GeneName) -> overlapping_genes
    
    
    ###############################################################################################################################
    # Check pairs/groups of overlapping genes to identify any reads that are likely to originate from an overlapping gene
    #
    # Also, flag types of splice match, or for single exon, the gene feature(s) overlapped (could be extended to 3'/5'UTR)
    # This gives similar info to eg SQANTI
    ###############################################################################################################################

    nexons_output_refined %>%
      tidyr::nest(data = -c(GeneID, GeneName, Chr, Strand)) %>% 
      dplyr::mutate(data = purrr::pmap(list(GeneName, Chr, Strand, data), function(gene, chrm, strand, data) {
        gtf %>%
          dplyr::filter(gene_name == gene & type == "exon" & strand == strand & seqnames == chrm) %>%
          dplyr::group_by(transcript_id) %>%
          dplyr::mutate(
            acceptor = case_when(
              exon_number == 1 ~ NA_integer_,
              strand == "+" ~ start, 
              TRUE ~ end
            ),
            donor = case_when(
              exon_number == max(exon_number) ~ NA_integer_,
              strand == "+" ~ end, 
              TRUE ~ start
            )
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(acceptor, donor) %>%
          dplyr::distinct() -> gene_exons
        gtf %>%
          dplyr::filter(gene_name == gene & type == "gene" & strand == strand & seqnames == chrm) %>%
          with(
            GRanges(seqnames, IRanges(start, end), strand = strand)
            ) -> gene_GR
        gtf %>%
          dplyr::filter(gene_name == gene & type == "exon" & strand == strand & seqnames == chrm) %>%
          group_by(transcript_id) %>%
          group_split() %>%
          map(~ {
            .x %>%
              with(
                GRanges(seqnames, IRanges(start, end), strand = strand)
              )
          }) -> txList_exon_GR
        gtf %>%
          dplyr::filter(gene_name == gene & type == "exon" & strand == strand & seqnames == chrm) %>%
          distinct(seqnames, start, end, strand) %>%
          with(
            GRanges(seqnames, IRanges(start, end), strand = strand)
          ) -> exon_GR
        GenomicRanges::setdiff(gene_GR, exon_GR) -> intron_GR
        gtf %>%
          dplyr::filter(gene_name == gene & type == "transcript" & strand == strand & seqnames == chrm) %>%
          dplyr::group_by(gene_name, transcript_id) %>%
          dplyr::summarise(start = min(start), end = max(end), .groups = "drop") %>%
          dplyr::add_row(gene_name = "none", transcript_id = "none", start = 0, end = Inf) %>%
          dplyr::mutate(width = 1+end-start) %>%
          tidyr::unite("gene_transcript", gene_name, transcript_id, remove = F) %>%
          dplyr::arrange(width) -> gene_range
        if(gene %in% overlapping_genes$gene1) {
          gtf %>%
            dplyr::filter(gene_id %in% overlapping_genes$gene2_id[overlapping_genes$gene1 == gene] & type == "exon" & strand == strand & seqnames == chrm) %>%
            dplyr::group_by(transcript_id) %>%
            dplyr::mutate(
              acceptor = case_when(
                exon_number == 1 ~ NA_integer_,
                strand == "+" ~ start, 
                TRUE ~ end
              ),
              donor = case_when(
                exon_number == max(exon_number) ~ NA_integer_,
                strand == "+" ~ end, 
                TRUE ~ start
              )
            ) %>%
            dplyr::ungroup() %>%
            dplyr::filter(!(is.na(acceptor) & is.na(donor))) %>%
            dplyr::select(gene_name, transcript_id, exon_number, start, end, acceptor, donor) %>%
            dplyr::distinct() -> overlap_exons
          gtf %>%
            dplyr::filter(gene_id %in% overlapping_genes$gene2_id[overlapping_genes$gene1 == gene] & type == "transcript" & strand == strand & seqnames == chrm) %>%
            dplyr::group_by(seqnames, gene_name, transcript_id, strand) %>%
            dplyr::summarise(start = min(start), end = max(end), .groups = "drop") %>%
            dplyr::add_row(gene_name = "none", transcript_id = "none", start = 0, end = Inf) %>%
            dplyr::mutate(width = 1+end-start) %>%
            tidyr::unite("gene_transcript", gene_name, transcript_id, remove = F) %>%
            dplyr::arrange(width) -> overlap_range
        } else {
          tibble(
            seqnames = chrm, gene_name = "none", transcript_id = "none", strand = strand, gene_transcript = "none_none", start = 0, end = Inf, width = Inf
          ) -> overlap_range
          tibble(gene_name = character(), acceptor = integer(), donor = integer()) -> overlap_exons
        }
        overlap_range %>% 
          dplyr::filter(gene_name != "none") %>%
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) -> overlap_GR
        
        data %>%
          tidyr::nest(data2 = -c(SplicePattern, TranscriptID)) %>%
          dplyr::mutate(exons = stringr::str_split(paste0("NA-",SplicePattern, "-NA"), ":")) %>%
          dplyr::mutate(read_start = map_int(data2, ~ stringr::str_split(.x$Starts, ",") %>% unlist() %>% tibble::as_tibble_col() %>% dplyr::group_by(value) %>% count() %>% dplyr::ungroup() %>% dplyr::arrange(desc(n)) %>% dplyr::slice(1) %>% dplyr::pull(value) %>% as.integer())) %>%
          dplyr::mutate(read_end = map_int(data2, ~ stringr::str_split(.x$Ends, ",") %>% unlist() %>% tibble::as_tibble_col() %>% dplyr::group_by(value) %>% count() %>% dplyr::ungroup() %>% dplyr::arrange(desc(n)) %>% dplyr::slice(1) %>% dplyr::pull(value) %>% as.integer())) %>%
          dplyr::mutate(n_exons = purrr::map_int(exons, length)) -> data
        data %>% 
          dplyr::slice(0) %>%
          add_column(match_overlapping = list()) -> data_new
        for(i in 1:nrow(data)) {
          data_new %>%
            bind_rows(
              data %>%
                dplyr::slice(i) %>%
                dplyr::mutate(match_overlapping = purrr::map(exons, function(x) {
                  purrr::map(x, function(y) {
                    which(
                      purrr::map_lgl(overlap_exons$donor, function(z) abs(z-as.integer(sub("-NA", "", sub("^.*?-", "", y)))) < flexibility) | purrr::map_lgl(overlap_exons$acceptor, function(z) abs(z-as.integer(sub("NA-", "", sub("-.*?$", "", y)))) < flexibility)# 2nd sub is for single exon (will check this coord vs both donor and acceptor...neither should really match but to do it properly will need to consider single exon coordinates separately)
                    ) -> exon_indices
                    overlap_exons %>%
                      dplyr::slice(exon_indices) %>%
                      dplyr::distinct(gene_name) %>%
                      dplyr::pull(gene_name)
                  }) %>% unlist() %>% unique()
                }))
            ) -> data_new
        }
        data_new %>% 
          dplyr::mutate(match_overlapping = purrr::map_chr(match_overlapping, ~ if_else(length(.x) == 0, "none", paste(.x, collapse = "_")))) %>%
          dplyr::mutate(match_overlap_n = if_else(match_overlapping == "none", 0, purrr::map_int(exons, function(x) {
            purrr::map_lgl(x, function(y) {
              (
                (sum(purrr::map_lgl(unique(overlap_exons$donor), function(z) abs(z-as.integer(sub("-NA", "", sub("^.*?-", "", y)))) < flexibility), na.rm = T) > 0 | is.na(as.integer(sub("-NA", "", sub("^.*?-", "", y))))) &
                  (sum(purrr::map_lgl(unique(overlap_exons$acceptor), function(z) abs(z-as.integer(sub("NA-", "", sub("-.*?$", "", y)))) < flexibility), na.rm = T) > 0 | is.na(as.integer(sub("NA-", "", sub("-.*?$", "", y)))))
              )
            }) %>%
              sum()
          }))) %>% 
          dplyr::mutate(match_gene = purrr::map_int(exons, function(x) {
            if(nrow(gene_exons) < 2) { 0 } else {
              purrr::map_lgl(x, function(y) {
                (
                  (sum(purrr::map_lgl(unique(gene_exons$donor), function(z) abs(z-as.integer(sub("-NA", "", sub("^.*?-", "", y)))) < flexibility), na.rm = T) > 0 | is.na(as.integer(sub("-NA", "", sub("^.*?-", "", y))))) &
                    (sum(purrr::map_lgl(unique(gene_exons$acceptor), function(z) abs(z-as.integer(sub("NA-", "", sub("-.*?$", "", y)))) < flexibility), na.rm = T) > 0 | is.na(as.integer(sub("NA-", "", sub("-.*?$", "", y)))))
                )
              }) %>%
                sum()
            }  
          })) %>% 
          dplyr::mutate(range_gene_tx = purrr::map2_chr(exons, read_start, function(x, y) {
            if(length(x) > 1) {
              purrr::map2_lgl(gene_range$start, gene_range$end, function(start, end) dplyr::between(as.integer(sub("NA-", "", x[1])), start, end) & dplyr::between(as.integer(sub("-NA", "", x[length(x)])), start, end)) %>% # splice junctions should all be within gene coord
                which() -> gene_indices
            } else {
              testend <- as.integer(sub("-NA", "", sub("NA-", "", x[1])))
              purrr::map2_lgl(gene_range$start, gene_range$end, function(start, end) dplyr::between(y, start-250, end+250) & dplyr::between(testend, start-250, end+250)) %>% # for single exon give some flexibility
                which() -> gene_indices
            }
            gene_range %>%
              dplyr::slice(gene_indices) %>%
              dplyr::slice(1) %>%
              dplyr::pull(gene_transcript)
          })) %>%
          dplyr::mutate(range_overlap_tx = purrr::map2_chr(exons, read_start, function(x, y) {
            if(length(x) > 1) {
              purrr::map2_lgl(overlap_range$start, overlap_range$end, function(start, end) dplyr::between(as.integer(sub("NA-", "", x[1])), start, end) & dplyr::between(as.integer(sub("-NA", "", x[length(x)])), start, end)) %>%
                which() -> overlap_indices
            } else {
              testend <- as.integer(sub("-NA", "", sub("NA-", "", x[1])))
              purrr::map2_lgl(overlap_range$start, overlap_range$end, function(start, end) dplyr::between(y, start-250, end+250) & dplyr::between(testend, start-250, end+250)) %>%
                which() -> overlap_indices
            }
            overlap_range %>%
              dplyr::slice(overlap_indices) %>%
              dplyr::slice(1) %>%
              dplyr::pull(gene_transcript)
          })) %>% 
          dplyr::mutate(range_overlapping = overlap_range$gene_name[match(range_overlap_tx, overlap_range$gene_transcript)]) %>%
          dplyr::mutate(range_gene = gene_range$gene_name[match(range_gene_tx, gene_range$gene_transcript)]) %>%
          dplyr::mutate(overlap_read = map2_chr(read_start, read_end, ~ {
            overlap_GR %>%
              subsetByOverlaps(GRanges(chrm, IRanges(min(.x, .y), max(.x, .y)), strand = strand)) %>% 
              as_tibble() %>% 
              distinct(gene_name, .keep_all = T) %>% 
              select(gene_name, transcript_id) -> overlapping_read
            if(nrow(overlapping_read) > 0) {
              overlapping_read %>%
                pivot_wider(names_from = transcript_id, values_from = gene_name) %>% 
                unite("overlap_genes", everything()) %>% 
                pull()
            } else { "none" }
          })) %>%
          dplyr::group_by(SplicePattern) %>% 
          dplyr::mutate(Overlapping = dplyr::case_when(
            !gene %in% overlapping_genes$gene1 ~ "none",
            overlap_read == "none" ~ "none",
            match_gene > 0 & match_overlap_n == 0 ~ "none",
            match_gene > match_overlap_n ~ "possible_overlapping",
            match_overlap_n > match_gene ~ "likely_overlapping",
            range_gene != "none" & range_overlapping == "none" ~ "none",
            range_gene == "none" & range_overlapping != "none" & grepl("^Variant", TranscriptID) ~ "likely_overlapping",
            range_gene != "none" & gene_range$width[gene_range$gene_transcript == range_gene_tx] <= overlap_range$width[overlap_range$gene_transcript == range_overlap_tx] ~ "possible_overlapping",
            range_gene != "none" & grepl("^Variant", TranscriptID) ~ "likely_overlapping",
            TRUE ~ "possible_overlapping"
          )) %>%
          dplyr::mutate(OverlappingGene = dplyr::case_when(
            Overlapping == "none" ~ "none",
            match_overlapping == "none" & range_overlapping == "none" ~ overlap_read,
            match_overlapping == "none" ~ range_overlapping,
            range_overlapping == "none" ~ match_overlapping,
            grepl(range_overlapping, match_overlapping) ~ range_overlapping, # if a choice that have matching coord, want the (shortest) gene that contains the full splice pattern
            TRUE ~ paste(match_overlapping, range_overlapping, sep = "_")
          )) %>% 
          dplyr::mutate(SpliceMatch = case_when(
            n_exons == 1 ~ "single_exon",
            match_gene == 0 ~ "none",
            match_gene == n_exons ~ "all",
            match_gene > 0 ~ "partial"
          )) %>%
          dplyr::mutate(ReadLocation = purrr::map2_chr(read_start, read_end, function(x, y) {
            GRanges(chrm, IRanges(min(x,y), max(x,y)), strand = strand) -> read_GR
            GRanges(chrm, IRanges(c(x)), strand = strand) -> readStart_GR
            GRanges(chrm, IRanges(c(y)), strand = strand) -> readEnd_GR
            if(sum(purrr::map_lgl(txList_exon_GR, function(txGR) {
              IRanges::countOverlaps(read_GR, txGR) == length(txGR) & IRanges::overlapsAny(readStart_GR, txGR) & IRanges::overlapsAny(readEnd_GR, txGR)
            })) > 0) {
              if(sum(purrr::map_lgl(txList_exon_GR, function(txGR) {
                IRanges::countOverlaps(read_GR, txGR, type = "within") == 1 & length(txGR) == 1
              })) > 0) {
                "complete_tx_singleExon" # does not necessarily mean the read is single exon, just that it overlaps and is fully within an annotated single exon transcript (so don't flag as RI isoform)
              } else {
                "complete_tx"
              }
            } else if(IRanges::overlapsAny(read_GR, exon_GR, type = "within")) {
              "exon"
            } else if(IRanges::overlapsAny(read_GR, exon_GR) & IRanges::overlapsAny(read_GR, gene_GR, type = "within")) {
              if(IRanges::overlapsAny(readStart_GR, intron_GR) & IRanges::overlapsAny(readEnd_GR, intron_GR)) {
                "intron_exon_intron"
              } else if(IRanges::overlapsAny(readStart_GR, exon_GR) & IRanges::overlapsAny(readEnd_GR, exon_GR)) {
                "exon_intron_exon"
              } else {
                "exon_intron"
              }
            } else if(IRanges::overlapsAny(read_GR, gene_GR, type = "within")) {
              "intron"
            } else if(IRanges::overlapsAny(read_GR, gene_GR) & IRanges::overlapsAny(read_GR, intron_GR)) {
              "exon_intron_outside"
            } else if(IRanges::overlapsAny(read_GR, gene_GR)) {
              "exon_outside"
            } else { "outside" }
          })) %>% 
          dplyr::ungroup() %>%
          dplyr::select(-exons, -n_exons, -starts_with("match"), -starts_with("range"), -read_start, -read_end, -overlap_read) %>%
          tidyr::unnest(data2)
      })) %>%
      tidyr::unnest(data) %>%
      select(File, everything()) -> nexons_output_refined
    
    
    ###############################################################################
    # Identify potential truncations
    ###############################################################################
    
    # needs to be done on one gene at a time - the way nexons works at the moment, each txt file only contains
    # 1 gene, so don't need to iterate, but in case that changes, will filter for one gene at a time.
    
    p(message = paste0("Identifying potential truncations: ", txtFile, "; ", Sys.time()), class = "sticky")
    
    purrr::map(unique(nexons_output_refined$GeneName), function(gene) {
      nexons_output_refined %>% 
        dplyr::filter(GeneName == gene) %>%
        identifyPotentialTruncations(gtf = gtf, flexibility = flexibility)
    }) %>%
      dplyr::bind_rows() -> nexons_output_refined
    
    
    
    ###############################################################################
    # Find commonly occurring start/ends of transcripts
    ###############################################################################
    
    p(message = paste0("Start/end sites: ", txtFile, "; ", Sys.time()), class = "sticky")
    
    nexons_output_refined %>%
      dplyr::mutate(
        Starts = stringr::str_split(Starts, ","),
        Ends = stringr::str_split(Ends, ",")
      ) %>%
      tidyr::nest(geneData = -GeneName) %>%
      dplyr::mutate(reduced_starts = purrr::map2(geneData, GeneName, ~ {
        unlist(.x$Starts) %>%
          as.integer() -> starts
        starts %>%
          density(bw = 2.5, n=max(20,(1+max(starts)-min(starts)))) %>%
          broom::tidy() -> sDens
        sDens %>%
          dplyr::pull(y) %>%
          diff() %>% 
          sign() %>% 
          diff() -> sDiff
        
        sDens %>%
          dplyr::slice(which(sDiff < 0)+1) %>% 
          dplyr::filter(y > 0.0001) %>%
          dplyr::mutate(x = round(x)) %>%
          dplyr::distinct(x) %>%
          dplyr::group_by(x) %>%
          dplyr::mutate(coord = seq(x-2,x+2)[which.max(map(seq(x-2,x+2), function(y) sum(starts == y)))]) %>%
          dplyr::mutate(count_ext = sum(dplyr::between(starts, coord-2, coord+2)), count = sum(starts == coord)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(count_ext >= start_end_minreads) %>%
          dplyr::select(coord, count) -> found_starts
        gtf %>%
          dplyr::filter(gene_name == .y) %>%
          dplyr::filter(type == "transcript") %>%
          dplyr::arrange(is.na(ccdsid), as.numeric(transcript_support_level)) %>%
          dplyr::mutate(coord = if_else(strand == "+", start, end)) %>%
          dplyr::select(coord) %>%
          tibble::add_column(annotated = T) %>%
          dplyr::group_by(coord) %>%
          dplyr::mutate(count = sum(starts == coord)) %>%
          dplyr::ungroup() -> anno_starts
        found_starts %>% 
          dplyr::full_join(anno_starts, by = c("coord", "count")) %>% 
          dplyr::group_by(coord) %>%
          dplyr::mutate(count = if_else(is.na(annotated), count, max(10, 1.3*count))) %>% # an unannotated start will still be prioritised over annotated if it has at least 30% more supporting reads (and a minimum of 10 reads)
          dplyr::ungroup() %>% 
          dplyr::arrange(desc(count)) %>%
          dplyr::distinct(coord) %>%
          dplyr::pull() -> final_starts
        for(i in 1:length(final_starts)) {
          if(i > 1) {
            if(sum(abs(final_starts[i] - final_starts[seq(1, i-1)]) <= flexibility_ends, na.rm = T) > 0) {
              final_starts[i] <- NA
            }
          }
        }
        final_starts[!is.na(final_starts)]
      })) %>% 
      dplyr::mutate(reduced_ends = purrr::map2(geneData, GeneName, ~ {
        unlist(.x$Ends) %>%
          as.integer() -> ends
        ends %>%
          density(bw = 2.5, n = max(20,(1+max(ends)-min(ends)))) %>%
          broom::tidy() -> eDens
        eDens %>%
          dplyr::pull(y) %>%
          diff() %>% 
          sign() %>% 
          diff() -> eDiff
        
        eDens %>%
          dplyr::slice(which(eDiff < 0)+1) %>% 
          dplyr::filter(y > 0.0001) %>%
          dplyr::mutate(x = round(x)) %>%
          dplyr::distinct(x) %>%
          dplyr::group_by(x) %>%
          dplyr::mutate(coord = seq(x-2,x+2)[which.max(map(seq(x-2,x+2), function(y) sum(ends == y)))]) %>%
          dplyr::mutate(count_ext = sum(dplyr::between(ends, coord-2, coord+2)), count = sum(ends == coord)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(count_ext >= start_end_minreads) %>%
          dplyr::select(coord, count) -> found_ends
        gtf %>%
          dplyr::filter(gene_name == .y) %>%
          dplyr::filter(type == "transcript") %>%
          dplyr::arrange(is.na(ccdsid), as.numeric(transcript_support_level)) %>%
          dplyr::mutate(coord = if_else(strand == "+", end, start)) %>%
          dplyr::select(coord) %>%
          tibble::add_column(annotated = T) %>%
          group_by(coord) %>%
          mutate(count = sum(ends == coord)) %>%
          ungroup() -> anno_ends
        found_ends %>% 
          dplyr::full_join(anno_ends, by = c("coord", "count")) %>% 
          dplyr::group_by(coord) %>%
          dplyr::mutate(count = if_else(is.na(annotated), count, max(10, 1.3*count))) %>%
          dplyr::ungroup() %>% 
          dplyr::arrange(desc(count)) %>%
          dplyr::distinct(coord) %>%
          dplyr::pull() -> final_ends
        for(i in 1:length(final_ends)) {
          if(i > 1) {
            if(sum(abs(final_ends[i] - final_ends[seq(1, i-1)]) <= flexibility_ends, na.rm = T) > 0) {
              final_ends[i] <- NA
            }
          }
        }
        final_ends[!is.na(final_ends)]
      })) -> all_res_reduced_ends
    
    
    ########################################################################################################################################
    # Assign each start/end to closest common start/end (or unclear if > flexibility_ends2 away); collate counts for each start/end combo
    ########################################################################################################################################
    
    
    
    all_res_reduced_ends %>%
      tidyr::unnest(geneData) %>% 
      dplyr::group_by(File, GeneName, SplicePattern) %>%
      dplyr::mutate(
        starts_ends = list(tibble(
          start = (as.numeric(unlist(Starts))),
          end = (as.numeric(unlist(Ends)))
        ) %>%
          dplyr::group_by(start, end) %>%
          dplyr::summarise(count = n(), .groups = "drop") %>%
          dplyr::mutate(start_new = unlist(reduced_starts)[map_int(start, function(y) which.min(abs(y - unlist(reduced_starts))))]) %>%
          dplyr::mutate(end_new = unlist(reduced_ends)[map_int(end, function(y) which.min(abs(y - unlist(reduced_ends))))]) %>%
          dplyr::mutate(
            start_dist = abs(start_new - start),
            end_dist = abs(end_new - end)
          ) %>%
          dplyr::mutate(
            Start = if_else(start_dist > flexibility_ends2, "unclear", as.character(start_new)),
            End = if_else(end_dist > flexibility_ends2, "unclear", as.character(end_new))
          )
        )) %>% 
      dplyr::select(-c(Starts, Ends, reduced_starts, reduced_ends, Count)) %>%
      tidyr::unnest(starts_ends) %>% 
      dplyr::ungroup() %>%
      dplyr::group_by(File, GeneName, SplicePattern, TranscriptID, Start, End) %>%
      dplyr::summarise(Count = sum(count), .groups = "drop") -> all_res_complete_starts_ends
    
    
    ########################################################################################################################################
    # Identify start/end combinations that are within flexibility from the start/end for the TranscriptID(s)
    ########################################################################################################################################
    
    
    all_res_complete_starts_ends %>%
      tidyr::nest(data = -c(GeneName, TranscriptID, Start, End)) %>%
      dplyr::mutate(TranscriptEndMatch = if_else(Start == "unclear" | End == "unclear" | grepl("^Variant", TranscriptID), "", pmap_chr(list(TranscriptID, Start, End), function(transcript,start, end) {
        read_start <- min(as.integer(start), as.integer(end))
        read_end <- max(as.integer(start), as.integer(end))
        gtf %>%
          dplyr::filter(type == "transcript" & transcript_id %in% unlist(str_split(transcript, ":"))) %>%
          dplyr::filter(dplyr::between(start, read_start-flexibility_ends2, read_start+flexibility_ends2)) %>%
          dplyr::filter(dplyr::between(end, read_end-flexibility_ends2, read_end+flexibility_ends2)) %>%
          dplyr::pull(transcript_id) %>%
          paste(collapse = ":")
      }))) %>%
      dplyr::mutate(TranscriptEndMatch = if_else(TranscriptEndMatch == "", NA_character_, TranscriptEndMatch)) %>%
      tidyr::unnest(data) %>%
      dplyr::group_by(File, GeneName, SplicePattern, TranscriptID) %>%
      dplyr::mutate(proportion = Count/sum(Count)) -> all_res_complete_starts_ends
    
    
    
    ########################################################################################################################################
    # Collate common/annotated start/end combinations onto a single row per splice pattern; join to overall output
    ########################################################################################################################################
    
    
    
    all_res_complete_starts_ends %>%
      dplyr::filter(!is.na(TranscriptEndMatch)) %>%
      dplyr::group_by(File, GeneName, SplicePattern, TranscriptID, TranscriptEndMatch) %>%
      dplyr::summarise(proportion = signif(sum(proportion), 3)) %>%
      tidyr::unite("TranscriptProportions", TranscriptEndMatch, proportion, sep = ":") %>%
      dplyr::summarise(TranscriptProportions = paste(TranscriptProportions, collapse = "_"), .groups = "drop") -> ends_transcripts
    
    all_res_complete_starts_ends %>% 
      dplyr::filter(proportion > start_end_proportion) %>%
      dplyr::arrange(desc(proportion)) %>%
      dplyr::mutate(proportion = signif(proportion, 3)) %>%
      tidyr::unite("Start_End", Start, End, sep = "-", remove = F) %>%
      tidyr::unite("Start_End", Start_End, proportion, sep = ":", remove = F) %>%
      dplyr::left_join(distinct(nexons_output_refined, TranscriptID, TruncStartsAnno, TruncEndsAnno)) %>%
      dplyr::mutate(
        minDist_truncStart = purrr::map2_dbl(TruncStartsAnno, Start, ~ {
          min(abs(as.numeric(.y) - as.numeric(unlist(stringr::str_split(.x, "_")))))
        }),
        minDist_truncEnd = purrr::map2_dbl(TruncEndsAnno, End, ~ {
          min(abs(as.numeric(.y) - as.numeric(unlist(stringr::str_split(.x, "_")))))
        })
      ) %>% 
      dplyr::mutate(minDist_truncStart = if_else(is.na(minDist_truncStart), Inf, minDist_truncStart)) %>% # means will get this for non-truncations as well - does this fix NAs too?? - think it should as will give truncStart(/End) to those that don't have an annotated splice site following
      dplyr::mutate(minDist_truncEnd = if_else(is.na(minDist_truncEnd), Inf, minDist_truncEnd)) %>%
      dplyr::summarise(
        StartEndProportions = paste(Start_End, collapse ="_"), 
        DefinedStart = sum(proportion[Start != "unclear" & minDist_truncStart > 10]), # on basis that if the read contains a small portion of the next exon, it will map to the end (more or less) of the current exon -> may look like a defined start/end
        DefinedEnd = sum(proportion[End != "unclear" & minDist_truncEnd > 10]),
        .groups = "drop"
      ) %>%
      dplyr::left_join(ends_transcripts) -> all_res_common_starts_ends
    
    nexons_output_refined %>%
      dplyr::select(-TruncStartsAnno, -TruncEndsAnno) %>%
      dplyr::left_join(all_res_common_starts_ends) -> nexons_output_refined
    
    
    ########################################################################################################################################
    # Return output files as list (summary and complete ends)
    ########################################################################################################################################
    
    p(message = paste0("Finished processing ", txtFile, ";", Sys.time()), class = "sticky")
    return(list(nexons_output_refined, all_res_complete_starts_ends))
    
    
  }) -> postprocessing_output
})

postprocessing_output %>%
  map(~ .x[[1]]) %>%
  dplyr::bind_rows() %>%
  readr::write_tsv(output_file, col_names = F)

postprocessing_output %>%
  map(~ .x[[2]]) %>%
  dplyr::bind_rows() %>%
  readr::write_tsv(ends_outfile, col_names = F)


print(Sys.time())
print("Finished!")

write(FINISHED, FINISHED)
