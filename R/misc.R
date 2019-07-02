#' Add sequence column to GRanges object
#'
#' Add sequence column to GRanges object
#'
#' @param BSgenome BSgenome. BSgenome to query sequences, e.g. Hsapiens.
#' @param intervals GRanges. GRanges object containing coordinates to query.
#' @return Returns GRanges object containing seqs column
#' @export
addSeq <- function(BSgenome = NULL, 
                   intervals)
{
  if (is.null(BSgenome)) { stop ('Please provide valid BSgenome object') }
  # add genomic seq
  message('Adding sequences')
  seqlevelsStyle(intervals) <- 'UCSC'
  intervals <- as(intervals, 'GRanges')
  
  n_threads <- parallel::detectCores()
  
  if (n_threads > 5)
  {
    intvs_chunked <- 
      intervals %>% 
        mutate(cum_sum = cumsum(as.numeric(width)),
               chunks = cut(cum_sum, seq(1, max(cum_sum), l = n_threads + 1), labels = FALSE, include.lowest = TRUE)) %>%
        as_tibble
    
    seqs <- plyr::dlply(intvs_chunked, "chunks", 
                        .fun = function(x) getSeq(BSgenome, as_granges(x)), 
                        .parallel = TRUE)
     
    intervals$seqs <- unlist(DNAStringSetList(seqs))
  } else {
    intervals$seqs <- getSeq(BSgenome, intervals)
  }
  
  return(intervals)
}

#' Randomly mutate read sequences
#'
#' Introduces random mutations in read sequences given a specified error rate. Useful to simulate read mismapping under given error constrains.
#'
#' @param BSgenome BSgenome. BSgenome to query sequences, e.g. Hsapiens.
#' @param bam_path Character. Path to input BAM file.
#' @param outdir Character. Specifies where to store output fastq files. If NULL (the default), same directory as input bam file. 
#' @param tes GRanges. GRanges object containing TE coordinates
#' @param multi_mapper Logical. if FALSE (currently only supported), throws reads mapping to multiple position in the genome (NH:i tag).
#' @return Returns GRanges object containing seqs column
#' @export
mutateReads <- function(BSgenome = NULL,
                        reads,
                        tes = NULL,
                        paired = NULL,
                        n_reads = NULL,
                        outdir = NULL,
                        prefix = 'aligned_reads',
                        which = GRanges(),
                        error_rate = NULL,
                        seed = 19)
{
  if (is.null(BSgenome))    { stop('Please provide a valid BSgenome object, e.g. Hsapiens') }
  if (is.null(error_rate))  { stop ('Please provide error rate' ) }
  if (is.null(paired))      { stop ('Define whether BAM input contains single- or paired-end reads!') }
  if (is.null(outdir))      { outdir <- dirname(reads) }
  
  set.seed(seed)
  
  if (class(reads) != 'GRanges')
  {
    # import reads
    reads <- 
      importBAM(reads, 
                paired = paired,
                multi  = FALSE,
                what   = c('flag', 'qname', 'cigar', 'seq', 'qual'),
                anchor = 'fiveprime')
    
    names(reads) <- reads$qname                
  }

  # subset (by 1st mate 5' overlap) and annotate with TEs if provided
  if (!is.null(tes))
  {
    hits <- findOverlaps(reads %>% filter(first), tes)
    
    reads <- reads[names(reads) %in% names(reads %>% filter(first))[from(hits)]] # throw mates where first mate doesn't overlap TEs
  }
  
  # adjust read seq and qual (BAM stores 5'3' from + strand alignment and adjust qual accordingly)
  mcols(reads[strand(reads) == '-'])$seq  <- reverseComplement(mcols(reads[strand(reads) == '-'])$seq)
  mcols(reads[strand(reads) == '-'])$qual <- reverse(mcols(reads[strand(reads) == '-'])$qual)

  # split CIGAR string into intervals, label, and combine
  cigar_intvs_list   <- extractAlignmentRangesOnReference(reads$cigar, pos = start(reads))
  n_elements         <- elementNROWS(cigar_intvs_list)
  cigar_intvs        <- unlist(cigar_intvs_list)
  names(cigar_intvs) <- rep(names(reads), times =  n_elements)
  cigar_intvs        <- GRanges(ranges   = cigar_intvs,
                                seqnames = rep(as.character(seqnames(reads)), times = n_elements),
                                first    = rep(as.character(mcols(reads)$first), times = n_elements),
                                strand   = rep(as.character(strand(reads)), times = n_elements)) 
  
  # throw reads that exceed BSgenome (may happen because of CIGAR)
  reads_beyond     <- names(subsetByOverlaps(cigar_intvs, GRanges(seqinfo(BSgenome)), invert = TRUE)) 
  cigar_intvs_filt <- cigar_intvs[!names(cigar_intvs) %in% reads_beyond]
  message ('Threw ',  length(unique(reads_beyond)), ' reads exceeding BSgenome chrom boundaries')
  
  # sample reads if specified
  if (!is.null(n_reads))
  {
    n_reads_samp       <- n_reads
    n_reads            <- length(unique(names(cigar_intvs_filt)))
    n_reads_samp       <- ifelse(n_reads_samp > n_reads, n_reads, n_reads_samp)
    cigar_intvs_filt   <- cigar_intvs_filt[names(cigar_intvs_filt) %in% sample(unique(names(cigar_intvs_filt)), n_reads_samp)]
  }
  
  # add  seq to reads
  cigar_seqs  <- addSeq(BSgenome = BSgenome, cigar_intvs_filt)
  
  # concatenate CIGAR seqs per read (slow, can be improved)
  cigar_seqs_c <-  
    data.table(read_id = names(cigar_seqs),
               first   = as.logical(cigar_seqs$first),
               seqs    = as.character(cigar_seqs$seqs))[, .(seqs_c = paste0(seqs, collapse = '')), by = c('read_id', 'first')]
  
  # add phred quality and repname (very slow, can be improved)
  reads_sim <- 
    merge(cigar_seqs_c,
          data.table(read_id   = names(reads), 
                     first     = mcols(reads)$first, 
                     #coord     = mcols(reads)$id_unique,
                     qual      = as.character(mcols(reads)$qual)), all.x = TRUE, by = c('read_id', 'first'))
  
  # trim phred quality or seq to match in length. Seq now contains aligned only part, may extend/shorten based on indels of alignment.
  reads_sim[, qual := substring(qual, 1, nchar(seqs_c))]   # trim qual when longer than seq
  reads_sim[, seqs_c := substring(seqs_c, 1, nchar(qual))] # trim seqs when longer than phred
                     
  # introduce read errors and export
  for (errr in error_rate)
  {
    reads_sim$seqs_mut <- mutateSeqs(reads_sim$seqs_c, error_rate = errr, seed = seed)
    
    # bring to fastq format (consider make Util function)
    reads_sim$row_id      <- 1:nrow(reads_sim)
    reads_sim$dummy       <- '+'
    reads_sim$read_id_mod <- paste0('@', reads_sim$read_id)#, '_', reads_sim$id_unique)
     
    fastq_df         <- 
      melt(reads_sim[, c('read_id_mod', 'seqs_mut', 'dummy', 'qual', 'row_id', 'first')], 
           id.vars = c('first', 'row_id'))[order(row_id),]
    
    fwrite(fastq_df[which(first), 'value'], paste0(outdir, '/', prefix, '_', errr, 'perc_error_R1.fastq'), col.names = FALSE)
    fwrite(fastq_df[which(!first), 'value'], paste0(outdir, '/', prefix, '_', errr, 'perc_error_R2.fastq'), col.names = FALSE)
  }
}

#' Randomly mutate DNA sequences
#'
#' Introduces random mutations in DNA sequences given a specified error rate.
#'
#' @param seqs DNAStringSet or data.frame with strings.
#' @param error_rate Integer. Percentage of base substituations to introduce.
#' @param seed Integer. Random seed for sampling.
#' @return Returns DNAStringSet or data.frame with mutated sequences.
#' @export
mutateSeqs <- function(seqs,
                       error_rate = 0,
                       seed = 19)
{
  message ('Mutating sequences')
  set.seed(seed)
  stringset <- class(seqs) == "DNAStringSet"
  
  seqs_col <- glue_collapse(seqs)
  seqs_exp <- stringr::str_split_fixed(seqs_col, '', n = nchar(seqs_col))
  
  n_muts <- error_rate / 100 * length(seqs_exp)
  mut_i  <- sample(1:length(seqs_exp), min(length(seqs_exp), round(n_muts + n_muts / 4)))
  
  # replace bases with mutated ones
  message ('Replacing ', length(mut_i), ' positions')
  seqs_exp_mut <- seqs_exp
  seqs_exp_mut[mut_i] <- sample(c('A', 'T', 'G', 'C'), length(mut_i), replace = TRUE)
  
  # collapse vector into string
  seqs_mut <- glue::glue_collapse(seqs_exp_mut, sep = '')
  
  res <- stringr::str_sub(seqs_mut, ((cumsum(width(seqs))) - (width(seqs) -1)), cumsum(width(seqs)))
  
  if (stringset)
  {
    res <- DNAStringSet(res)
  }
  
  return(res)
}    

