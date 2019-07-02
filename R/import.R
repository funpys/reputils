#' Parallel import of BAM files
#'
#' Imports reads stored in BAM files
#'
#' @param path Character. Path to file.
#' @param main_assembly Logical. If TRUE (the default), only returns TE intervals on the main chromosomes. 
#' @return Returns a GRanges object with TE intervals
#' @export
importBAM <- function(bams, 
                      bai_pos = 'suffix', 
                      paired  = NULL, 
                      multi   = TRUE,
                      what    = NULL,
                      mate    = NA,
                      anchor  = NULL,
                      filt    = TRUE,
                      barcode = 'CB')
{
  if (class(bams) == 'data.frame')
  {
    bam_files <- as.character(bams$paths)
    if (!is.null(bams$paired))  { paired  <- as.character(bams$paired) }
    if (!is.null(bams$mate))    { mate    <- as.character(bams$mate) }
    if (!is.null(bams$barcode)) { barcode <- as.character(bams$barcode) }
  } else {
    bam_files <- bams
    paired    <- rep(paired,  length(bams))
    mate      <- rep(mate,    length(bams))
    barcode   <- rep(barcode, length(bams))
  }
  
  if (sum(dirname(bam_files) == '.') > 0) { stop ('Please provide full path to files') }
  
  # check paired status
  if (is.null(paired))
  { 
    # get matching indicies
    if (bai_pos == 'suffix')
    {
      bam_indices <- gsub('.bai$', '', bam_files)
    } else {
      bam_indices <- gsub('.bam$', '', bam_files) 
    }
    
    paired <- Rsamtools::testPairedEndBam(bam_files)
  }
  
  # actual reading function
  readBamF <- function(bam_file,
                       #bam_index,
                       paired,
                       barcode, 
                       mate,
                       what)
  {
    message(glue("Reading {bam_file}"))
    
    # create param
    param <- ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired         = as.logical(paired), 
                                                       isProperPair     = as.logical(paired), 
                                                       isFirstMateRead  = mate == 'first',
                                                       isSecondMateRead = mate == 'second'))
                                                     
    if (barcode == 'CB')
    {
      param@tag <- c('NH', 'CB')
    } else {
      param@tag <- c('NH')
    }
    
    if (!is.null(what))
    {
      param@what <- what
    } 
    
    res <- 
      GenomicAlignments::readGAlignments(bam_file, 
                                         #index = bam_index,
                                         use.names = FALSE,
                                         param = param)
                                         
    if (class(res) == 'list') { res <- unlist(res) }
    res <- granges(res, use.mcols = TRUE, use.names = TRUE)      
    
    if (barcode != 'CB')
    {
      res$barcode <- barcode
    }
    
    # add mate status
    if (!is.null(res$flag))
    {
      res$first <- as.logical(bamFlagAsBitMatrix(res$flag)[, 'isFirstMateRead'])
    }
    
    return(res)
  }
  
  # import in parallel
  message ('Importing BAM file(s)')
  results <- list()
  for (i in 1:length(bam_files))
  {
    results[[i]] <- 
      future({ readBamF(bam_files[i], 
                        #bam_indices[i], 
                        barcode = barcode[i], 
                        paired  = paired[i], 
                        mate    = mate[i],
                        what    = what) 
            })
  }
  reads <- unlist(as(lapply(results, value), 'GRangesList'))
  message ('Done')

  # add NH weighted column and rename CB to barcode
  vect <- c('unique', rep('multi', 50))
  reads$NH_weight <- 1 / reads$NH
  reads$NH_flag   <- vect[reads$NH]
             
  if (!is.null(reads$CB))
  {
    reads$barcode <- reads$CB
    reads <- reads %>% select(-CB)
  }
  
  if (!multi)
  {
    reads <- reads[reads$NH == 1]
  }
  
  if (filt)
  {
    reads <-
      reads %>%
      filter(!is.na(barcode)) #%>%   
      #filter(!is.na(UB))
  }
  
  # make UCSC seqlevelsStyle and filter chroms
  GenomeInfoDb::seqlevelsStyle(reads) <- 'UCSC'       
  reads <- GenomeInfoDb::keepStandardChromosomes(reads, pruning.mode = 'coarse')
  
  # anchor coordinates
  if (!is.null(anchor))
  {
    if (anchor == 'fiveprime')
    {
      reads <- mutate(anchor_5p(reads), width = 1)
    }
    if (anchor == 'threeprime')
    {
      reads <- mutate(anchor_3p(reads), width = 1)
    }
  }
  
  return(reads)
}

#' Import TE intervals from DFAM
#'
#' Imports TE interval coordinates from DFAM database
#'
#' @param path Character. Path to file.
#' @param curate Logical. If TRUE, calls `curateTEs` to resolve overlapping TE intervals.
#' @param main_assembly Logical. If TRUE (the default), only returns TE intervals on the main chromosomes. 
#' @return Returns a GRanges object with TE intervals
#' @export
importDFAM <- function(path,
                       curate = FALSE,
                       main_assembly = TRUE)          
{
  dfam_hits <- fread(path)
  colnames(dfam_hits) <- gsub('-', '_', colnames(dfam_hits))
  colnames(dfam_hits) <- gsub('#', '', colnames(dfam_hits))
  
  res <- 
    dfam_hits %>% 
      mutate(start = ifelse(strand == '+', ali_st, ali_en),
             end = ifelse(strand == '+', ali_en, ali_st),
             seqnames = seq_name, 
             name = family_name,
             id_unique = paste(name, 1:nrow(.), sep = '|'))  
      as_granges()
      
  if (main_assembly)
  {
    res <- GenomeInfoDb::keepStandardChromosomes(res, pruning.mode = 'coarse')
  }
  
  if (curate)
  {
    res <- curateTEs(res)
  }
    
  return(res)
}

importMSA <- function(fnames,
                      long = TRUE,
                      n_jobs = 500)
{
  # f_sizes <- 
    # tibble(fname = fnames, size = file.info(fnames)$size) %>%
    # arrange(size) %>%
    # mutate(chunk = cut(cumsum(size), breaks = seq(1, max(cumsum(size)) + 1e9, by = 1e9), labels = FALSE, include.lowest = TRUE))
  # test <- future.apply::future_by(f_sizes[1:5,], f_sizes$chunk[1:5], function(x) readDNAStringSet(x$fname))
  
  if (length(fnames) > 1)
  {
    # create commands for distribution
    cmds <- list()
    for (fname in fnames)
    {
      if (long)
      {
        cmds[[basename(fname)]] <- glue("Repexpr:::msaToLong(readDNAStringSet('{fname}'))")
      } else {
        cmds[[basename(fname)]] <- glue("readDNAStringSet('{fname}')")
      }
    }
  
    # create new empty env and fill with relevant
    empty <- new.env(parent = emptyenv())
    
    # distribute, compute and combine
    res_gclust <- 
      gcluster.run3(command_list = cmds,  
                    packages = c("Biostrings", "base", "stats"), 
                    job_names = names(cmds),
                    max.jobs = n_jobs, 
                    envir = empty, 
                    io_saturation = FALSE)
    
    res <- rbindlist(lapply(res_gclust, function(x) x$retv), use.names = FALSE, fill = FALSE, idcol = NULL)
  } else {
    if (long)
    {
      res <- msaToLong(readDNAStringSet(fnames))
    } else {
      res <- readDNAStringSet(fnames)
    }
  }
  
  return(res)
}

#' Import TE intervals from Repeatmasker
#'
#' Imports TE interval coordinates from Repeatmasker
#'
#' @param path Character. Path to file.
#' @param curate Logical. If TRUE, calls `curateTEs` to resolve overlapping TE intervals.
#' @param main_assembly Logical. If TRUE (the default), only returns TE intervals on the main chromosomes. 
#' @param main_classes Logical. If TRUE (the default), only returns TE intervals unambiguously belonging to DNA, LINE, SINE, and LTR classes.
#' @param proper_alignments Logical. If TRUE (the default), only returns TE intervals with proper start and end coordinates.
#' @return Returns a GRanges object with TE intervals
#' @export  
importRMSK <- function(path, 
                       curate = FALSE,
                       main_assembly = TRUE,
                       main_classes = TRUE,
                       proper_alignments = TRUE)          
{
  rmsk_hits <- 
    fread(path, skip = 3, 
          fill = TRUE,
          header = FALSE,
          data.table = TRUE,
          col.names = c('sw_score', 'perc_div', 'perc_del', 'perc_ins', 'seqnames', 'start', 'end', 'left_chrom', 'strand', 'name', 'repclass_family', 
								        'position_start', 'position_end', 'left_rep', 'id_unique'))
  
  rmsk_hits[, c("repclass", "repfamily") := tstrsplit(repclass_family, "/", fixed=TRUE)]
  rmsk_hits[, strand := gsub('C', '-', strand, fixed = TRUE)]
  
  # keep only major repclasses
  if (main_classes)
  {
    rmsk_hits <- rmsk_hits[which(repclass %in% c('DNA', 'SINE', 'LTR', 'LINE')), ]
  }

  # adjust position on alignment based on strand
  rmsk_hits <- rmsk_hits[which(strand == '-'), position_start := (left_rep)][, !'left_rep']
  rmsk_hits$position_start <- as.integer(rmsk_hits$position_start)
  
  if (proper_alignments)
  {
    rmsk_hits <- rmsk_hits[which(position_start < position_end), ]
  }
   
  res <-
    rmsk_hits %>%
    mutate(id_unique = paste(name, 1:nrow(.), sep = '|')) %>%
    select(-left_chrom, -repclass_family) %>%    
    as_granges()
    
  if (main_assembly)
  {
    res <- GenomeInfoDb::keepStandardChromosomes(res, pruning.mode = 'coarse')
  }
  
  if (curate)
  {
    res <- curateTEs(res)
  }
    
  return(res)
}

#' Import TE intervals from UCSC table browser
#'
#' Import TE intervals from UCSC table browser
#'
#' @param path Character. Path to file.
#' @param curate Logical. If TRUE, calls `curateTEs` to resolve overlapping TE intervals.
#' @param main_assembly Logical. If TRUE (the default), only returns TE intervals on the main chromosomes. 
#' @param main_classes Logical. If TRUE (the default), only returns TE intervals unambiguously belonging to DNA, LINE, SINE, and LTR classes.
#' @return Returns a GRanges object with TE intervals
#' @export  
importUCSC <- function(path = 'path to UCSC RMSK file',
                       curate = FALSE,
                       main_classes = TRUE,
                       main_assembly = TRUE)
{
  # Import
  ucsc_rmsk <-
    fread(path,
          data.table = TRUE,
          col.names = c('bin', 'sw_score', 'milli_div', 'milli_del', 'milli_ins', 'seqnames', 'start', 'end', 'geno_left', 'strand', 'name', 
								        'repclass', 'repfamily', 'position_start', 'position_end', 'left_rep', 'id_unique'))
  
  # adjust position on alignment based on strand
  ucsc_rmsk                <- ucsc_rmsk[which(strand == '-'), position_start := (left_rep)][, !'left_rep']
  ucsc_rmsk$position_start <- as.integer(ucsc_rmsk$position_start)
  
  # keep only major repclasses
  if (main_classes)
  {
    ucsc_rmsk <-  ucsc_rmsk[which(repclass %in% c('DNA', 'SINE', 'LTR', 'LINE')), ]
  }
  
  res <-
    ucsc_rmsk%>%
    mutate(id_unique = paste(name, 1:nrow(.), sep = '|')) %>%
    as_granges()
    
  if (main_assembly)
  {
    res <- GenomeInfoDb::keepStandardChromosomes(res, pruning.mode = 'coarse')
  }

  if (curate)
  {
    res <- curateTEs(res)
  }
  return(res)
}