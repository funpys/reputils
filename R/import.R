#' Parallel import of BAM files
#'
#' Imports reads stored in BAM files
#'
#' @param path Character. Path to file.
#' @param main_assembly Logical. If TRUE (the default), only returns TE intervals on the main chromosomes. 
#' @return Returns a GRanges object with TE intervals
#' @export
importBAM <- function(bams, 
                      paired       = NA, 
                      multi        = TRUE,
                      what         = NA,
                      tag          = NA,
                      mate         = NA,
                      anchor       = NULL,
                      meta         = NA,
                      granges      = FALSE,
                      barcode      = NA,
                      data.table   = TRUE,
                      sorted       = TRUE,
                      use_gcluster = FALSE)
{
  # construct data.frame
  if (sum(class(bams) == 'data.frame') > 0)
  {
    if (is.null(bams$barcode))  { bams$barcode <- barcode }
    if (is.null(bams$chunk))    { bams$chunk   <- 1       }
  } else {
    bams <-
      data.frame(paths   = bams,
                 paired  = rep(paired,  length(bams)),
                 mate    = rep(mate,    length(bams)),
                 barcode = rep(barcode, length(bams)),
                 meta    = rep(meta,    length(bams)),
                 chunk   = 1,
                 stringsAsFactors = FALSE)
  }
  
  # character to factor columns
  bams$paths   <- as.character(bams$paths)
  bams$barcode <- as.character(bams$barcode)
  
  # adjust read in fields according to barcode
  if (sum(nchar(na.omit(bams$barcode))) > 0)
  {
    tag <- unique(as.character(na.omit(c(tag, bams$barcode[nchar(bams$barcode) == 2]))))
  }
  
  # pre-checking
  if (is.na(paired) & is.na(bams$paired)) { stop ('Please specify paired status of bam files')  }
  if (is.na(mate)   & is.na(bams$mate))   { stop ('Please specify mate to read from bam files') }
  
  if (sum(dirname(bams$paths) == '.') > 0) { stop ('Please provide full path to files') }
  
  # import in parallel (nested)
  message ('Importing BAM file(s)')
  if (use_gcluster)
  {
     # commands for mafft alignment and conversion
    cmds <- list()
    for (i in unique(bams$chunk))
    {
      cmds[[i]] <- 
        glue("Reputils:::.readBamF(bams[bams$chunk == {i}, ], what = {what}, tag = {tag}, data.table = {data.table})")
    }
    # create new empty env and fill with relevant
    empty <- new.env(parent = emptyenv())
    empty$bams <- bams
    
    # distribute, compute and combine
    res <- 
      gcluster.run3(command_list = cmds,  
                    packages = c("data.table", "base", "stats"), 
                    job_names = names(cmds),
                    max.jobs = 50, 
                    envir = empty, 
                    io_saturation = FALSE)
    if (data.table)
    {
      reads <- rbindlist(lapply(res, function(x) x$retv))
    } else {
      reads <- unlist(GAlignmentsList(lapply(res, function(x) x$retv)))
    }
  } else {
    res <- listenv::listenv()
      for (i in unique(bams$chunk))
      {
        print(i)
        res[[i]] %<-% Reputils:::.readBamF(bams[bams$chunk == i, ], what = what, tag = tag, data.table = data.table) %packages% "data.table"
      }
    if (data.table)
    {
      reads <- rbindlist(as.list(res))
    } else {
      reads <- unlist(GAlignmentsList(as.list(res)))
    }
  }
  message ('Done')

  # throw multi mapping reads if desired
  if (!multi)
  {
    if (data.table)
    {
      reads <- reads[which(NH == 1), ]
    } else {
      reads <- reads[mcols(reads)$NH == 1]
    }
  }
  
  # anchor coordinates
  if (!is.null(anchor))
  {
    message ('Anchoring reads at ', anchor)
    if (anchor == 'fiveprime')
    {
      #reads <- mutate(anchor_5p(reads), width = 1)
      reads[which(strand == '+'), end   := start]
      reads[which(strand == '-'), start := end]
    }
    if (anchor == 'threeprime')
    {
      #reads <- mutate(anchor_3p(reads), width = 1)
      reads[which(strand == '+'), start := end]
      reads[which(strand == '-'), end   := start]
    }
  }
  
  if (!data.table)
  {
    if (GenomeInfoDb::seqlevelsStyle(reads) != 'UCSC')
    {
      GenomeInfoDb::seqlevelsStyle(reads) <- 'UCSC' 
      if (!is.null(mcols(reads)$mrnm))
      {
        mcols(reads)$mrnm <- paste0('chr', mcols(reads)$mrnm)
        mcols(reads)$mrnm <- gsub('^chrMT', 'chrM', mcols(reads)$mrnm)
      }
    }
  } else {
    if (seqlevelsStyle(as_granges(reads[sample(1:nrow(reads), min(nrow(reads), 1000)),])) != 'UCSC')
    {
      reads <- reads[, seqnames := paste0('chr', seqnames)]
      reads <- reads[which(seqnames == 'chrMT'), seqnames := 'chrM']
      if (!is.null(reads$mrnm))
      {
        reads <- reads[, mrnm := paste0('chr', mrnm)]
        reads <- reads[which(mrnm == 'chrMT'), mrnm := 'chrM']
      }
    }  
  }
  
  if (sorted)
  {
    if (data.table)
    {
      reads <- reads[order(seqnames, start, end), ]
    } else {
      reads <- sort(reads, ignore.strand = TRUE)
    }  
  }
  
  return(reads)
}

#' Import TE intervals in DFAM format
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
      Reputils:::gcluster.run3(command_list = cmds,  
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

#' Import TE intervals in Repeatmasker format
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
          col.names = c('sw_score', 'perc_div', 'perc_del', 'perc_ins', 'seqnames', 'start', 'end', 'left_chrom', 'strand', 'name', 'class_family', 
								        'position_start', 'position_end', 'left_rep', 'id_unique'))
  
  rmsk_hits[, c("class", "repfamily") := tstrsplit(class_family, "/", fixed=TRUE)]
  rmsk_hits[, strand := gsub('C', '-', strand, fixed = TRUE)]
  
  # keep only major classes
  if (main_classes)
  {
    rmsk_hits <- rmsk_hits[which(class %in% c('DNA', 'SINE', 'LTR', 'LINE')), ]
  }

  # adjust position on alignment based on strand
  rmsk_hits <- rmsk_hits[which(strand == '-'), ':=' (position_start = as.character(position_end), position_end = as.integer(left_rep))][, !'left_rep']
  rmsk_hits$position_start <- as.integer(rmsk_hits$position_start)
  
  if (proper_alignments)
  {
    rmsk_hits <- rmsk_hits[which(strand == '+' & position_start < position_end | strand == '-' & position_start > position_end), ]
  }
   
  res <-
    rmsk_hits %>%
    mutate(id_unique = paste(name, 1:nrow(.), sep = '|')) %>%
    select(-left_chrom, -class_family) %>%    
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

#' Import TE intervals in UCSC table browser format
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
								        'class', 'repfamily', 'position_start', 'position_end', 'left_rep', 'id_unique'))
  
  # adjust position on alignment based on strand
  ucsc_rmsk                <- ucsc_rmsk[which(strand == '-'), position_start := (left_rep)][, !'left_rep']
  ucsc_rmsk$position_start <- as.integer(ucsc_rmsk$position_start)
  
  # keep only major classes
  if (main_classes)
  {
    ucsc_rmsk <-  ucsc_rmsk[which(class %in% c('DNA', 'SINE', 'LTR', 'LINE')), ]
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

import10x <- function(path,
                      long = TRUE)
{
  if (!requireNamespace('hdf5r', quietly = TRUE)) 
  {
    stop("hdf5r package required to read file, please install!")
  }
  if (!file.exists(path)) 
  {
    stop("File not found")
  }

  h5f      <- hdf5r::H5File$new(path, mode = 'r+')
  assembly <- names(x = h5f)
  
  if (length(assembly) > 1)
  {
    stop('More than 1 genome in file')
  }
  
  # adjust depending on cell ranger version
  if (h5f$attr_exists("PYTABLES_FORMAT_VERSION")) {
    gene_names <- 'gene_names'
    ensembl    <- 'genes'
  } else {
    gene_names <- 'features/name'
    ensembl    <- 'features/id'
  }
  
  counts   <- h5f[[paste0(assembly, '/data')]]
  indices  <- h5f[[paste0(assembly, '/indices')]]
  indptr   <- h5f[[paste0(assembly, '/indptr')]]
  dims     <- h5f[[paste0(assembly, '/shape')]]
  genes    <- h5f[[paste0(assembly, '/', gene_names)]][]
  ensembl  <- h5f[[paste0(assembly, '/', ensembl)]][]
  barcodes <- h5f[[paste0(assembly, '/barcodes')]][]
  smat <- 
    sparseMatrix(
                 i = indices[] + 1,
                 p = indptr[],
                 x = as.numeric(x = counts[]),
                 dims = dims[],
                 giveCsparse = TRUE
                )             
  
  # resolve ambigious gene names
  gene_ids_dt <- 
    data.table(genes = genes, 
               ensembl = ensembl, 
               duplicated = duplicated(genes))[which(duplicated), genes := paste(genes, ensembl, sep = '|')]

  rownames(smat) <- gene_ids_dt$genes
  colnames(smat) <- barcodes
  
  if (long)
  {
    res <- sparseToLong(smat)
  } else {
    res <- smat
  }
  
  h5f$close_all()
  return(res)
}
  
# actual reading function
  .readBamF <- function(df, what = NA, tag = NA, data.table = TRUE)
  {
    res <-
      future.apply::future_apply(df, 1, future.packages = 'data.table', function(x)
      {
        path    <- x['paths']
        paired  <- x['paired']
        barcode <- x['barcode']
        mate    <- x['mate']
        meta    <- x['meta']
        
        mate_1 <- mate == 'first'
        mate_2 <- mate == 'second' | mate == 'last'
        
        # check file type
        if (sum(substring(path, nchar(path) - 2, nchar(path)) %in% c('bam', 'BAM', 'Bam')) > 0) { file_type <- 'bam' }
        if (sum(substring(path, nchar(path) - 2, nchar(path)) %in% c('sam', 'SAM', 'Sam')) > 0) { file_type <- 'sam' }
        
        message(glue("Reading {path}"))

        # create param
        param <- ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired          = as.logical(paired), 
                                                            isProperPair      = as.logical(paired), 
                                                            isFirstMateRead   = mate_1,
                                                            isSecondMateRead  = mate_2))

        if (!is.na(tag))
        {
          param@tag <- tag
        }
        
        if (!is.na(what))
        {
          param@what <- what
        } 
        
        # import reads from bam
        if (file_type == 'bam')
        {
          res <- 
            GenomicAlignments::readGAlignments(path, 
                                               #index = bam_index,
                                               use.names = FALSE,
                                               param = param)
        }
                                             
        # import reads from sam (MARS-seq specific)
        if (file_type == 'sam')
        {
          header <- .headerToSegInfo(fread(cmd       = glue('head -n 1000 {path} | grep "^@SQ"'), header = FALSE))
          
          reads <-
            fread(cmd       = glue('grep -v "^@" {path}'), # skip header
                  fill      = TRUE,
                  select    = c(1, 2, 3, 4, 6, 10),
                  col.names = c('qname', 'flag', 'seqnames', 'start', 'cigar', 'NH'))
                  
          # filter unmapped and barcode/UMI missing reads
          reads_f <- 
            reads[which(seqnames != '*'),     # throw unmapped reads
              ][, qname_nchar := nchar(qname)
                ][which(qname_nchar > 100), ]  # throw reads without barcode/umi suffix
              
          # add strand and barcode/UMI from read name
          strand              <- bamFlagAsBitMatrix(reads_f$flag)[, 'isMinusStrand']
          strand[strand == 1] <- '-'
          strand[strand == 0] <- '+'
          reads_f <- 
            reads_f[, ':=' (end = start, 
                            strand = strand, 
                            barcode = substring(qname, qname_nchar - 20, qname_nchar - 9), 
                            UR = substring(qname, qname_nchar - 7, qname_nchar), 
                            NH = as.integer(substring(NH, 6, 10)))][, !'qname_nchar']  
          
          # filter incomplete UMIs/barcodes
          res <- reads_f[-which(grepl('N', UR) | grepl('N', barcode)),]
          
          # convert to GAlignments
          if (!data.table)
          {
            res            <- as(as_granges(res), 'GAlignments')
            seqlevels(res) <- seqlevels(header)
            seqinfo(res)   <- header
          }
        }
        
        # transform to data.table and only keep relevant columns
        if (data.table)
        {
          res <- as.data.table(res)
          res <- res[, colnames(res) %in% c('seqnames', 'strand', 'start', 'end', 'NH', 'UR', 'barcode', barcode, what, tag), with = FALSE]
        }
        
        # modify barcode column names
        if (!is.na(barcode) & file_type != 'sam') 
        {
          if (barcode == 'CB')
          {
            if (data.table)
            {
              res <- res[, barcode := CB][, !'CB']
            } else {
              colnames(mcols(res)) <- gsub('CB', 'barcode', colnames(mcols(res)))
            }
          } else {
            if (data.table)
            {
              res[, 'barcode' := barcode]
            } else {
              mcols(res)$barcode = barcode
            }
          }
        }

        if (!is.na(meta))
        {
          if (data.table)
          {
            res[, ':=' ('meta' = meta, 'barcode' = paste(barcode, meta, sep = '|'))]
          } else {
            mcols(res)$barcode <- paste(mcols(res)$barcode, meta, sep = '|')
          }
        }

        return(res)
      })
    
    if (data.table)
    {
      res <- rbindlist(res)
    } else {
      res <- unlist(GAlignmentsList(res))
    }
    return(res)
  }