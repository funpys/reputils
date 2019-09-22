#' Compute multiple sequence alignment
#'
#' Computes multiple sequence alignment using mafft.
#'
#' @param msa DNAStringSet. 
#' @return Returns a long data.frame
#' @export
mafft = function(sequences, 
                 ids = NULL,
                 subgroup = FALSE,
                 family = 'te', 
                 save = FALSE, 
                 outdir = NULL,
                 threads = -1, # uses all
                 epsilon = 0.123, 
                 memsave = FALSE,
                 fft = TRUE,
                 genafpair = FALSE,
                 overwrite = FALSE,
                 op = 1.53, # default
                 parttree = FALSE,
                 retree = 2,
                 maxiterate = 0,
                 long = FALSE,
                 local = FALSE)
{
  # if not outdir saves to Repexpr package
  if (is.null(outdir)) { outdir <- system.file('extdata', 'msa', package = 'Repdata') }
  # create flag vector
  memsave   <- ifelse(memsave, '--memsave', '--nomemsave')
  local     <- ifelse(local, '--localpair', '')
  genafpair <- ifelse(genafpair, '--genafpair', '')
  fft       <- ifelse(fft, '--fft', '--nofft')
  parttree  <- ifelse(parttree, '--parttree', '')
  
  flags <- glue("--thread {threads} --ep {epsilon} --op {op} {memsave} {local} {parttree} --retree {retree} {genafpair} --maxiterate {maxiterate} {fft}")
  
  if (!is.null(ids))
  {
    names(sequences) <- ids
  }
  
  input_file <- tempfile()
  writeXStringSet(DNAStringSet(sequences), input_file, format = 'fasta')
  	
  outfile <- paste0(outdir, '/', family, '_msa.fasta')

  if (!file.exists(outfile) | overwrite)
  {
    cmd = glue('mafft {flags} {input_file} > {outfile}')
    message(cmd)
    
    system(cmd)
  } else {
    message('file exists already')
  }
  
  file.remove(input_file)
  
  msa <- importMSA(outfile, long = long)
  
  if (!save)
  {
    file.remove(outfile)
  }
  
  return(msa)
}

#' Deduplicates SAM/BAM files using umi-tools
#' @param bam path to single or multiple bam files for combined deduplication
#' @export
deduplicateBAM <- function(bam,
                           samtools   = FALSE,
                           suffix     = 'dedup',
                           paired     = NULL,
                           outdir     = NULL,
                           align_dist = NULL,
                           ncores     = 1)
{
  if (sum(file.exists(bam)) != length(bam)) { stop ('Some bam files are not found, please provide full path') }
  if (is.null(paired))                      { stop ('Please specify whether SAM/BAM files contain paired-end alignments') }
  if (is.null(align_dist))                  { stop ('Please specify min distance between read alignments for contig splitting') }
  if (length(bam) > 1)                      { stop ('Please provide path to single bam file or to directory containing bams to be merged') }
  
  # check if directory and retrieve files
  if (dir.exists(bam))
  {
    bn   <- basename(bam)
    bam_out <- paste0(bn, '_merged.bam')
    if (!is.null(outdir)) { bam_out <- paste0(outdir, '/', basename(bam_out)) } else { bam_out <- paste0(bam, '/', bam_out) }
    
    # create output files
    bam_contig <- gsub('.bam$|.sam$', '_contig.bam', bam_out)
    bam_dedup  <- gsub('.bam$|.sam$', paste0('_', suffix, '.bam'), bam_contig)
    bai_dedup  <- paste0(bam_dedup, '.bai')
    log_dedup  <- gsub('.bam$', paste0('_', suffix, '.log'), bam_contig)
    
    # get input bams
    bam <- dir(bam, full.names = TRUE, pattern = '.bam$|.sam$')
  } else {

    bam_orig <- bam
    if (!is.null(outdir))
    {
      bam = paste0(outdir, '/', basename(bam))
    }
    bam_contig <- gsub('.bam$|.sam$', '_contig.bam', bam)
    bam_dedup  <- gsub('.bam$|.sam$', paste0('_', suffix, '.bam'), bam_contig)
    bai_dedup  <- paste0(bam_dedup, '.bai')
    log_dedup  <- gsub('.bam$', paste0('_', suffix, '.log'), bam_contig)
    bam        <- bam_orig
  }
  
  if (!dir.exists(outdir))      { stop ('Please provide valid outdir')                 }
  if (file.exists(bam_contig)) { stop ('Cannot overwrite existing file ', bam_contig) }

  # import reads from BAM  
  if (paired)
  {
    what <- c('qname', 'flag', 'mpos', 'isize', 'mrnm')
  } else {
    what <- c('qname', 'flag')
  }

  reads <- 
    importBAM(bam, 
              paired = paired, 
              mate = NA, 
              what = what, 
              tag = c('CB', 'UR', 'NH'), 
              data.table = FALSE)                                          
  
  # change barcode column to CB
  colnames(mcols(reads)) <- gsub('barcode', 'CB', colnames(mcols(reads)))
  
  compContigs <- function(reads,
                          paired = FALSE,
                          align_dist = 1e4)
  {
    # add contig id based on read alignment distance  
    if (paired)
    {
      dt <- data.table(seqnames = as.character(seqnames(reads)),
                       start    = start(reads),
                       end      = end(reads),
                       strand   = as.character(strand(reads)),
                       mpos     = mcols(reads)$mpos,
                       index    = 1:length(reads))
      # adjust start end coordinates according to mate
      dt[which(strand == '+'), end   := pmax(mpos, end)]
      dt[which(strand == '-'), start := pmin(mpos, start)]
      
      # reorder
      dt <- dt[order(seqnames, start, end) ]
    } else {
      dt <- data.table(seqnames = as.character(seqnames(reads)),
                       start    = start(reads),
                       end      = end(reads),
                       index    = 1:length(reads))
    }
    
    # compute delta to previous
    contigs <-
      dt[, delta := start - shift(end, n = 1L, type = 'lag', fill = -1e9), by = 'seqnames'
                   ][, 'CT' := as.character(cumsum(delta >= align_dist))
                      ][order(index), ]
    
    mcols(reads)$CT <- contigs$CT
    return(reads)
  }
  
  # add contig
  reads <- compContigs(reads, paired, align_dist)
  
  # export modified BAM
  exportBAM(reads, file = bam_contig, ncores = ncores)
  
  # execute umi-tools
  if (paired)
  {
    paired_status <- '--paired'
  } else {
    paired_status <- ''
  }
  message ('Deduplicating with umi-tools')
  system(glue("umi_tools dedup --log2stderr --per-cell --cell-tag=CB --per-gene --gene-tag=CT --extract-umi-method=tag --umi-tag=UR {paired_status} --multimapping-detection-method=NH -I {bam_contig} -S {bam_dedup} &> {log_dedup}"))
  
  if (TRUE)
  {
    system(glue("samtools index {bam_dedup} {bai_dedup}"))
  }

  plotContig <- function(contigs)
  {
    contig_stats <- contigs[, .(cont_size = max(starts) - min(starts), cont_nreads = .N), by = 'contig']
  
    p_size   <- contig_stats %>% ggplot(aes(cont_size)) + geom_density() + scale_x_log10()
    p_nreads <- contig_stats %>% ggplot(aes(cont_nreads)) + geom_density() +scale_x_log10()
    
    p_size_vs_nreads <- contig_stats %>% ggplot(aes(cont_size, cont_nreads)) + geom_point()
  }

  return(TRUE)
}

#' Splits BAM file per chromosome
#' @export
splitBAM <- function(bam,
                     samtools  = TRUE,
                     main_only = TRUE,
                     future_plan = multicore)
{
  chroms_orig <- sapply(strsplit(fread(cmd = glue("samtools view -H {bam} | grep ^@SQ"), header = FALSE)$V2, ':', fixed = T), function(x) x[2])
  
  if (sum(grepl('chr', chroms_orig)) == 0)
  {
    chroms_mod <- paste0('chr', chroms_orig)
  }
  
  if (main_only)
  {
    chroms_orig <- chroms_orig[gsub('chr', '', chroms_mod) %in% c(as.character(1:10000), 'MT', 'X', 'Y', 'M')]
    chroms_mod  <- chroms_mod[gsub('chr', '', chroms_mod) %in% c(as.character(1:10000), 'MT', 'X', 'Y', 'M')]
  }
  
  res <- listenv::listenv()
  future::plan(future_plan)
  for (i in (1:length(chroms_orig)))
  {
    bam_chrom <- gsub('.bam$', paste0('_', chroms_mod[i], '.bam'), bam)
    res[[chroms_orig[i]]] %<-% system(glue("samtools view -bh {bam} {chroms_orig[i]} > {bam_chrom}"))
  }
  
  return(as.list(res))
}