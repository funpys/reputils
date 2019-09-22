#' Export (sparse) matrix in matrix market format
#' Much faster than writeMM from Matrix R package
#' @export
writeMM2 <- function(smat, file = NULL)
{
  to_export <- Matrix::summary(smat)
  type <- class(to_export$x)
  type <- ifelse(type == 'numeric', 'real', 'integer')
  colnames(to_export) <- c(nrow(smat), ncol(smat), nrow(to_export))
  
  cat(paste0('%%MatrixMarket matrix coordinate ', type, ' general\n'), file = file)
  fwrite(to_export, file, col.names = TRUE, row.names = FALSE, sep = ' ', append = TRUE)
}

#' @export
exportBAM <- function(GA,
                      file     = NULL,
                      samtools = TRUE,
                      ncores   = 1)
{
  if (is.null(file))              { stop ('Please provide path to output file') }
  if (file.exists(file))          { stop ('Cannot overwrite existing file') }
  if (class(GA) != "GAlignments") { stop ('Input must be GAlignments class') }
  
  file_sam <- gsub('.bam$', '.sam', file)
  file_bam <- file
  file_bai <- gsub('.bam$', '.bai', file)
  
  if (samtools)
  {
    message ('Exporting BAM file using samtools')
    header <- .seqinfoToHeader(seqinfo(GA))
    dt     <- as.data.table(GA)
    
    # add missing but required columns
    if (is.null(dt$mapq))
    {
      dt[, mapq := 0]
    }
    if (is.null(dt$mrnm))
    {
      dt[, mrnm := '*']
    } else {
      dt[which(seqnames == mrnm), mrnm := '=']
    }
    if (is.null(dt$mpos))
    {
      dt[, mpos := 0]
    }
    if (is.null(dt$tlen))
    {
      dt[, tlen := 0]
    }
    if (sum(colnames(dt) == 'seq') == 0) # don't change to dt$seq (finds seqnames)
    {
      dt[, ':=' (seq = '*', qual = '*')]
    }
        
    # create column universe in right order
    tags <- colnames(dt)[nchar(colnames(dt)) == 2]
    cols <- c('qname', 'flag', 'seqnames', 'start', 'mapq', 'cigar', 'mrnm', 'mpos', 'tlen', 'seq', 'qual', tags)
    
    # modify tag columns to indicate class
    types       <- c('integer' = 'i', 'character' = 'Z')
    classes     <- sapply(dt, class)[tags]
    tag_types   <- types[classes]
    for (i in 1:length(tags))
    {
      prefix <- paste(tags[i], tag_types[i], sep = ':')
      dt[, tags[i] := paste(prefix,  dt[, get(tags[i])], sep = ':')]
    }
  
    # arrange columns to match sam format
    dt_sam <- dt[, ..cols, with = TRUE]
    
    # export
    fwrite(header, file_sam, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)                # write header
    fwrite(dt_sam, file_sam, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE, append = TRUE) # append alignments
    
    # convert to BAM
    message ('Converting SAM to BAM using ', ncores, ' cores (modify by setting ncores parameter)')
    system(glue("samtools view -bS -@{ncores} {file_sam} > {file_bam}")) # sam to bam
    system(glue("samtools index {file_bam} {file_bai}"))                 # index bam
    
    if (TRUE)
    {
      file.remove(file_sam)
    }
  } else {
    rtracklayer::export(GA, BamFile(file))
  }
}