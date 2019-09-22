#' Transforms long data.frame into sparse matrix
#' 
#' 
#' 
#' 
#' @export
longToSparse <- function(df)
{
	df <- as.data.table(df)
  
  i <- df[[1]]
  j <- df[[2]]
  n <- df[[3]]
	
  if (class(i) != 'factor') { i <- factor(i) }
  if (class(j) != 'factor') { j <- factor(j) }
		
	smat <- 
    Matrix::sparseMatrix(i        = as.integer(i), 
                         j        = as.integer(j), 
                         x        = n,
                         dims     = c(length(levels(i)), length(levels(j))),
                         dimnames = list(levels(i), levels(j)))
	
	return(smat)
}

#' Transforms multiple sequence alignment into long data.frame
#'
#' Transforms multiple sequence alignment into long data.frame
#'
#' @param msa DNAStringSet. 
#' @return Returns a long data.frame
#' @export
msaToLong <- function(msa)
{
  n_alignments <- as.numeric(length(msa))
  n_positions  <- as.numeric(width(msa)[1])
  msa_size     <- n_alignments * n_positions
  
  if (msa_size < 1.5e7)
  {
    mat  <- as.matrix(msa)
    colnames(mat) <- 1:ncol(mat)
    
    res <- 
      data.table(locus = rep(rownames(mat), ncol(mat)), 
                 pos   = as.numeric(rep(colnames(mat), each = nrow(mat))), 
                 nuc = as.vector(mat))[which(nuc != '-'), ]
  } else {  
    chunk_size <- max(1e4, round(n_alignments / 20))
    results    <- list()
    chunks     <- split(1:length(msa), ceiling(seq_along(1:length(msa)) / chunk_size))
    for (i in 1:length(chunks)) 
    {
      chunk <- msa[chunks[[i]]]
      results[[i]] <- future::future({
                                      rbindlist(lapply(chunk, function(locus) {data.table(pos=1:length(locus), nuc=as.vector(locus))[which(nuc!='-'),]}), idcol = 'locus')
                                     })
    }

    res <- rbindlist(lapply(results, future::value))
  }
  
  return(res)
}

.headerToSegInfo <- function(header)
{
  res <- Seqinfo(seqnames = substring(header$V2, 4, 99999), seqlengths = as.integer(substring(header$V3, 4, 99999)))
  return(res)
}

.seqinfoToHeader <- function(seqinfo)
{
  dt <- data.table(seqnames   = seqnames(seqinfo),
                   seqlengths = seqlengths(seqinfo))
                   
  header <- dt[, ':=' (tag = '@SQ', seqnames = paste0('SN:', seqnames), seqlengths = paste0('LN:', seqlengths))][, c('tag', 'seqnames', 'seqlengths')]
  return(header)
}

sparseToLong <- function(smat)
{
  if (is.null(rownames(smat)))
  {
    rownames(smat) <- 1:nrow(smat)
  }
  if (is.null(colnames(smat)))
  {
    colnames(smat) <- 1:ncol(smat)
  }
  
  dt <- as.data.table(Matrix::summary(smat))
  
  dt[, ':=' (i = rownames(smat)[i], j = colnames(smat)[j])]
  colnames(dt) <- c('name', 'barcode', 'n')
  
  return(dt)
}

toMat <- function(df)
{
  rowlabs <- pull(df, 1)
  mat     <- as.matrix(df[,-1])
  rownames(mat) <- rowlabs
  return(mat)
}

