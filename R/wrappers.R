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
                 overwrite = FALSE,
                 op = 1.53, # default
                 parttree = FALSE,
                 retree = 2,
                 maxiterate = 0,
                 long = FALSE,
                 local = FALSE)
{
  # if not outdir saves to Repexpr package
  if (is.null(outdir)) { outdir <- system.file('extdata', 'msa', package = 'Repexpr') }
  
  # create flag vector
  memsave  <- ifelse(memsave, '--memsave', '--nomemsave')
  local    <- ifelse(local, '--localpair', '')
  fft      <- ifelse(fft, '--fft', '--nofft')
  parttree <- ifelse(parttree, '--parttree', '')
  
  flags <- glue("--thread {threads} --ep {epsilon} --op {op} {memsave} {local} {parttree} --maxiterate {maxiterate} {fft}")
  
  names(sequences) <- ids
  
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

