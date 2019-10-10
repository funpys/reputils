#' Check format of TE annotation file
#' @export
checkFormat <- function(path)
{
  if (!file.exists(path)) { stop ('Path to file does not exist') }
  
  header <- colnames(fread(path, nrow = 2, fill = TRUE))
  
  rmsk <- paste0('V', 1:14)
  dfam <- c("#seq_name","family_acc","family_name","bits","e-value","bias","hmm-st","hmm-en","strand","ali-st","ali-en","env-st","env-en","sq-len","kimura_div")
  ucsc <- c("#bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart","repEnd","repLeft","id")
  
  type <- NULL
  
  if (!any(!header %in% rmsk)) { type <- 'rmsk'    }
  if (!any(!header %in% dfam)) { type <- 'dfam'    }
  if (!any(!header %in% ucsc)) { type <- 'ucsc'    }
  if (is.null(type))           { type <- 'unknown' }

  return(type)
}

