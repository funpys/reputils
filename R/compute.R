#' Downsamples data in long data.frame format
#' 
#' 
#' 
#' 
#' @export
downsamp <- function(df, 
                     i = 'id_unique', 
                     j = 'barcode', 
                     n = 'n', 
                     ds = NULL,
                     full = TRUE,
                     seed = 19)
{
  set.seed(seed)
   
  # create data,tabke
  dt <- as.data.table(df)
  
  # select i j n and data.table
  dt_subs <- dt[, c(i, j, n), with = FALSE]
  colnames(dt_subs) <- c('i', 'j', 'n')
  
  # get cell sizes
  cell_size <- dt_subs[, .(n_tot = sum(n)), by = 'j']
  
  # get downsample threshold if ds is NULL
  if(is.null(ds)) 
  {
    ds <- min(median(cell_size$n_tot), max(750, round(quantile(cell_size$n_tot, 0.05))))
  }
  message( 'Downsampling counts to ', ds)
  
  # only retain cells above ds threshold
  dt_filt <- dt_subs[dt_subs$j %fin% cell_size[which(n_tot >= ds), ]$j, ]
  
  # expand dt by number of reads and sample
  dt_filt[, id_row := 1:nrow(dt_filt)]
  dt_exp <- dt_filt[rep(1:nrow(dt_filt), times = n), ]
  
  indeces_samp <- dt_exp[, .(id_row = sample(id_row, replace = FALSE, size = ds)), by = 'j']
  
  # sampled data.frame
  dt_sampled <- dt_filt[indeces_samp[['id_row']], ]
  
  # count ds reads
  dt_samp_n <- dt_sampled[, .(n_ds = .N), by = c('i', 'j')]
  colnames(dt_samp_n) <- c(i, j, paste0(n, '_ds'))
  
  if (full)
  {
    res <- merge(dt, dt_samp_n, all.x = TRUE, by = c(i, j))[which(is.na(n_ds)), n_ds := 0]
  } else {
    res <- merge(dt_samp_n, dt, all.x = TRUE, by = c(i, j))
  }
  
  if (length(unique(res[, .(tot = sum(n_ds, na.rm = TRUE)), by = j]$tot)) > 2) { warning ('Full joining caused problems. Downsampled counts are not equal accross barcodes') }
  
  return(res)
}

#' Computes variance over mean
#' 
#' 
#' 
#' 
#' @export
varmean <- function(counts,
                    min_expr = 10,
                    reg = 10,
                    future_plan = multicore)
{
  message ('Computing variance over mean')
  # aggregate counts per locus and barcode and downsample
  # counts_aggr    <- counts[, .(n = sum(n)), by = c('id_unique', 'barcode')]
  # counts_aggr_ds <- downsamp(counts_aggr, ds = ds, full = FALSE)
  
  if (is.null(counts$n_ds))
  {
    stop ('Please downsample counts first')
  }
  
  # throw lowly expressed features and empty barcodes
  counts_f <- 
    counts[, gene_tot := sum(n), by = 'id_unique'
      ][, cell_ds_tot := sum(n_ds), by = 'barcode'
        ][which(gene_tot >= min_expr & cell_ds_tot > 0), ]
  
  n_cells <- length(unique(counts_f$barcode))
  reg     <- reg / n_cells
  
  # vm_stats using matrix operations
  smat   <- longToSparse(counts_f[, c('id_unique', 'barcode', 'n_ds')])
  f <- function(i)
  {
    mat    <- as.matrix(smat[chunks[[i]], ])
    means  <- rowMeans(mat)
    tot_ds <- rowSums(mat)
    vars   <- matrixStats::rowVars(mat)
    
    vm_stats <-
      data.table(id_unique = rownames(mat),
                 mean_ds = means,
                 var_ds  = vars,
                 tot_ds  = tot_ds,
                 varmean = (reg + vars) / (reg + means))
    return(vm_stats)
  }
  chunks <- chunk(1:nrow(smat), chunk_size = 5e4)

	res       <- plyr::alply(1:length(chunks), 1, .fun = f, .parallel=TRUE)
	vm_stats  <- rbindlist(res)[order(mean_ds), ]
	message("done computing basic gstat, will compute trends")

  # vm_stats using data.table
  # smat_long <- as.data.table(melt(dcast.data.table(counts_f[, c('id_unique', 'barcode', 'n_ds')], id_unique ~ barcode, fun.aggregate = sum, fill = 0, value.var = 'n_ds')))
  # colnames(smat_long) <- c('id_unique', 'barcode', 'n_ds')
  
  # vm_stats <- 
    # smat_long[, .(mean_ds = sum(n_ds) / n_cells, var_ds = sum((n_ds - sum(n_ds) / n_cells)^2) / (n_cells - 1)), by = 'id_unique'
      # ][, varmean := (reg + var_ds) / (reg + mean_ds)
          # ][order(mean_ds), ]

  # calc trend
  vm_stats <- vm_stats[, vm_trend := zoo::rollmedian(log2(varmean), k = 101, fill = 'extend')]
  
  # add type anno if present
  if (!is.null(counts$type))
  {
    vm_stats <-
      merge(vm_stats,
            unique(counts[, c('id_unique', 'type')]),
            all.x = TRUE,
            by = 'id_unique')
  }
  
  # add name anno if present
  if (!is.null(counts$name))
  {
    vm_stats <-
      merge(vm_stats,
            unique(counts[, c('id_unique', 'name')]),
            all.x = TRUE,
            by = 'id_unique')
  }
  
  # add class anno if present
  if (!is.null(counts$class))
  {
    vm_stats <-
      merge(vm_stats,
            unique(counts[, c('id_unique', 'class')]),
            all.x = TRUE,
            by = 'id_unique')
  }
  
  # add counts
  res <- vm_stats
    # merge(vm_stats, 
          # counts[, .(n = sum(n), n_ds = sum(n_ds)), by = 'id_unique'],
          # all.x = TRUE,
          # by = 'id_unique')
  
  return(res)
}