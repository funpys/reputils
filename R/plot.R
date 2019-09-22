# plot split violin plot (copied from internet, acknowledge source!)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

plotEmpty <- function(string = 'NA', size = 3, border = TRUE)
{
  p <- 
    ggplot() + 
      annotate('text', x = 4, y = 25, size = size, label = string) +
      theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
      
  if (!border)
  {
    p <- p + theme(panel.border = element_blank())
  }
  return(p)
}

plotMSA <- function(msa, 
                    cluster = TRUE,
                    type = 'base',
                    max_gap = 0.95,
                    threads = 1,
                    ds = NULL)
{
  if (!is.null(ds))
  {
    msa <- sample(msa, ds)
  }
  
  msa <- DNAMultipleAlignment(msa)
  
  # remove gapped columns
  gap_cols <- which(consensusMatrix(msa, as.prob = T)['-',] > max_gap)
  colmask(msa) <- IRanges(start = gap_cols, end = gap_cols)
  
  msa <- DNAStringSet(msa)
  
  # Remove short loci after gap removal
  msa <- msa[width(DECIPHER::RemoveGaps(msa)) >= 6]
  
  msa <- ape::as.DNAbin(msa)
  
  if (cluster)
  {
    message ('Clustering loci')
    # msa <- DNAStringSet(msa)
    # d <- DECIPHER::DistanceMatrix(msa, processors = 70)
    # clusts <- IdClusters(myDistMatrix = d, method = "UPGMA", processors = 70)
    # msa <- msa[order(clusts[,1])]
    
    #d <- tgs_dist(oligonucleotideFrequency(msa, width = 6))
    
    #dist.dna(msa[1:1])
    #d <- DECIPHER::DistanceMatrix(msa, type = 'dist')
    d <- dist.dna(as.DNAbin(msa), model = 'indelblock')
    if (length(msa) > 100)
    {
      #clusts <- hclust(d, 'ward.D2')
      
      #test = DECIPHER::IdClusters(myDistMatrix = d, method = "UPGMA")
      hc <- hclust(d, 'ward.D2')
      msa <- msa[hc$order]
    } else {
      clusts <- ape::fastme.bal(d)
      #msa <- msa[names(msa)[as.numeric(clusts$tip.label)]]
      msa <- msa[clusts$tip.label]
    }
  }
    
  if (type == 'base')
  {
    #msa_bin <- as.DNAbin(msa)
    
    print(ape::image.DNAbin(ape::as.DNAbin(msa), col = c("#5DA731", "#D63317", "#DFD93C", "#2F4C9B", "lightgrey", "white"), show.labels = TRUE, cex.lab = 0.5))
  } else {
    msa_mat <- as.matrix(msa)
    msa_dt <- as.data.table(msa_mat)
    msa_dt$id_unique <- rownames(msa)
    msa_long <- melt(msa_dt, id.vars = 'id_unique')
    
    # plot ggplot
    msa_long %>% 
      ggplot(aes(variable, id_unique, fill = value)) + 
        geom_raster() +
        scale_fill_manual(values = c("white", "#5DA731", "#D63317", "#DFD93C", "#2F4C9B")) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              panel.background = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
  }
  
  return(clusts$tip.label)
}