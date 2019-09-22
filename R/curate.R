#' Curation of gene intervals
#' @export
curateGenes <- function(genes, 
                        extend = 0)
{
  message ('Curating gene intervals')
 
  # delete bio from column names for consistency between 10x and gencode gtf files
  colnames(mcols(genes)) <- gsub('bio', '', colnames(mcols(genes)))
    
  # modify duplicated gene names
  genes_dup <- 
    genes %>% 
    as_tibble %>% 
    select(gene_id, gene_name) %>% 
    distinct %>% 
    count(gene_name) %>% 
    filter(n > 1) %>% 
    pull(gene_name)
  genes <- 
    genes %>%
      mutate(gene_name = ifelse(gene_name %in% genes_dup, paste(gene_name, gene_id, sep = '_'), gene_name))  
  
  exons <- genes[genes$type == 'exon']
  
  # extend intervals at 5'
  exons <- stretch(anchor_3p(exons), extend)
  
  # reduce and disjoin exon intervals per gene
  exons_uniq <- 
    exons %>% 
    group_by(gene_name) %>% 
    reduce_ranges_directed %>% 
    disjoin_ranges_directed %>%
    mutate(id_unique = 1:length(.))
    
  # find annotation for new intervals
  anno <-
    join_overlap_inner_directed(exons_uniq, exons) %>%
      group_by(id_unique) %>%
      summarise(name  = paste(unique(gene_name), collapse = '|'),
                gene_id    = paste(unique(gene_id), collapse = '|'),
                class      = paste(unique(gene_type), collapse = '|'))
              
  # combine intervals and anno
  exons_uniq_anno <-
    merge(as.data.table(exons_uniq), 
          as.data.table(anno),
          all=FALSE,
          by = 'id_unique')[, !'id_unique'
            ][, id_unique := paste(name, 1:.N, sep = '|'), by = 'name']
    
  # add position start/end column (similar to TEs) to map reads on model (corresponds to distance of first exon pos to most upstream TSS)
  exons_uniq_anno <- 
    exons_uniq_anno[which(strand == '+'), ':=' (position_start = as.integer(cumsum(width) - width + 1), position_end = as.integer(cumsum(width))),  by = 'name'
      ][which(strand == '-'), ':=' (position_start = as.integer(rev(cumsum(rev(width))))), by = 'name'][which(strand == '-'), position_end := as.integer(position_start - width + 1)]
  
  # add external annotation(TSS, UTR, etc.)
  exons_uniq_anno[, ':=' (tss = NA, poly = NA, tes = NA)]
  # tss_nested      <- as.data.table(join_overlap_left_directed(getTSS(genes), as_granges(exons_uniq_anno)) %>% select(id_unique))[, list(tss = list(.SD)), by = 'id_unique']
  # exons_uniq_anno <- merge(exons_uniq_anno, tss_nested, all.x = TRUE, by = 'id_unique')
 
  res <- as_granges(exons_uniq_anno)
 
  return(res)
}

#' Curation of TE intervals
#' @export
curateTEs <- function(tes, 
                      genes = NULL, # should be curated
                      extend_3p = 0,
                      filt_exonic = FALSE,
                      max_e_value = 0.05, 
                      min_size = 25)
{
  message ('Curating TE intervals')
  if (class(tes) != 'GRanges') { tes <- as_granges(tes) }
  
  if (!is.null(tes$e_value))
  {
    # filter based on e_value threshold
    tes <- 
      tes %>% 
      filter(e_value < max_e_value) %>%
      select(-bits, -e_value, -bias, -sq_len, -kimura_div, -family_acc)
  }
  
  message ('Resolving overlaps')
  # throw overlapping/nested repeats
  tes_unnest <-
    tes %>% 
    disjoin_ranges_directed %>% 
    mutate(n_overl = countOverlaps(., tes)) %>% 
    filter(n_overl == 1) %>% 
    join_overlap_left_directed(., tes) %>%
    select(-n_overl)
    
  # filter short intervals
  tes_unnest_filt <-
    tes_unnest %>% 
    filter(width >= min_size)
    
  message ('Labeling modified intervals')
  # add column to indicate if interval was modified or not
  tes_unnest_filt <-
    tes_unnest_filt %>% 
    mutate(orig = 1:length(.) %in% findOverlaps(tes_unnest_filt, tes, type = 'equal')@from,
           position_start = ifelse(orig, position_start, NA),
           position_end   = ifelse(orig, position_end, NA))
    
  # add new unique id (alphabetical prefix for locus chunks)
  tes_unnest_filt_df <- as.data.table(tes_unnest_filt)
  
  id_conv_df <-
    tes_unnest_filt_df[, fragmented := duplicated(id_unique)
     ][which(id_unique %in% id_unique[fragmented]), 
      ][, .(id_unique2 = paste0(letters[1:.N], id_unique)), by = 'id_unique']
      
  tes_unnest_filt[tes_unnest_filt$id_unique %in% id_conv_df$id_unique]$id_unique <- id_conv_df$id_unique2
  
  if (length(unique(tes_unnest_filt$id_unique)) != length(tes_unnest_filt)) { stop ('new id_unique is not unique') }
  if (!identical(tes_unnest_filt_df$start, start(tes_unnest_filt))) { stop ('data.table conversion changed ordering of rows') }
    
  # mark exonic TEs and add overlapping genes as concatenated string
  if (!is.null(genes))
  {
    message ('Labeling gene overlaps')
    hits <- findOverlaps(tes_unnest_filt, genes)
    
    tes_unnest_filt <- tes_unnest_filt %>% mutate(exonic = 1:length(.) %in% from(hits))
    
    gene_hits <-
      hits %>% 
          as_tibble %>% 
          mutate(subjectHits = genes[subjectHits]$name) %>%
          group_by(queryHits) %>%
          summarize(name = paste(unique(subjectHits), collapse = '_'))
    tes_unnest_filt$gene <- NA
    tes_unnest_filt[gene_hits$queryHits]$gene <- gene_hits$name
    
    # update id unique
    tes_unnest_filt <- tes_unnest_filt %>% mutate(id_unique = ifelse(exonic, paste(id_unique, gene, sep = '|'), id_unique))
  } else {
    tes_unnest_filt$genes  <- NA
    tes_unnest_filt$exonic <- NA
  }
  
  if (filt_exonic)
  {
    tes_unnest_filt <- tes_unnest_filt %>% filter(!exonic)
  }
  
  # throw unwanted columns
  #res <- tes_unnest_filt %>% select(which(!colnames(mcols(tes_unnest_filt)) %in% c('hmm_st', 'hmm_en', 'orig', 'position_start', 'position_end', 'left_rep')))
  res <- tes_unnest_filt
  return(res)
}