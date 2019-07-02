devtools::load_all('~/tgdata/src/Reputils/')
devtools::load_all('~/tgdata/src/Repsc/')
library(BSgenome.Hsapiens.UCSC.hg38)

# Introduce sequence errors

bams    <- grep('deduplicated', dir('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/aligned/chunked_bams/', pattern = '.bam$', full.names = TRUE), value = TRUE)
tes     <- importRMSK('~/tgdata/src/Repsc/inst/extdata/hg38.fa.out.gz', curate = TRUE)

mutateReads(BSgenome   = Hsapiens,
            reads      = bams,
            paired     = TRUE,
            tes        = tes,
            n_reads    = 1e6,
            outdir     = '~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/fastqs/',
            error_rate = c(0, 0.25, 0.5, 1, 2, 4, 8))

# STAR alignment
 
r1 <- dir('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/fastqs/', pattern = 'perc_error_R1', full.names = TRUE)
r2 <- dir('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/fastqs/', pattern = 'perc_error_R2', full.names = TRUE)

setwd('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/aligned/')
cmds <- list()
for (i in 1:length(r1))
{
  cmds[[i]] <- 
    glue("system('/net/mraid14/export/data/tools/star/STAR-2.5.1b/source/STAR --runThreadN 10 --genomeDir /net/mraid14/export/data/users/davidbr/tools/cellranger/references/refdata-cellranger-GRCh38-1.2.0/star/ --outSAMtype BAM SortedByCoordinate --outFilterMultimapScoreRange 2 --outFileNamePrefix {gsub('.fastq', '.', basename(r1[i]))} --outSAMunmapped Within --readFilesIn {r1[i]} {r2[i]}')")
}

# create new empty env and fill with relevant
empty <- new.env(parent = emptyenv())

# distribute, compute and combine
res <- 
  gcluster.run3(command_list = cmds,  
                  max.jobs = 50, 
                  envir = empty, 
                  io_saturation = FALSE)

# Create scSet

# path to Gencode gtf file (provided)
gene_path <- system.file(package = 'Repsc', 
                        'extdata', 
                        'gencode.v29.annotation.gtf.gz')

# path to RepeatMasker hg38 repeat annotation (provided)
rmsk_path <- system.file(package = 'Repsc', 
                         'extdata', 
                         'hg38.fa.out.gz')
                         
# creating the scSet
sc <- createScSet(genome   = Hsapiens,
                  protocol = 'fiveprime',
                  tes      = rmsk_path,
                  genes    = gene_path)

# Count reads

# path to bam files containing mapped reads                         
bam_paths <- 
  grep('error', 
       dir('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/aligned/', 
           pattern = 'bam$', 
           full.names = TRUE), 
       value = TRUE)
                 
# create a data.frame specifying import parameters                 
bam_df    <- data.frame(paths   = bam_paths,
                        paired  = TRUE,       # use FALSE for single-end libraries
                        mate    = 'first',    # only imports the first mate of properly aligned read pairs, set to NA when using single-end libraries
                        barcode = c('0.5', '0', '1', '2', '4', '8'),       # 10x barcode included in BAM flag
                        stringsAsFactors = FALSE)
sc <- addCounts(sc,
                reads    = bam_df,
                bin_size = 20,
                msa_dir  = NULL)

# Mapping QC

ltr12c <- 
  counts(sc) %>% 
    filter(name == 'LTR12C') %>%
    group_by(id_unique, barcode) %>%
    summarize(counts = sum(n)) %>%
    ungroup
                



                  
  
  
  
  
bams  <- grep('error', dir('~/davidbr/proj/epitherapy/data/h1299/10x/dacsb/aligned/', pattern = 'bam$', full.names = TRUE), value = TRUE)
reads <- importBAM(bams,
                  paired  = TRUE, 
                  mate    = 'first',
                  anchor  = 'fiveprime',
                  multi   = FALSE,
                  barcode = c('0.25', '0.5', '0', '1', '2', '4', '8'))

# add overlapping TE locus
reads_anno <- join_overlap_left_directed(reads, tes)   

reads_df <- reads_anno %>% mutate(read_id = names(reads_anno)) %>% as_tibble()        

n_reads    <- reads_dt %>% count(barcode)
fam_counts <- reads_dt[which(barcode == 0), ] %>% count(name)

bla = 
  reads_dt[, same := abs(start - start[barcode == 0]) < 10, by = 'read_id'
    ][, .(correctly_mapped = sum(same, na.rm = TRUE) / length(same) * 100, repclass = repclass[1]), by = c('barcode', 'name')]
bla %>% 
  filter(name %in% (fam_counts %>% filter(n > 100) %>% pull(name))) %>%
  ggplot(aes(x = as.numeric(barcode), y = correctly_mapped, group = name, col = repclass, alpha = 0.25)) + 
    geom_point() +
    geom_line() +
    facet_wrap(~repclass) +
    theme(legend.position = 'none')
     
                  

                  