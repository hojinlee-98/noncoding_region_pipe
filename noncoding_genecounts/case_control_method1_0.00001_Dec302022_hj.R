##################
# method1. syn
# *.geneCounts.txt
# 20221230
# hojinlee
##################

setwd("/Volumes/hjdrive/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8/results/20221230/")
library(dplyr)

### make count file using cases ###
cases <- read.delim("/Volumes/hjdrive/thyroiditis/thyroiditis_v2_LofSpliceai_DmisNonframeshift/data/dom_v2_001/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_0.001_metaSVM_plDiff8_synonymous_variant_removeDup.txt",sep="\t",stringsAsFactors = F)

### GT and PLDIFF/DP ###
cases <- cases %>% filter(GT %in% c("0/1", "0|1"))
cases$PLDIFF.DEPTH <- as.numeric(cases$PLDIFF.DEPTH)
cases <- cases %>% filter(PLDIFF.DEPTH >= 8)
cases <- as_tibble(cases)

### set threshold ###
MAF <- 0.00001
asian <- 0.001

### total alleles ###
N <- 17*2 # het

### file name ###
filename <- paste0('thy_n17_case_syn_', format(MAF, scientific=FALSE), '_asian_', format(asian, scientific=FALSE))

# filter for gnomad_combined
cases$gnomAD_r2.1.1_merged_all <- as.numeric(cases$gnomAD_r2.1.1_merged_all)
cases <-cases[which(cases$gnomAD_r2.1.1_merged_all < MAF | is.na(cases$gnomAD_r2.1.1_merged_all) == TRUE), ]

# filter for bravo
cases$bravo_freeze8 <- as.numeric(cases$bravo_freeze8)
cases <-cases[which(cases$bravo_freeze8 < MAF | is.na(cases$bravo_freeze8) == TRUE), ]

# filter for MQ>=40
cases$MQ <- as.numeric(cases$MQ)
cases <- cases %>% filter(MQ >= 40)

# filter for DP >= 8
cases$DP <- as.numeric(cases$DP)
cases <- cases %>% filter(DP >= 8)

# filter for GQ >= 20
cases$GQ <- as.numeric(cases$GQ)
cases <- cases %>% filter(GQ >= 20)

# filter for gnomad asian
cases$gnomad_exome_AF_eas <- as.numeric(cases$gnomad_exome_AF_eas)
cases <-cases[which(cases$gnomad_exome_AF_eas < asian | is.na(cases$gnomad_exome_AF_eas) == TRUE), ]

### Sum of All Genes
count_table<-function(cases) {
  count.table<-NULL
  gene.list<-unique(cases$Gene.refGeneWithVer) 
  pb = txtProgressBar(0, length(gene.list), style=3)
  for(a in 1:length(gene.list)) {
    gene<-gene.list[a]
    line<-cases[which(cases$Gene.refGeneWithVer==gene.list[a]),] # filtering rows matched with gene.list[a]
    #    count<-nrow(line) ##multi hits
    count<-length(unique(line$SAMPLE)) ##unique hits (count samples for each gene using unique() function)
    temp<-t(as.data.frame(c(as.character(gene),count,N))) # df <= col1:gene name, col2:allele count, col3:total allele
    count.table<-rbind(count.table,temp) # combine all df
    row.names(count.table)[a] <- gene # change rownames 
    setTxtProgressBar(pb, a)
  }
  colnames(count.table) <- c("gene","mut","nor")
  count.table<-as.data.frame(count.table)
  return(count.table)
}

count.table <- count_table(cases)
count.table$mut <- as.numeric(count.table$mut)
count.table$nor <- as.numeric(count.table$nor)
count.table <- count.table %>% arrange(desc(mut))
write.table(count.table, paste(filename,'_geneCounts.txt',sep=''),row.names = F, quote = F, sep="\t")