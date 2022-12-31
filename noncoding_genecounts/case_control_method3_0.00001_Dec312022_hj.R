##################
# method3. LoF+Dmis
# *.geneCounts.txt
# 20221231
# hojinlee
##################

setwd("/Volumes/hjdrive/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8/results/20221230/")
library(tidyverse)
library(dplyr)
source("/Volumes/hjdrive/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8/scripts/20221230/dom/multigene_process_Dec312022_hj.R")

### set threshold ###
MAF <- 0.00001
asian <- 0.001

### total alleles ###
N <- 17*2 # het

### file name ###
filename <- paste0('thy_n17_case_method3_', format(MAF, scientific=FALSE), '_asian_', format(asian, scientific=FALSE))

########### coding ############
### make count file using cases ###
cases <- read.delim("/Volumes/hjdrive/thyroiditis/thyroiditis_v2_LofSpliceai_DmisNonframeshift/data/dom_v2_001/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_0.001_metaSVM_plDiff8_dom_variant_removeDup.txt",sep="\t",stringsAsFactors = F)

### GT and PLDIFF/DP ###
cases <- cases %>% filter(GT %in% c("0/1", "0|1"))
cases$PLDIFF.DEPTH <- as.numeric(cases$PLDIFF.DEPTH)
cases <- cases %>% filter(PLDIFF.DEPTH >= 8)
cases <- as_tibble(cases)

### Dmis + lof ###
cases1 <- cases %>% filter(ExonicFunc.refGeneWithVer == "nonsynonymous_SNV" & MetaSVM_pred == "D")
cases2 <- cases[which(cases$ExonicFunc.refGeneWithVer == "stoploss"),]
cases3 <- cases[which(cases$ExonicFunc.refGeneWithVer %in% c("nonframeshift_deletion", "nonframeshift_insertion")),]
cases4 <- cases[which(cases$EFFECTIVE.VARIANT.TYPE == "lof"),]
# spliceAI
cname <- c("spliceAI_DS_AG.1", "spliceAI_DS_AL.1", "spliceAI_DS_DG.1", "spliceAI_DS_DL.1")
cases[,cname] <- apply(cases[,cname], MARGIN = 2, as.numeric)
cases5 <- cases %>% filter(EFFECTIVE.VARIANT.TYPE == "deleterious") %>%
  filter((spliceAI_DS_AG.1 >= 0.8) | (spliceAI_DS_AL.1 >= 0.8) | (spliceAI_DS_DG.1 >= 0.8) | (spliceAI_DS_DL.1 >= 0.8))
cases <- rbind(cases1, cases2, cases3, cases4, cases5)
# remove duplicate variants 
cases <- cases[which(!duplicated(cases[,c("SAMPLE","POSITION","CHROM")])),]

### filter ###
# filter for spanning variants
cases <- cases %>% filter((ALTERNATE.ALLELES) != "*")

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

### subset table 
cases[c("encode_PLS", "encode_pELS", "encode_dELS")] <- NA_character_ # add empty field
cases <- cases %>% select('Gene.refGeneWithVer', 'SAMPLE','AD', 'CHROM', 'POSITION', 'REFERENCE.ALLELE', 'ALTERNATE.ALLELES', 'AAChange.refGeneWithVer', 'Func.refGeneWithVer', 'ExonicFunc.refGeneWithVer', 'gnomAD_r2.1.1_merged_all', 'bravo_freeze8', 'gnomad_exome_AF_eas', 'gnomad_genome_AF_eas', 'MetaSVM_pred', 'REVEL_rankscore', 'CADD13_PHRED.1', 'spliceAI_DS_AG.1', 'spliceAI_DS_AL.1', 'spliceAI_DS_DG.1', 'spliceAI_DS_DL.1', 'pLi_gnomad.refGeneWithVer', 'oe_lof_upper_gnomad.refGeneWithVer', 'mis_z_gnomad.refGeneWithVer', 'oe_mis_upper_gnomad.refGeneWithVer', "encode_PLS", "encode_pELS", "encode_dELS", "GT", "DP", "GQ", "MQ")
colnames(cases) <- c("Gene", "Sample", "AD", "CHROM", "POSITION", "REF", "ALT", "AAChange", "func.region", "VariantClass", "gnomAD", "bravo", "gnomAD.eas.exome", "gnomAD.eas.genome", "MetaSVM", "REVEL", "CADD13", 'spliceAI_AG', 'spliceAI_AL', 'spliceAI_DG', 'spliceAI_DL', "pLI", "LOEUF", "mis_z", "MOEUF", "PLS", "pELS", "dELS", "GT", "DP", "GQ", "MQ")

### mutigene process using multigene_process_Dec312022_hj.R
cases <- expand_df(cases)

### Sum of All Genes ###
count_table <- function(cases) {
  count.table <- NULL
  gene.list <- unique(cases$split.Gene)
  pb = txtProgressBar(0, length(gene.list), style=3)
  for(a in 1:length(gene.list)) {
    gene <- gene.list[a]
    line <- cases[which(cases$split.Gene == gene.list[a]),] # filtering rows matched with gene.list[a]
    count <- length(unique(line$Sample)) ##unique hits (count samples for each gene using unique() function)
    temp <- t(as.data.frame(c(as.character(gene),count,N))) # df <= col1:gene name, col2:allele count, col3:total allele
    count.table <- rbind(count.table,temp) # combine all df
    row.names(count.table)[a] <- gene # change rownames 
    setTxtProgressBar(pb, a)
  }
  colnames(count.table) <- c("gene","mut","nor")
  count.table<-as.data.frame(count.table)
  return(count.table)
}

count.table <- count_table(cases)
cname <- c("mut", "nor")
count.table[,cname] <- apply(count.table[,cname], MARGIN = 2, as.numeric)
count.table <- count.table %>% arrange(desc(mut))

write.table(count.table, paste(filename,'_geneCounts.txt',sep=''),row.names = F, quote = F, sep="\t")