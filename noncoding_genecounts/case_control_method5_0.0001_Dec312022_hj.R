##################
# method5. LoF+Dmis+PLS(indel+snp)
# *.geneCounts.txt
# 20221231
# hojinlee
##################

setwd("/Volumes/hjdrive/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8/results/20221230/")
library(tidyverse)
library(dplyr)
source("/Volumes/hjdrive/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8/scripts/20221230/dom/multigene_process_Dec312022_hj.R")

### set threshold ###
MAF <- 0.0001
asian <- 0.001

### total alleles ###
N <- 17*2 # het

### file name ###
filename <- paste0('thy_n17_case_method5_', format(MAF, scientific=FALSE), '_asian_', format(asian, scientific=FALSE))

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


########## noncoding ###########
cases_noncoding <- read.delim("/Volumes/hjdrive/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8/data/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.variant.table", stringsAsFactors = FALSE, quote="")

### locus split ###
# hail uses locus field instead of chrom and position fields. 
# as "cases" object is made from hail environment, locus field must be split using strsplit func. 
cases_noncoding <- cases_noncoding %>% separate(locus, c("chrom", "position"), ":", remove = FALSE)

### GT and PLDIFF/DP ###
cases_noncoding <- cases_noncoding %>% filter(GT %in% c("0/1", "0|1"))
cases_noncoding$PLDIFF.DEPTH <- as.numeric(cases_noncoding$PLDIFFDP)
cases_noncoding <- cases_noncoding %>% filter(PLDIFF.DEPTH >= 8)

### functional region ###
cases_noncoding <- cases_noncoding %>% filter(!is.na(encode_PLS))

### indel_finder
indel_finder <- function(variant_table) {
  ref_nchar = nchar(variant_table$Ref)
  alt_nchar = nchar(variant_table$Alt)
  indel_index = which((ref_nchar > 1) | (alt_nchar > 1))
  return(variant_table[indel_index,])
}
cases_indel <- indel_finder(cases_noncoding) # indel
cases_noncoding$jarvis <- as.numeric(cases_noncoding$jarvis)
cases_noncoding <- cases_noncoding %>% filter((jarvis >=0.9) & !(is.na(jarvis))) # snp
# snp + indel 
cases_noncoding <- rbind(cases_noncoding, cases_indel)


### filter ###
# filter for spanning variants
cases_noncoding <- cases_noncoding %>% filter(Alt != "*")

# filter for gnomad_combined
cases_noncoding$gnomAD_r2.1.1_merged_all <- as.numeric(cases_noncoding$gnomAD_r2.1.1_merged_all)
cases_noncoding <- cases_noncoding[which(cases_noncoding$gnomAD_r2.1.1_merged_all < MAF | is.na(cases_noncoding$gnomAD_r2.1.1_merged_all) == TRUE), ]

# filter for bravo
cases_noncoding$bravo_freeze8 <- as.numeric(cases_noncoding$bravo_freeze8)
cases_noncoding <- cases_noncoding[which(cases_noncoding$bravo_freeze8 < MAF | is.na(cases_noncoding$bravo_freeze8) == TRUE), ]

# filter for MQ>=40
cases_noncoding$MQ <- as.numeric(cases_noncoding$MQ)
cases_noncoding <- cases_noncoding %>% filter(MQ >= 40)

# filter for DP >= 8
cases_noncoding$DP <- as.numeric(cases_noncoding$DP)
cases_noncoding <- cases_noncoding %>% filter(DP >= 8)

# filter for GQ >= 20
cases_noncoding$GQ <- as.numeric(cases_noncoding$GQ)
cases_noncoding <- cases_noncoding %>% filter(GQ >= 20)

# filter for gnomad asian (genome)
cases_noncoding$gnomad_genome_AF_eas <- as.numeric(cases_noncoding$gnomad_genome_AF_eas)
cases_noncoding <-cases_noncoding[which(cases_noncoding$gnomad_genome_AF_eas < asian | is.na(cases_noncoding$gnomad_genome_AF_eas) == TRUE), ]

### mutigene process using multigene_process_Dec312022_hj.R
cases_noncoding$gnomad_exome_AF_eas <- NA_integer_ # add NA 
cases_noncoding <- cases_noncoding %>% select('Gene_refGeneWithVer', 's','AD', 'chrom', 'position', 'Ref', 'Alt', 'AAChange', 'Func_refGeneWithVer',  'ExonicFunc_refGeneWithVer', 'gnomAD_r2.1.1_merged_all', 'bravo_freeze8', 'gnomad_exome_AF_eas', 'gnomad_genome_AF_eas','MetaSVM_pred', 'REVEL_rankscore','CADD13', 'spliceAI_DS_AG', 'spliceAI_DS_AL', 'spliceAI_DS_DG', 'spliceAI_DS_DL', 'pLi_gnomad_refGeneWithVer', 'oe_lof_upper_gnomad_refGeneWithVer', 'mis_z_gnomad_refGeneWithVer', 'oe_mis_upper_gnomad_refGeneWithVer', 'encode_PLS', 'encode_pELS', 'encode_dELS', "GT", "DP", "GQ", "MQ")
colnames(cases_noncoding) <- c("Gene", "Sample", "AD", "CHROM", "POSITION", "REF", "ALT", "AAChange", "func.region", "VariantClass", "gnomAD", "bravo", "gnomAD.eas.exome", "gnomAD.eas.genome", "MetaSVM", "REVEL", "CADD13", "spliceAI_AG", "spliceAI_AL", "spliceAI_DG", "spliceAI_DL", "pLI", "LOEUF", "mis_z", "MOEUF", "PLS", "pELS", "dELS", "GT", "DP", "GQ", "MQ")

### merge coding and noncoding 
cases <- rbind(cases, cases_noncoding)

### remove duplicate variants 
cases <- cases[which(!duplicated(cases[,c("Sample","POSITION", "CHROM")])),]

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