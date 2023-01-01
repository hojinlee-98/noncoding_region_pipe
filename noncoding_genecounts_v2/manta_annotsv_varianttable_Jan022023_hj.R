##############
# make sv variant table from manta
# hojinlee
# 20230102
##############


library(dplyr)
library(tidyverse)


manta_concat <- function(filelist) {
  # initialize list
  annotsv_list <- list()
  for (pathtofile in filelist) {
    # data load
    annotsv_table_temp <- read.table(pathtofile, sep = "\t", quote = "", header = T, stringsAsFactors = F)
    filename <- basename(pathtofile)
    # get sampleID
    sampleID <- str_split(filename, "_")[[1]][1]
    # when some sampleID are composed of only numeric, ex) 178898
    # they are changed to "X${smpleID}" automatically in R dataframe column name
    # so, we need to see the data type of sampleID 
    detector <- as.numeric(sampleID)
    # if NA is assigned to detector, 
    # it means sampleID has only alphabet.
    # on the other hands, if there are not any problems, as.numeric returns NA, 
    # the sampleID is totally composed of numeric.
    # then we need to paste "X" chr before the sampleID. 
    if (!is.na(detector)) {
      sampleID <- paste0("X", sampleID)
    }
    # split format values to each column
    # sampleID field automatically removed, because of "keep" parameter
    annotsv_table_temp <- annotsv_table_temp %>% separate(sampleID, into=c("GT", "FT", "GQ", "PL", "PR", "SR"), sep=":")
    # split ref,alt to each column
    annotsv_table_temp <- annotsv_table_temp %>% separate(PL, into=c("PL_homref", "PL_het", "PL_homalt"), sep=",") %>% separate(PR, into=c("PR_ref", "PR_alt"), sep=",") %>% separate(SR, into=c("SR_ref", "SR_alt"), sep=",")
    
    if (!is.na(detector)) {
      sampleID <- paste0("X", sampleID)
      sampleID <- gsub(pattern = "X", replacement = "", sampleID)
    }
    # add sample field
    annotsv_table_temp <- annotsv_table_temp %>% mutate(Sample = sampleID)
    
    ### calculate fraction of alt alleles
    depthcol <- c("PR_ref", "PR_alt", "SR_ref", "SR_alt")
    annotsv_table_temp[,depthcol] <- apply(annotsv_table_temp[,depthcol], 2, function(x) as.numeric(x))
    # fr=alt/(ref+alt)
    annotsv_table <- annotsv_table_temp %>% mutate(PR_fraction=PR_alt/(PR_ref+PR_alt)) %>% mutate(SR_fraction=SR_alt/(SR_ref+SR_alt))
    
    annotsv_list[[pathtofile]] <- annotsv_table
  }
  
  # variant table for structural variants
  allsamples_annotsv <- do.call(rbind, annotsv_list)
  return(allsamples_annotsv)
}

#####  Gene units  #####
### Sum of All Genes
# @ Gene_name : column name include gene symbols
# @ cases : variant table
count_table <- function(cases, Gene_name) {
  count.table <- NULL
  gene.list <- unique(cases[[Gene_name]]) # make vector including uniq Gene_name
  pb = txtProgressBar(0, length(gene.list), style=3)
  for(a in 1:length(gene.list)) { # explore gene.list 
    gene <- gene.list[a] # gene name assign 
    line <- cases[which(cases[[Gene_name]]==gene.list[a]),] # filtering rows matched with gene.list[a]
    #    count <- nrow(line) ##multi hits
    count <- length(unique(line$Sample)) ##unique hits (count samples for each gene using unique() function)
    temp <- t(as.data.frame(c(as.character(gene),count,N))) # df <= col1:gene name, col2:allele count, col3:total allele
    count.table <- rbind(count.table,temp) # combine all df
    row.names(count.table)[a] <- gene # change rownames 
    setTxtProgressBar(pb, a)
  }
  colnames(count.table) <- c("gene","mut","nor")
  count.table <- as.data.frame(count.table)
  return(count.table)
}