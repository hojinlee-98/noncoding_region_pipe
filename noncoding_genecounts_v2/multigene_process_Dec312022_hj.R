### rep_row
# replicate dataframe only consisted of a line 
# called by expand_df
# @ line  dataframe only consisted of a line 
# @ time  number of iterations
rep_row<- function(line, time, gene) {
  line.seed <- list()
  for (ii in 1:time) {line.seed <- rbind(line.seed, line)}
  # add split.gene[[i]] to split.refGeneWithVer field
  line.seed$split.Gene <- gene
  return(line.seed)
}

### expand_df
# merge two datafames
# one does not have multi-genes,the other has multi-gen
# @ df  case or case_noncoding
expand_df <- function(df) {
  ### split
  df.list <- list() # empty list
  # df not having multi-genes
  df.list$temp1 <- df %>% filter(!str_detect(Gene, "\\\\x3b")) %>% mutate(split.Gene = Gene) 
  # df having multi-genes
  df.list$temp2 <- df %>% filter(str_detect(Gene, "\\\\x3b"))
  # every each row has has sub-list included in split.gene
  split.gene <- str_split(df.list$temp2$Gene, "\\\\x3b") # split.gene is list
  df.list$df.seed <- list() # empty list
  
  for (i in 1:nrow(df.list$temp2)) {
    line <- df.list$temp2[i,] # get each line from temp2 df
    time <- length(split.gene[[i]]) # get iteration time for rep_row func
    df.list$df.seed <- rbind(df.list$df.seed, rep_row(line, time, split.gene[[i]]))
  }
  df.list$temp2 <- NULL
  df <- do.call(rbind, df.list) # merge the two dataframes
  return(df)
}