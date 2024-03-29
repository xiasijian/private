rm(list = ls())
gc()

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)

df = read.csv(file = "./mm9_ucsc_tss20240329.bed",sep = "\t")

transcripts <- df %>%
  dplyr::select(X.name, chrom, txStart, txEnd, strand, name2) %>%
  dplyr::rename(transcript = X.name, gene = name2)

table(duplicated(transcripts$transcript))
# FALSE  TRUE 
# 25922   518

## get TSS
# if on - strand, start and end must be flipped
tss <- transcripts %>% 
  dplyr::mutate(tss_start = ifelse(strand == '+', txStart, txEnd - 1),
                tss_end = tss_start + 1) %>%
  dplyr::select(-txStart, -txEnd) %>%
  unique()


## eliminate redundant TSS, so that at a given position and strand, there is no more than one TSS
tss$id <- with(tss, paste(chrom, tss_start, strand, sep = ':'))
table(duplicated(tss$id,fromLast =T))

tss_dup =  tss[!duplicated(tss$id,fromLast =T),]


## keep chr1:19, chrX and chrY 
tss_dup_filt1 = tss_dup %>% subset(chrom %in% c(paste0("chr",1:19),"chrX","chrY"))

## order
tss_final <- tss_dup_filt1[order(tss_dup_filt1$chrom, tss_dup_filt1$tss_start, tss_dup_filt1$strand),]

## export bed file 

write.table(tss_final[,c("chrom","tss_start","tss_end","strand","gene","transcript")],
            file = paste0("mm9_ucsc_tss20240329_clean.bed"),sep = "\t",col.names = F,row.names = F,quote = F)
