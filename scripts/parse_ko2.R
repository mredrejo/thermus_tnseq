#this code is not very efficient and it takes a while
#not intended to reproduce in every knitting, just run once and save the tables
library(data.table)
library(tidyverse)
library(dplyr)

ko <- fread("scripts/ko.csv",sep=";",header=TRUE)

datos <- read.csv2("tnSeq_pangenome2.csv")
trna <- read.csv2("scripts/tRNAs.csv",header = FALSE)
datos$KEGG_ko <- gsub("ko:","",datos$KEGG_ko)
#tRNAs annotation is not very good, so we corrected for statistics
#all tRNAs not annotated were masked as K14228 (tRNA-leu)
datos$KEGG_ko[datos$Genes %in% trna$V1] <- "K14228"

duplicate_ko_rows <- function(data) {
  data %>% 
    mutate(reaction_ko = strsplit(KEGG_ko, ",")) %>%  # Split ko annotations
    unnest(reaction_ko) %>%                         # Unnest each annotation into separate rows
    select(-KEGG_ko)                            # Keep original columns except split annotations
}
datos_unnested <- duplicate_ko_rows(datos)

stats <- summarise(group_by(datos_unnested,reaction_ko,cluster,tt_i4c5,thermus_i4c5,thermaceae_i4c5,deinococcota_i4c5),count =n())
stats <- plyr::join(stats,ko[,c(1:3,5)])


