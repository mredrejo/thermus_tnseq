#this code is not very efficient and it takes a while
#not intended to reproduce in every knitting, just run once and save the tables
library(data.table)

ko <- fread("scripts/ko.csv",sep=";",header=TRUE)

datos <- read.csv2("tnSeq.csv")
trna <- read.csv2("scripts/tRNAs.csv",header = FALSE)
datos$KEGG_ko <- gsub("ko:","",datos$KEGG_ko)
#tRNAs annotation is not very good, so we corrected for statistics
#all tRNAs not annotated were masked as K14228 (tRNA-leu)
datos$KEGG_ko[datos$Genes %in% trna$V1] <- "K14228"
stats <- as.data.frame(table(datos$KEGG_ko,datos$cluster))
categories <- c()
for (i in 1:nrow(ko)){
  categories <- grep(ko[i,5],stats$Var1)
  categories_D <- subset(stats[categories,], Var2=="Highly permissive")
  categories_R <- subset(stats[categories,], Var2=="Intermediate")
  categories_T <- subset(stats[categories,], Var2=="Less permissive")
  ko$Dispensable[i] <- sum(categories_D$Freq)
  ko$Relevant[i] <- sum(categories_R$Freq)
  ko$Top[i] <- sum(categories_T$Freq)
}

kegg <- cbind(ko[,c(1,2,3,5)],stack(ko[,9:11]))

#xtabs
xkegg <- as.data.frame(xtabs(values~a_class+b_class+pathway+ind,data=kegg))
xkegg <- xkegg[!xkegg$Freq==0,]

for (i in 1:nrow(xkegg)){
  xkegg$ratio_b[i] <- xkegg$Freq[i]*100/sum(xkegg$Freq[xkegg$b_class==xkegg$b_class[i]])
}
write.table(xkegg,"xkegg.txt",quote=FALSE,sep=";",col.names = TRUE,row.names=FALSE)

