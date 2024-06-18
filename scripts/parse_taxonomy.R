#1: Load & prepare the data
library(data.table)
deinococcota <- fread("pangenomes/deinococcota_dataset.txt") #created for dataformat: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/
names(deinococcota) <- c("Accession","Organism","TaxId")


#2: Get full taxonomy
library(taxonomizr)
prepareDatabase('accessionTaxa.sql',getAccessions=FALSE)
taxonomy <- getTaxonomy(unique(deinococcota$TaxId),'accessionTaxa.sql')
dataset <- as.data.frame(taxonomy)
dataset$TaxId <- as.numeric(row.names(dataset))
deinococcota_full <- merge(deinococcota,dataset)

#3: save full table
write.table(deinococcota_full, "pangenomes/deinococcota_full.tsv", quote=FALSE,row.names = FALSE, sep="\t")

#4: save tables for PPanGGOLIN

#Deinococcota

deinococcota <- data.frame(deinococcota_full$Accession, paste0("bakta_results/",deinococcota_full$Accession,"/",deinococcota_full$Accession,".gbff"))
write.table(deinococcota, "pangenomes/deinococcota.gbff.list", quote=FALSE,col.names =FALSE,row.names = FALSE, sep="\t")

#thermaceae (=Thermales)
tmp <- subset(deinococcota_full,deinococcota_full$order=="Thermales")
thermaceae <- data.frame(tmp$Accession, paste0("bakta_results/",tmp$Accession,"/",tmp$Accession,".gbff"))
write.table(thermaceae, "pangenomes/thermaceae.gbff.list", quote=FALSE,col.names =FALSE, row.names = FALSE, sep="\t")
#Thermus
tmp <- subset(deinococcota_full,deinococcota_full$genus=="Thermus")
thermus <- data.frame(tmp$Accession, paste0("bakta_results/",tmp$Accession,"/",tmp$Accession,".gbff"))
write.table(thermus, "pangenomes/thermus.gbff.list", quote=FALSE,col.names =FALSE,row.names = FALSE, sep="\t")
#T. thermophilus
tmp <- subset(deinococcota_full,deinococcota_full$species=="Thermus thermophilus")
tt <- data.frame(tmp$Accession, paste0("bakta_results/",tmp$Accession,"/",tmp$Accession,".gbff"))
write.table(tt, "pangenomes/tt.gbff.list", quote=FALSE,col.names =FALSE,row.names = FALSE, sep="\t")





