############################################################
# 1 # Create a GTF file with all the genes and COGs
############################################################
library(data.table)
library(readxl)
#load & merge GTFs
#RefSeq assemblies except for Mycoplasm & Strepto (to allow the locus_tag nomenclature used in reference papers)
GTFs <- list.files("bacteria_essential_genes", pattern = "*.gtf")
GTFs <- paste0("bacteria_essential_genes/",GTFs)
names(GTFs) <- substr(GTFs,26,40)
GTF_list <- lapply(GTFs, fread) 
GTF_list<- mapply(cbind, GTF_list, "assembly"=names(GTF_list), SIMPLIFY=F)
colnombres <- c("ID","DB","feature","start","stop","V6","strand","V8","annotations","assembly")
GTF_list <- lapply(GTF_list, `names<-`, colnombres)
GTF <- do.call("rbind", GTF_list) 
GTF <- GTF[, c(10,3,4,5,7,9)]

#extract annotation attributes
#modified from: https://www.biostars.org/p/272889/#418833
extract_attributes <- function(gtf_attributes, att_of_interest){
  attribute <- c()
  for (i in 1:length(gtf_attributes)){
    att <- unlist(strsplit(gtf_attributes[i], " "))
    if(att_of_interest %in% att){
      attribute[i] <- gsub("\"|;","", att[which(att %in% att_of_interest)+1])
    }else{
      attribute[i] <- NA}
  }
  return(attribute)
}
GTF <- GTF[GTF$feature=="CDS",]
GTF$gene_id <- extract_attributes(GTF$annotations,"gene_id")
GTF$locus_tag <- extract_attributes(GTF$annotations,"locus_tag")
GTF$protein_id <- extract_attributes(GTF$annotations,"protein_id")
GTF$gene <- extract_attributes(GTF$annotations,"gene")
#fix gene names
GTF$gene[setdiff(which(nchar(GTF$gene)>7),grep("PE_",GTF$gene))] <- ""



#add COGs
nog <- read_xlsx("bacteria_essential_genes/out.emapper.annotations.xlsx",sheet=1,skip=2)
nog <- nog[,c(1,5,8,9)]
colnames(nog)[1] <- "protein_id"

GTF <- merge(GTF,nog,by="protein_id")

#take GTF missing 4 letters name from out.emapper (if available)
for (i in 1:length(GTF$gene)){
  if (c(is.na(GTF$gene) & GTF$Preferred_name != "-")[i]==TRUE){
    GTF$gene[i] <- GTF$Preferred_name[i]
  }
}
#write table
fwrite(GTF,"bacteria_essential_genes/all_annotations_dd.GTF", sep=";", quote=FALSE)


############################################################
# 2 # Parse the essential genes tables
############################################################

#2.1# E. coli
coli <- read_xlsx("bacteria_essential_genes/mbo001183726st1.xlsx",sheet=1,skip=1)
coli$Ecoli <- ifelse(coli$Essential==TRUE,"essential",ifelse(coli$`Non-essential`==TRUE,"dispensable","unknown"))
names(coli)[1] <- "gene"
GTF[,11] <-substr(GTF$eggNOG_OGs,1,7)
#merge
essential <- merge(na.omit(GTF[,c(10,11,12)]),coli[,c(1,7)],by="gene", all.x=TRUE)

#2.2# M. tuberculosis
tub <- read_xlsx("bacteria_essential_genes/mbo002173137st3.xlsx",sheet=1,skip=1)
names(tub)[1] <- "locus_tag"
names(tub)[2] <- "gene"
#recode
tub$mtuber <- ifelse(tub$`Final Call`=="ES","essential",ifelse(tub$`Final Call`=="NE","dispensable","intermediate"))
#improve the gene names with the locus_tag
tub <- merge(tub[,c(1,14)],GTF[,c(8:12)], by="locus_tag")
#merge
essential <- merge(essential,tub[,c(2,4:6)],by=c("gene","eggNOG_OGs","Description"), all.x=TRUE)


#2.3# Mycoplasma
myco <- read_xls("bacteria_essential_genes/pnas_0510013103_10013Table2.xls",sheet=1,skip=0)
names(myco)[7] <- "mycoplasm" 
#consider as non-essential those that were disrupted in https://pubmed.ncbi.nlm.nih.gov/16407165/
myco$mycoplasm <- ifelse(is.na(myco$mycoplasm)==TRUE,"essential","dispensable")
#fix locus name
myco$Locus <- gsub("MG","MG_",myco$Locus)
#improve the gene names with the locus_tag
names(myco)[2] <- "locus_tag" 
myco <- merge(myco[,c(2,5,7)],GTF[,c(8:12)], by="locus_tag", all.x=TRUE)
#merge
essential <- merge(essential,myco[,c(3,5:7)],by=c("gene","eggNOG_OGs","Description"), all.x=TRUE)

#2.4# Pseudomonas
pseudo <- read_xlsx("bacteria_essential_genes/pnas.1422186112.sd01.xlsx",sheet=1,skip=1)
names(pseudo)[1:2] <- c("locus_tag","gene")
pseudo$pseudomonas <- "essential"
pseudo <- merge(pseudo[,c(1,22)],GTF[,c(8:12)], by="locus_tag")
essential <- merge(essential,pseudo[,c(4,2,5,6)],by=c("gene","eggNOG_OGs","Description"), all.x=TRUE)

#2.5# Staph
staph_E <- read_xlsx("bacteria_essential_genes/12864_2015_1361_MOESM2_ESM.xlsx",sheet=1,skip=0)
staph_Dom <- read_xlsx("bacteria_essential_genes/12864_2015_1361_MOESM2_ESM.xlsx",sheet=2,skip=0)
staph_NE <- read_xlsx("bacteria_essential_genes/12864_2015_1361_MOESM2_ESM.xlsx",sheet=3,skip=0)
staph_E$staph <- "essential"
staph_Dom$staph <- "essential"
staph_NE$staph <- "dispensable"
staph <- rbind(staph_E[,c(1,4)],staph_Dom[,c(1,4)])
staph <- rbind(staph,staph_NE[,c(1,4)])
names(staph)[1] <- "locus_tag" 
staph$locus_tag <- substr(staph$locus_tag,1,13)
staph <-merge(staph,GTF[,c(8:12)], by="locus_tag", all.x=TRUE)
essential <- merge(essential,staph[,c(4,2,5,6)],by=c("gene","eggNOG_OGs","Description"), all.x=TRUE)

#2.6# Strepto
strep <- read_xlsx("bacteria_essential_genes/strepto_lebreton.xlsx",sheet=1,skip=0)
strep$essentiality <- ifelse(strep$essentiality=="E","essential",ifelse(strep$essentiality=="C"|strep$essentiality=="NC", "unknown","dispensable"))
names(strep)[3] <- "streptococcus"
strep$locus_tag <- gsub("Spy_","Spy",strep$locus_tag)
strep <- merge(strep[,c(1,3)],GTF[,c(8:12)], by="locus_tag")
essential <- merge(essential,strep[,c(4,2,5,6)],by=c("gene","eggNOG_OGs","Description"), all.x=TRUE)


#2.7# Synechococcus
Syne <- read_xlsx("bacteria_essential_genes/pnas.1519220112.sd03.xlsx",sheet=1,skip=0)

names(Syne)[c(3,10)] <- c("locus_tag" ,"Synechococcus")
Syne <-merge(Syne[,c(3,8,10)],GTF[,c(8:12)], by="locus_tag", all.x=TRUE)
Syne$gene <- apply(Syne[,c(2,5)], 1, function(x) x[!is.na(x)][1])
Syne$Synechococcus[Syne$Synechococcus=="non-essential"] <- "dispensable"
essential <- merge(essential,Syne[,c(3,5:7)],by=c("gene","eggNOG_OGs","Description"), all.x=TRUE)
#remove duplicates and missing
essential <- unique(essential)
essential <- essential[-which(is.na(essential$Ecoli) & is.na(essential$mtuber) & is.na(essential$mycoplasm) & is.na(essential$staph) & is.na(essential$streptococcus) & is.na(essential$Synechococcus) ), ]

#merge by COG 
kk <- essential %>% 
  group_by(eggNOG_OGs) %>% 
  summarise(gene = paste(na.omit(gene), collapse = ","), Description = paste(Description, collapse=","),
            Ecoli = paste(na.omit(Ecoli), collapse=","),
            mtuber = paste(na.omit(mtuber), collapse=","),
            mycoplasm = paste(na.omit(mycoplasm), collapse=","),
            pseudomonas = paste(na.omit(pseudomonas), collapse=","),
            staph = paste(na.omit(staph), collapse=","),
            streptococcus = paste(na.omit(streptococcus), collapse=","),
            Synechococcus = paste(na.omit(Synechococcus), collapse=","),
            )  


kk <- kk %>%
  mutate(essential = rowSums(
    sapply(select(., Ecoli, mtuber, mycoplasm, pseudomonas, staph, streptococcus, Synechococcus ),
           function(x) grepl("essential", x, ignore.case = TRUE)
    ))  
  )
kk <- kk %>%
  mutate(dispensable = rowSums(
    sapply(select(., Ecoli, mtuber, mycoplasm, pseudomonas, staph, streptococcus, Synechococcus ),
           function(x) grepl("dispensable", x, ignore.case = TRUE)
    )) 
  )

#write table
fwrite(kk,"bacteria_essential_genes/all_essentials.csv", sep=";", quote=FALSE)

