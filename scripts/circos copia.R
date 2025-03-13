library(BioCircos)
library(circlize)

#https://www.royfrancis.com/beautiful-circos-plots-in-r/

#pTT27  
# coverage pTT27
cov <- list(data.frame(matrix(NA, nrow = 1, ncol = 4)),data.frame(matrix(NA, nrow = 1, ncol = 4)),data.frame(matrix(NA, nrow = 1, ncol = 4)),data.frame(matrix(NA, nrow = 1, ncol = 4)))
for (v in 1:4){
  for (i in 1:233){
    cov[[v]][i,]<- coverage %>% filter(Sample == samples[v]) %>% slice((((i-1)*1000)+1):(i*1000)) %>% summarize(chr="AE017222.1", start=(((i-1)*1000)+1),end=(i*1000),value1 = mean(Pos+Neg))
  }
}
save(cov, file = "cov_pTT27.Rdata")
#circos plot
circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=0, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=c("AE017222.1"), 
                  xlim=matrix(c(0, 232605), ncol=2))

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.6, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.06)

brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)*10^5
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^5, 1), labels.cex=0.6, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

#cov tracks
circos.genomicTrack(data=cov[[1]], panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", col="#1F78B4", lwd=1)
},track.height=0.08, bg.border=F)

circos.track(factors=cov[[2]]$X1, x=cov[[2]]$X2, y=cov[[2]]$X4, panel.fun=function(x, y) {
  circos.lines(x, y, col="#1F78B4", lwd=1)
}, ylim=range(cov[[2]]$X4), track.height=0.08, bg.border=F)

circos.track(factors=cov[[3]]$X1, x=cov[[3]]$X2, y=cov[[3]]$X4, panel.fun=function(x, y) {
  circos.lines(x, y, col="#33A02C", lwd=1)
}, ylim=range(cov[[3]]$X4), track.height=0.08, bg.border=F)


circos.track(factors=cov[[4]]$X1, x=cov[[4]]$X2, y=cov[[4]]$X4, panel.fun=function(x, y) {
  circos.lines(x, y, col="#33A02C", lwd=1)
}, ylim=range(cov[[4]]$X4), track.height=0.08, bg.border=F)


#chr

# coverage chr
cov_chr <- list(data.frame(matrix(NA, nrow = 1, ncol = 4)),data.frame(matrix(NA, nrow = 1, ncol = 4)),data.frame(matrix(NA, nrow = 1, ncol = 4)),data.frame(matrix(NA, nrow = 1, ncol = 4)))
for (v in 1:4){
  for (i in 1:1895){
    cov_chr[[v]][i,]<- coverage %>% filter(Sample == samples[v]) %>% slice((((i-1)*1000)+1):(i*1000)) %>% summarize(chr="AE017221.1", start=(((i-1)*1000)+1),end=(i*1000),value1 = mean(Pos+Neg))
  }
}
save(cov_chr, file = "cov_chr.Rdata")

#circos plot
circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=0, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=c("AE017221.1"), 
                  xlim=matrix(c(0, 1894877), ncol=2))

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.6, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.06)

brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)*10^5
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^5, 1), labels.cex=0.6, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

#cov tracks
circos.genomicTrack(data=cov_chr[[1]], panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", col="#1F78B4", lwd=1)
},track.height=0.08, bg.border=F)

circos.track(factors=cov_chr[[2]]$X1, x=cov_chr[[2]]$X2, y=cov_chr[[2]]$X4, panel.fun=function(x, y) {
  circos.lines(x, y, col="#1F78B4", lwd=1)
}, ylim=range(cov_chr[[2]]$X4), track.height=0.08, bg.border=F)

circos.track(factors=cov_chr[[3]]$X1, x=cov_chr[[3]]$X2, y=cov_chr[[3]]$X4, panel.fun=function(x, y) {
  circos.lines(x, y, col="#33A02C", lwd=1)
}, ylim=range(cov_chr), track.height=0.08, bg.border=F)


circos.track(factors=cov_chr[[4]]$X1, x=cov_chr[[4]]$X2, y=cov_chr[[4]]$X4, panel.fun=function(x, y) {
  circos.lines(x, y, col="#33A02C", lwd=1)
}, ylim=range(cov_chr[[4]]$X4), track.height=0.08, bg.border=F)


