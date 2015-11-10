library(ggplot2)
library(cowplot)
library(reshape)
library(magrittr)

Jurkat <- 
structure(c(11, 14, 16, 17, 19, 13, 17, 22, 27, 31, 16, 23, 28, 
35, 42, 22, 32, 41, 53, 69, 27, 45, 59, 76, 114, 37, 70, 102, 
160, 301, 1793, 2965, 3915, 5968, 9877), .Dim = c(5L, 7L), .Dimnames = list(
    c("reg.FDR.1%", "reg.FDR.3%", "reg.FDR.5%", "reg.FDR.10%", 
    "reg.FDR.20%"), c("prom.FDR.1%", "prom.FDR.3%", "prom.FDR.5%", 
    "prom.FDR.10%", "prom.FDR.20%", "overlap", "total")))

GM12878 <- 
structure(c(72, 86, 91, 101, 107, 86, 109, 118, 133, 145, 96, 
124, 135, 156, 172, 106, 146, 164, 197, 229, 114, 162, 187, 237, 
295, 136, 214, 270, 407, 635, 2592, 4310, 5633, 8599, 13620), .Dim = c(5L, 
7L), .Dimnames = list(c("reg.FDR.1%", "reg.FDR.3%", "reg.FDR.5%", 
"reg.FDR.10%", "reg.FDR.20%"), c("prom.FDR.1%", "prom.FDR.3%", 
"prom.FDR.5%", "prom.FDR.10%", "prom.FDR.20%", "overlap", "total"
)))

efet <- function(y,nm) {
# Called in [i,j] | called in i=20%
  y <- y[-5,]/matrix(y[5,],ncol=ncol(y),nrow=nrow(y)-1,byrow=TRUE)  
# Called in [i,j] | called in i=20%

  z <- (y[,-7]/y[,7])  %>% as.data.frame()
  z$reg.fdr <- c(1,3,5,10)
  z$Cell.line <- nm
  return(z)
}

y1 <- efet(GM12878,"GM12878")
y2 <- efet(Jurkat,"Jurkat")
y <- rbind(y1,y2)
m <- melt(y,c("reg.fdr","Cell.line"))
#m$reg.fdr <- as.numeric(gsub("reg.FDR.|%","",m$reg.fdr))
  m$prom.fdr <- sub("prom.FDR.","",m$variable)
  m$prom.fdr[m$prom.fdr=="overlap"] <- "--"
m$prom.fdr <- factor(m$prom.fdr,levels=c("1%","3%","5%","10%","20%","--"))
p <- ggplot(m,aes(x=reg.fdr,y=value,col=prom.fdr)) +
  geom_hline(yintercept=1,linetype="dashed") +
  geom_point() +
  geom_path() +
  ylab("Fold enrichment") +
  xlab("Region FDR (%)") +
  facet_grid(. ~ Cell.line,scales="free") + scale_colour_discrete("Promoter\nFDR") +
  background_grid(major = 'xy', minor = "none") + panel_border() 
print(p)

tiff("FDR_overlap.tiff",height=4,width=6,units="in",res=300)
print(p)
dev.off()
