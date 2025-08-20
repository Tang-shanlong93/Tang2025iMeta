if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")

library(mixOmics)
library(readxl)

Protein <- read.csv("D:/R workpalce/pro.csv", header = T, row.names = 1)

MicroL <- read.csv("D:/R workpalce/Microl.csv", header = T, row.names = 1)
MicroM <- read.csv("D:/R workpalce/Microm.csv", header = T, row.names = 1)

Metabolite <- read.csv("D:/R workpalce/metab.csv", header = T, row.names = 1)

X <- list(MicrobiomeL = MicroL, MicrobiomeM = MicroM, Metabolomics = Metabolite,
          Proteomics = Protein)


Y<-factor(c('CON','CON','CON','CON','CON','CON','XOS','XOS','XOS','XOS','XOS','XOS'))
summary(Y)

list.keepX <- list(MicrobiomeL = c(15, 15), MicrobiomeM = c(15, 15), Metabolomics = c(5, 15), 
                   Proteomics = c(5, 15))

MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)

plotDiablo(MyResult.diablo, ncomp = 1)
circosPlot(MyResult.diablo, cutoff=0.7,line = T)

#plotLoadings(MyResult.diablo, contrib = "max")
plotLoadings(MyResult.diablo, comp = 2, contrib = "max")


plotIndiv(MyResult.diablo) ## sample plot

plotVar(MyResult.diablo) ## variable plot

plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2),
          title = 'BRCA with DIABLO')

plotVar(MyResult.diablo, var.names = c(FALSE, TRUE),
        legend=TRUE, pch=c(16,16))

