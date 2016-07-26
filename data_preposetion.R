try http:// if https:// URLs are not supported 
source("https://bioconductor.org/biocLite.R")
biocLite("PAA")
# add library for normalisation data
biocLite("vsn")
# add librery for Protein Microarray Data Analysis
library("PAA")
install.packages("Rcpp")

library(limma)
#Load GPR files
my_data <- loadGPR(gpr.path="Data", targets.path="Data/targets.txt",array.type = "ProtoArray")
# delete background
bg_matrix <- backgroundCorrect(my_data, method="normexp", normexp.method="saddle")


library("vsn")
# normalize data
plotMAPlots(elist=elist, idx=10)
elist_1 <- normalizeArrays(elist=elist, method="cyclicloess", cyclicloess.method="fast")


library("pheatmap")
# take matrix part
E <- elist_1$E
head(E)

# take names without trash
gn = sub(pattern = 'Hs~', rep = "", x = elist_1$genes$Name)
gn = sub(pattern = '~.*', rep = "", gn)
gn = sub(pattern = '.*:', rep = "", gn)
head(gn)

#add names to matrix
row.names(E) <- gn
aggr <- aggregate(E, list(gn), mean)
row.names(aggr) <- aggr$Group.1
aggr <- aggr[-1]
