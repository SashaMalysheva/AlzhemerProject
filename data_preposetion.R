library(PAA)
#Load GPR files
my_data <- loadGPR(gpr.path="C:/Users/alex/alzhemerproject/Alzhemerproject/data", targets.path="C:/Users/alex/alzhemerproject/Alzhemerproject/data/targets.txt",array.type = "ProtoArray")

elist = my_data
library(limma)
my_data <- backgroundCorrect(elist, method="normexp", normexp.method="saddle")


plotArray(elist=my_data, idx=3, data.type="bg", log=FALSE, normalized=FALSE, aggregation="min", colpal="topo.colors")






elist = my_data
lot1 <- elist$targets[elist$targets$Batch=='Batch1','ArrayID']
lot2 <- elist$targets[elist$targets$Batch=='Batch1','ArrayID']
elist.bF <- batchFilter(elist=elist, lot1=lot1, lot2=lot2, log=FALSE, p.thresh=0.001, fold.thresh=3)

