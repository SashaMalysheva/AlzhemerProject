#Load GPR files
library(PAA)
my_data <- loadGPR(
  gpr.path="C:/Users/alex/alzhemerproject/Alzhemerproject/data",
  targets.path = 
    "C:/Users/alex/alzhemerproject/Alzhemerproject/data/targets.txt",
  array.type = "ProtoArray");


#plotArray(elist=my_data,  idx=3, data.type="bg", log=FALSE, normalized=FALSE,aggregation="min", colpal="topo.colors")


#Filtring
library(limma)
elist <- backgroundCorrect(my_data, method="normexp",
                           normexp.method="saddle");


#Normalisation
norm_data <- normalizeArrays(elist = elist, 
                             method = "cyclicloess", cyclicloess.method = "fast");


#Clear names
name <- norm_data$genes$Name;
gn = sub(pattern = 'Hs~', rep = "", x = name)
gn = sub(pattern = '~.*', rep = "", gn)
gn = sub(pattern = '.*:', rep = "", gn)


#Naming
row.names(norm_data$E) <- gn;


#Meaning
E <- norm_data$E
mean_data <- aggregate(E, list(gn), mean);
row.names(mean_data) <- mean_data$Group.1
mean_data <- mean_data[-1];


#P-value
p_premult <- apply(mean_data,1,function(x) {t.test(x[1:50],x[51:90])$p.value});


#false discovery rate multy normalisation
adjusted <- p.adjust(p_premult, method = "fdr", n = length(p_premult));


#Just low p-value
adj <- sort(adjusted[adjusted < 10^-5]);


#Choise samples vith low p-value
re <- which(rownames(mean_data) %in% names(adj));
er <- as.matrix(mean_data[re, ]);


#hitmap for samples with minimum p-value
library(pheatmap)
heatmap(er, col = cm.colors(256))
er<-t(er)


#Убираем кариляцию
library(caret)
descrCor <-  cor(er)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
er <- er[,-highlyCorDescr]


##Анализ главных компонент.
my_pca=prcomp(er, scale=T, center=T)
imp=summary(my_pca)$importance


##Вносим принадлежность к группе. График каменистой осыпи
barplot(imp[2,]*100,ylab="Persentage of variance", 
        xlab="Principal Components",
        main = "Variance explained by individual PC", col="cyan3")


##Делаем таблицу данных из объекта - иначе его ggplot2 не распознает
my_pca1 <- as.data.frame(my_pca$x)
my_pca1$Group <- NA
my_pca1[1:50,'Group'] <- 'AD'
my_pca1[51:90,'Group'] <- 'NDC'


##График главных компонент
library(ggplot2)
ggplot(my_pca1, aes(x = PC1, y = PC2, colour = Group)) + 
  theme_bw(base_size = 8) + 
  geom_point()+
  labs(x = 'PC1', y = 'PC2') + 
  scale_colour_manual(values = c("#00ba38","#d4170a","#00ebf2")) + 
  theme(
    legend.position = c(1, 0),
    legend.direction = "vertical",
    legend.justification = c(1, 0)
  )+
  stat_ellipse()

