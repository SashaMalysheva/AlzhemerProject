
library(PAA)
my_data <- loadGPR(
  gpr.path="D:/MVS/Alzheimer/DATA_RAW",
  targets.path = 
    "D:/MVS/Alzheimer/DATA_RAW/targets.txt",
  array.type = "ProtoArray")

my_data1 <- loadGPR(
  gpr.path="D:/MVS/Alzheimer/DATA_RAW",
  targets.path = 
    "D:/MVS/Alzheimer/DATA_RAW/targets1.txt",
  array.type = "ProtoArray")


plotArray(elist=my_data,  idx=3, data.type="bg", log=FALSE, normalized=FALSE,aggregation="min", colpal="topo.colors")

plotArray(elist=elist1,  idx=3, data.type="bg", log=FALSE, normalized=FALSE,aggregation="min", colpal="topo.colors")


#Filtring
library(limma)
elist1 <- backgroundCorrect(my_data1, method="normexp",
                           normexp.method="saddle")

#Normalization

norm_data1 <- normalizeArrays(elist = elist1, 
                             method = "cyclicloess", cyclicloess.method = "fast")


##Грузим таблицу с нормализованными данными
load('normalized.Rdata')

##Извлекаем названия белков, пихаем в таблицу данных

gn = sub(pattern = 'Hs~', rep = "", x = norm_data$genes$Name)
gn = sub(pattern = '~.*', rep = "", gn)
gn = sub(pattern = '.*:', rep = "", gn)
colnames(df) <- gn

##Извлекаем матрицу, аггрегируем, транспонируем,сохраняем в виде таблицы данных

aggr1 <- aggregate(norm_data1$E, list(gn), mean)
df2 <- as.data.frame(t(aggr1[,-1]))
names(df2) <- aggr1[,1]

#Сохраняем/грузим данные

write.csv(df, file='Complete_table.csv')

##Смотрим, что различается между группами.
##Нормальности распределения не увидел в некоторых случаях,
##скодировать массовую проверку/бутстреп не умею, так что будет Вилкоксон

rr <- apply(df,2,function(x) {t.test(x[1:50],x[51:90])$p.value})

rr <- apply(df[,-9354],2,function(x) {wilcox.test(x[1:50],x[51:90])$p.value})

rr <- apply(df[,-9354],2,function(x) {oneway_test(x[1:50],x[51:90])$p.value})

rr <- apply(df2[,-9354],2,function(x){kruskal.test(x~groups)$p.value})

rr <- apply(df[,-9354],2,function(x) {pvalue(oneway_test(x~groups))})

rr0 <- apply(df[,-9354],2,function(x) {pvalue(kruskal_test(x~groups))})

oneway_test(df$AB065619.1~df$Group)

groups <- factor(c(rep('AD', 50), rep('NDC', 40)))

rr1 <- p.adjust(rr, method='bonferroni', n=length(rr))
rr2 <- rr1[rr1<10^-5]
rr3 <- rr1[rr1<10^-6]



##Оставляем в таблице только достоверно отличающиеся между группами белки

df1 <- df[,names(rr2)]
df_ad <- df[,names(rr3)]
df2 <- df[,names(rr3)]

##Очистка от коррелирующих
library(caret)
descrCor <-  cor(df1)
#highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .75)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)

highlyCorDescr <- findCorrelation(descrCor, cutoff = .75, exact=T)

df1 <- df1[,-highlyCorDescr]
#descrCor2 <- cor(df2)
#summary(descrCor2[upper.tri(descrCor2)])
##Corrplot for visualization
library(corrplot)
corrplot(descrCor, method='pie', tl.pos='n', order='hclust')


##Анализ главных компонент.

my_pca=prcomp(df, scale=T, center=T)
imp=summary(my_pca)$importance

##Вносим принадлежность к группе. График каменистой осыпи



barplot(imp[2,]*100,ylab="Persentage of variance", 
        xlab="Principal Components",
        main = "Variance explained by individual PC", col="cyan3")

##Делаем таблицу данных из объекта - иначе его ggplot2 не распознает

my_pca1 <- as.data.frame(my_pca$x)
my_pca1$Group <- factor(c(rep('AD', 50), rep('NDC', 40)))

##График главных компонент
library(ggplot2)
ggplot(my_pca1, aes(x = PC1, y = PC2, colour = Group)) + 
  theme_bw(base_size = 8) + 
  geom_point()+
  labs(x = 'PC1', y = 'PC2') + 
  scale_colour_manual(values = c("#d4170a","#00ba38")) + 
  theme(
    legend.position = c(1, 0),
    legend.direction = "vertical",
    legend.justification = c(1, 0)
  )+
  stat_ellipse()

##3D главные компоненты

library(pca3d)
p <- pca3d(my_pca, group=my_pca1$Group)

##Подбор моделей. Черновой вариант

df2$Group <- factor(c(rep('PD',20),rep('AD', 50), rep('NDC', 40)))

fitControl <- trainControl(## 7-fold CV
  method = "repeatedcv",
  number = 8,
  repeats = 15,
  classProbs = TRUE,
  savePredictions = T,
  summaryFunction = twoClassSummary)

confusionMatrix( data = glmnet_2$pred[,1], 
                        reference = glmnet_2$pred[,2])

##Random forest для начала
set.seed(900)
set.seed(1000)
inTraining <- createDataPartition(df1$Group, p = .75, list = FALSE)
training <- df1[ inTraining,]
testing  <- df1[-inTraining,]

df$Group <- groups

rfmodel_1 <- train(Group ~.,data=training, method='glm', trControl=fitControl)

rfmodel_2 <- train(Group~., data=training, method='rf', trControl=fitControl)
rfmodel_1$results
predict(rfmodel_1, newdata = testing, type='prob')

predictions <- predict(rfmodel_1, testing)
confusionMatrix(predictions, testing$Group)

importance <- varImp(rfmodel_1, scale=FALSE)
importance


##Разрыв мозга с glmnet

glmnet_1 <- train(Group ~.,data=training, method='glmnet', trControl=fitControl, family='multinomial')
glmnet_1$results
t <- predict(glmnet_1, newdata=testing, type='prob')

glmnet_2 <- train(Group ~.,data=training, method='glmnet', trControl=fitControl)

glmnet_3 <- train(Group ~.,data=training[,c(1:10, 26)], method='glmnet', trControl=fitControl)

glmnet_4 <- train(Group~., data=training[,c(6,5,7,26)], method='glmnet', trControl=fitControl)

predictions <- predict(glmnet_1, testing)
confusionMatrix(predictions, testing$Group)

predictors(glmnet_3)


glmodel_1 <- train(Group ~.,data=training, method='glm', family='binomial',trControl=fitControl)

predictions=predict(glmodel_1, testing)
confusionMatrix(predictions, testing$Group)


## А теперь визуализируем

importance <- varImp(glmnet_1, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

names(df2)

##Pheat map?

library(pheatmap)
heatmap(as.matrix(df1[,-9354]))


pheatmap(df1, scale='row', cluster_rows=F,
         cluster_cols=F, show_rownames = F,show_colnames = F)
pheatmap(t(df1), scale='row',cluster_rows=T, 
         cluster_cols=F, show_rownames = T,show_colnames = T)
