

library(affy)

library(WGCNA)

library(GEOquery)

library(tidyverse)

library(limma)

setwd("C:/Users/acer/Desktop/mehedi/GSE/GSE3790_CN")

gse <- getGEO("GSE3790", GSEMatrix = T)

pdata <- gse$`GSE3790-GPL96_series_matrix.txt.gz`@phenoData@data

write.table(pdata, "pdata_gse3790.txt", sep="\t",row.names = T)

fdata <- gse$`GSE3790-GPL96_series_matrix.txt.gz`@featureData@data

write.table(fdata,"fdata_gse3790.txt", sep="\t",row.names = T)


data <- read.delim("GSE3790_CN_avg.txt",header = T, row.names = 1)

tmp <- read.delim("pdata_gse3790.txt",header = T)

pdata <- AnnotatedDataFrame(tmp)

exp <- pdata@data


phdata <- gse$`GSE3790-GPL96_series_matrix.txt.gz`@phenoData@data

write.table(phdata,"phdata.txt", sep="\t", row.names = T) #it is necessory to write

phedata <- phdata[c(1:70),]

write.table(phedata,"phedata.txt",sep="\t", row.names = T)

tmp <- phedata

pdata <- AnnotatedDataFrame(tmp)

tmp <- fdata

fdata <- AnnotatedDataFrame(tmp)

eset <- new("ExpressionSet", exprs = as.matrix(data))

eset@phenoData <- pdata

eset@featureData <- fdata


str <- phedata$characteristics_ch1
row <- grep("HD", str)
row

gsms <- c(1:nrow(phedata))

gsms <-gsms%in%row
gsms <- as.integer(gsms)

gs <- factor(gsms)

groups <- make.names(c("Con","HD"))

levels(gs) <- groups

eset$group <- gs

design <- model.matrix(~group +0,eset)
View(design)

colnames(design) <- levels(gs)

fit <- lmFit(eset,design)

cts <- paste(groups[1],groups[2], sep ="-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)

fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2,0.01)

tT <- topTable(fit2,adjust ="fdr", sort.by = "B", number = 250)
tT <- subset(tT, select=c("ID","Gene.Title","Gene.Symbol","adj.P.VAl","t","logFC"))


tT2 <- topTable(fit2,adjust="fdr",sort.by = "B",number = "inf")
write.table(tT2,"tT3col.txt", sep="\t",row.names = T)


hist(tT2$adj.P.Val, col = "grey",border="white",xlab = "P-adj", ylab="Number of Genes", main="P-adj Value Distribution")

dT <- decideTests(fit2, adjust.method = "fdr", adj.p.val=0.05 )
vennDiagram(dT, circle.col = palette())

t.good <- which(!is.na(fit2$F))
qqt(fit2$t[t.good], fit2$df.total[t.good],main="Moderated Statistic")

colnames(fit2)
ct <- 1
volcanoplot(fit2, coef=ct, main= colnames(fit2)[ct],pch=20,highlight = length(which(dT[,ct]!=0)),names = rep('+',nrow(fit2)))
write.table(tT2,"tT2_gse3790cn.txt", sep="\t",row.names = T)

