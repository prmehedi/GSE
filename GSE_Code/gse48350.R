

library(affy)

library(WGCNA)

library(GEOquery)

library(tidyverse)

library(limma)

gse <- getGEO("GSE48350", GSEMatrix = T)

pdata <- gse$GSE48350_series_matrix.txt.gz@phenoData@data

write.table(pdata, "pdata_gse4835.txt", sep="\t",row.names = T)

fdata <- gse$GSE48350_series_matrix.txt.gz@featureData@data

write.table(fdata,"fdata_gse4835.txt", sep="\t",row.names = T)


data <- read.delim("GSE48350_HIPP_avg.txt",header = T, row.names = 1)

phdata <-read.delim("phdata.txt",header = T, row.names = 1)

tem <- phdata

phdata <- AnnotatedDataFrame(tem)

tem <- fdata

fdata <- AnnotatedDataFrame(tem)



eset <- new("ExpressionSet", exprs = as.matrix(data))

eset@phenoData <- phdata

eset@featureData <- fdata

phdata
str <- phdata[,1]

row <- grep("AD", str)
row

gsms <- c(1:nrow(phdata))
gsms <-gsms%in%row
gsms <- as.integer(gsms)
gsms
gs <- factor(gsms)

groups <- make.names(c("CN","HD"))
levels(gs) <- groups
eset$group <- gs

design <- model.matrix(~group +0,eset)
colnames(design) <- levels(gs)

fit <- lmFit(eset,design)

cts <- paste(groups[1],groups[2], sep ="-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)

fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2,0.01)

tT <- topTable(fit2,adjust ="fdr", sort.by = "B", number = 250)

tT <- subset(tT, select=c("ID","Gene.Title","Gene.Symbol","adj.P.VAl","t","logFC"))




tT2 <- topTable(fit2,adjust="fdr",sort.by = "B",number = "inf")
write.table(tT2,"tT2.txt", sep="\t",row.names = T)

hist(tT2$adj.P.Val, col = "grey",border="white",xlab = "P-adj", ylab="Number of Genes", main="P-adj Value Distribution")


dT <- decideTests(fit2, adjust.method = "fdr", adj.p.val=0.05 )

vennDiagram(dT, circle.col = palette())


t.good <- which(!is.na(fit2$F))

qqt(fit2$t[t.good], fit2$df.total[t.good],main="Moderated t statistic")


colnames(fit2)
ct <- 1
volcanoplot(fit2, coef=ct, main= colnames(fit2)[ct],pch=20,highlight = length(which(dT[,ct]!=0)),names = rep('+',nrow(fit2)))





