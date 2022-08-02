
library(affy)

library(WGCNA)

library(GEOquery)

library(tidyverse)

library(limma)



GSE <- getGEO("GSE5281",GSEMatrix = T)

# read and write phenoData
pdata <-GSE$GSE5281_series_matrix.txt.gz@phenoData@data
write.table(pdata, "pdata_gse5281ec.txt", sep="\t", row.names = T)

# read and write featureData  
fdata <-GSE$GSE5281_series_matrix.txt.gz@featureData@data
write.table(fdata,"fdata_gse5281ec.txt", sep = "\t", row.names = T)

# read procced file and anno
data <- read.csv("GSE5281-EC-avg.csv",header = T, row.names = 1)

# annotation phenodata # not necessory 
temp <- pdata
pdata <- AnnotatedDataFrame(temp)
exp <- pdata@data 

# Annotate feature data
temp <- fdata
fdata <- AnnotatedDataFrame(fdata)

#it is may not necessory
phdata <- GSE$GSE5281_series_matrix.txt.gz@phenoData@data
write.table(phdata,"phdata_gse5281ec.txt", sep ="\t",row.names = T)  # it is not necessory

# see the dataset and identify the sample
view(phdata)
view(data)

# create a new data view the phdata and data 
phedata <- phdata[c(1:13,75:84),]
write.table(phedata,"phedata_gse5281ec.txt", sep ="\t", row.names = T)
view(phedata)

#Annotated phedata in pdata
temp <- phedata
pdata <- AnnotatedDataFrame(temp)


temp <- read.delim("fdata_gse5281ec.txt",header = T, row.names = 1)    # or temp <-fdata@data
view(temp)

# create a new data matrix and store annotated pdada(proccesd Annotated data) and fdata(Annotated featured data)
eset <- new("ExpressionSet", exprs=as.matrix(data))
eset@phenoData <- pdata
eset@featureData <- fdata

# find the affected and control indicated column note the affected sample in a vec
str <- phedata$title
affected <- grep("affected",str)
affected

# create total number of sample vec
sample <- c(1:nrow(phedata))
sample

#store 1 for affected sample and 0 for control sample
sam_aff <- sample%in%affected
sam_aff <- as.integer(sam_aff)
sam_aff

#create a factor to separate sample
gs <- factor(sam_aff)
gs

# give factor gropus names in levels
groups <- make.names(c("Control","Affected"))
levels(gs) <- groups

# put the factor in eset data
eset$group <- gs

#create design named matrix assign value 1 for activate and 0 for deactivate the affected and control matrix.
design <- model.matrix(~group +0, eset)
design

#give the design matrix colnames
colnames(design) <- levels(gs)
design

# fit the design into eset
fit <- lmFit(eset,design)
fit

# create cts from groups name and separate them "-" and make cont.matrix with help of cts and design
cts <- paste(groups[1],groups[2],sep="-")
``
cont.matrix <-makeContrasts(contrasts = cts, levels = design )
cont.matrix


fit2 <- contrasts.fit(fit,cont.matrix)
fit2


fit2 <- eBayes(fit2, 0.01 )
fit2


tT <- topTable(fit2, adjust="fdr", sort.by = "B", number=250)

tT <- subset(tT, select=c("ID","Gene.Title","Gene.Symbol","adj.P.Val","t","logFC"))

write.table(tT, file=stdout(),row.names = F,sep="\t")

tT2 <- topTable(fit2, adjust="fdr",sort.by = "B", number="inf")
tT2 <- subset(tT2, select=c("ID","Gene.Title","Gene.Symbol","adj.P.Val","t","logFC"))

write.table(tT2,"tT2_gse5281ec.txt",row.names = T,sep="\t")

#histgram plot tT2$adj.P.Val
hist(tT2$adj.P.Val, col="grey",border="white",xlab="P-Adj", ylab="Number of genes", main="P-Adj Value Distribution")

#vendiagram plot fit2 with help of adj.p.vlue
dT <-decideTests(fit2, adjust.method = "fdr", adj.p.value=0.05)

View(dT)

vennDiagram(dT, circle.col = "red")


t.good <- which(!is.na(fit2$F))
qqt(fit2$df.total[t.good], main="Moderated Statistic")

#valcano plot fit2
colnames(fit2)
ct <- 1
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch = 20, highlight = length(which(dT[,ct]!=0)),names = rep('+', nrow(fit2)))


