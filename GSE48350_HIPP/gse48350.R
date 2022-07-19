# get the pheno and the feature data and write it down in a table because we have to create an expression set of the csv files before we can proceed to do a limma test




library(affy)

library(GEOquery)

library(WGCNA)

library(tidyverse)

library(limma)

#library(limma)

gse <- getGEO("GSE48350", GSEMatrix = T)

pdata <- gse$GSE48350_series_matrix.txt.gz@phenoData@data

write.table(pdata, "pdata.txt", quote = F, sep = "\t", row.names = T)

fdata <- gse$GSE48350_series_matrix.txt.gz@featureData@data

write.table(fdata, "fdata.txt", quote = F, sep = "\t", row.names = T)

# load the expression data, now this data is normalised and processed
# "GSE48350_HIPP_avg.txt" download from "https://drive.google.com/drive/folders/1KYfyOmwdpF26taYBxG7M6WPMuUmXm9Rh?usp=sharing"
data <- read.delim("GSE48350_HIPP_avg.txt", header = T, row.names = 1)

tmp <- read.delim("pdata.txt", header = T)

pdata <- AnnotatedDataFrame(tmp)

exp <- pdata@data

phdata <- gse$GSE48350_series_matrix.txt.gz@phenoData@data

write.table(phdata, "phdata.txt", quote = F, sep = "\t", row.names = T)

tmp <- read.delim("phedata.txt", header = T, row.names = 1)

pdata <- AnnotatedDataFrame(tmp)


# do the same with the feature data

tmp <- read.delim("fdata.txt", header = T, row.names = 1)

fdata <- AnnotatedDataFrame(tmp)


# now make an expression set using biobase generic functions

eset <- new("ExpressionSet", exprs=as.matrix(data))

eset@phenoData <- pdata

eset@featureData <- fdata

gsms <- paste0("00000000000000000001111111111111111111111111")

sml <- strsplit(gsms, split = "")[[1]]

# assign samples to groups and set up design matrix
gs <- factor(sml)

groups <- make.names(c("AD","CN"))

levels(gs) <- groups


eset$group <- gs

design <- model.matrix(~group + 0, eset)

colnames(design) <- levels(gs)

fit <- lmFit(eset, design)  # fit linear model


# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")

cont.matrix <- makeContrasts(contrasts=cts, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.Title"))

write.table(tT, file=stdout(), row.names=F, sep="\t")

tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", adj.p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes

qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")


# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))


