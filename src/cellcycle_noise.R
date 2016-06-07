library(statmod)
library(gdata)
library(genefilter)
library(EBImage)
library(rhdf5)
library(DESeq)
library(statmod)
library(hom.Hs.inp.db)
library(AnnotationDbi)
library(org.Hs.eg.db)


# Load data
data <- read.table("data/runs/BSF_0222_HEK293T_Cas9/digital_expression.100genes.assigned.tsv", sep="\t", row.names=1, header=TRUE)

#log-transformed counts 
LData <- log10(data + 1)
LmeansData <- rowMeans(LData)
LvarsData <- rowVars(LData)
Lcv2Data <- LvarsData / LmeansData ^ 2

#Fit baseline variability
LogNcountsList = list()
useForFitL = LmeansData > 0.3
LogNcountsList$mean = LmeansData[useForFitL]
LogNcountsList$cv2 = Lcv2Data[useForFitL]
fit_loglin = nls(cv2 ~ a * 10^(-k*mean), LogNcountsList, start=c(a=10, k=2))

# TECHNICAL noise
#variable genes
is_het = (coefficients(fit_loglin)["a"] * 10 ^ (-coefficients(fit_loglin)["k"] * LmeansData) < Lcv2Data) & LmeansData > 0.3
tech_noise <- coefficients(fit_loglin)["a"] * 10 ^ (-coefficients(fit_loglin)["k"] * LmeansData) * LmeansData ^ 2

#plot mean/cv2 relationship and variable genes
plot(LmeansData, Lcv2Data, log="y", col=1 + is_het, xlab='meansLogData', ylab='cv2LogData', ylim=c(1e-3,1e2))  
xg <- seq(0, 4.5, length.out=100)
lines(xg, coefficients(fit_loglin)["a"] * 10 ^ (-coefficients(fit_loglin)["k"]*xg),lwd=2,col='blue')
legend('topright',c('Variable genes'),pch=c(1),col=c('red'),cex=0.8)

#gene names in the Data-Seq mESC Data.
gene_names = rownames(data)
gene_names_het = gene_names[is_het]

#get cell cycle genes from GO 
xxGO <- as.list(org.Hs.egGO2EG)
cell_cycleEG <-unlist(xxGO['GO:0007049'])
#get ENSEMBLE ids
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xxE <- as.list(x[mapped_genes])
cc_genes<-unlist(xxE[cell_cycleEG])

# Get cell cycle genes from Cyclebase
# wget http://download.jensenlab.org/Cyclebase3/human_periodic.tar
# tar xfz human_periodic.tar

per = read.table("human_periodic.tsv", sep="\t", header=TRUE)
# map ensembl protein to index
x = org.Hs.egENSEMBLPROT2EG
xxE <- as.list(x[unlist(per["gene"])])
# map index to symbol
x = org.Hs.egSYMBOL
cc_genes_cyclebase <- unlist(as.list(x[unlist(xxE)]))

# Join both (GO and Cyclebase) together
cc_genes = unique(c(cc_genes_cyclebase, cc_genes))

h5save(cc_genes, tech_noise, gene_names_het, file='cellcycle_noise.h5f')
