# Compute signature scores for EMT6 RNAseq samples
library(dplyr)
library(edgeR)
############################################################
# Get logCPM values
############################################################
basedir <- "~/workspace/TGFbeta/EA17073/compute_scores/"
dge <- readRDS(paste0(basedir,"EMT6_dge_list.rds"))
cpmd <- cpm(dge)
# Filter out genes not having at least 1 CPM in at least N samples
# N is the number of samples in each treatment group
N <- 8
keep.exprs <- rowSums(cpmd > 1) >= N
d <- dge[keep.exprs,,keep.lib.sizes=TRUE]
d <- calcNormFactors(d)
dat <- cpm(d, log=TRUE,prior.count=1)
pheno <- d$samples
saveRDS(pheno,paste0(basedir,"pheno.rds"))
genes <- d$genes
############################################################
# Read in signatures
############################################################
allsigs <- read.csv(paste0(basedir,"Figure4signatures.csv"),header=T)
allsigs <- mutate(allsigs,combinedName = paste(Collection,Name,sep="."))
############################################################
# Scores will be computed using the mean Z approach 
# (the score for a sample is defined as the average of gene-specific z-scores for the genes
# in the signature)
# However, the Calon Ma-TBRS and T-TBRS have a large number of genes (1284 and 71 respectively,
# as opposed to 19 and 8 for the Pan-F-TBRS and CD8 Teff signatures).
# Therefore the Calon Ma-TBRS and T-TBRS will first be reduced by computing the 
# first principal component loadings, and retaining only genes with at least 0.8 Pearson
# correlation with these loadings.

small.sig <- allsigs[grep("Calon",allsigs$combinedName),]
tmp <- split(as.character(small.sig$GeneSymbol),small.sig$combinedName)
calonRes  <- vector("list",2)

THRESH <- 0.8
for(i in 1:length(tmp)){
  symbols <- unique(sort(tmp[[i]]))
  # GET MOUSE ORTHOLOGS for genes in the Calon signatures
  conversionDB <- human2mouse(symbols)
  mx <- match(conversionDB$MGI.symbol,genes$symbol)
  cdb2 <- cbind(conversionDB,res = genes$symbol[mx])
  
  dat.nonna <- dat[mx[!is.na(mx)],]
  pcx <- prcomp(t(dat.nonna))
  # correlation with PC1
  corx <- cor(pcx$x[,1],t(dat.nonna))
  colnames(corx) <- genes$symbol[mx[!is.na(mx)]]
  abs.corx <- abs(corx)
  # Get top genes
  top.mouse <- unique(colnames(abs.corx)[which(abs.corx > THRESH)])
  calonRes[[i]] <- cdb2$HGNC.symbol[match(top.mouse,cdb2$res)]
}#end for i
names(calonRes) <- names(tmp)
# 123 and 18 genes remained for the Ma-TBRS and T-TBRS signatures respectively

# NOW GET MOUSE ORTHOLOGS for genes in the other signatures
newsig <- filter(allsigs,!grepl("Calon",combinedName))
glist <- split(as.character(newsig$GeneSymbol),newsig$combinedName)
glist <- c(glist,calonRes)
symbols <- unique(sort(as.character(allsigs$GeneSymbol)))
conversionDB <- human2mouse(symbols) 

gprofile <- genes$symbol
scores <- score.signatures.DB.mouse(x=dat,glist=glist,gprofile=gprofile,conversionDB=conversionDB)
#saveRDS(scores,paste0(basedir,"scores.rds"))

