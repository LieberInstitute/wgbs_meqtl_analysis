###
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

library(MatrixEQTL)
library(sva)
library(minfi)
library(bsseq)

### full genome meQTLs
#load bsobj
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/DLPFC/PC1filtsmoothDLPFC.rda")
#load genotype info
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/DLPFC/Genotypes/Genotypes_processed.rda")

#put matrixes in same order
colnames(BSobj) <- colData(BSobj)$BrNum
mat <- match(colnames(snp),colnames(BSobj))
BSobj <- BSobj[,mat]

#extract phenotype data and M values
pd2 = colData(BSobj) 
p = as.matrix(getMeth(BSobj, type = "smooth"))

## snps
snp2 <- as.matrix(snp)
rownames(snpMap) = snpMap$SNP

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")

#only SNPs which have hg38 position data
snp2 <- snp2[!is.na(snpspos$pos),]
snpspos <- snpspos[!is.na(snpspos$pos),]

theSnps = SlicedData$new(snp2)
theSnps$ResliceCombined(sliceSize = 10000)

# model
mod = model.matrix(~pd2$MDS1 + pd2$MDS2 + 
                     pd2$MDS3 + pd2$MDS4 + pd2$MDS5)
#remove sex chromosomes, perform PCA:
oo = order(matrixStats::rowSds(p[1:27852739,]), decreasing=TRUE)[1:1000000]
pca = prcomp(t(p[oo,]))
nsv = num.sv(p[oo,], mod) # rawhippo: 64 smoothhippo:28 smoothdlpfc: 28 rawdlpfc:51
pca$x = pca$x[rownames(pca$x) %in% colnames(BSobj),]
covs = SlicedData$new(t(cbind(mod[,-1], pca$x[,1:nsv]))) 

## meth position
posMeth = data.frame(probeid = paste0(seqnames(BSobj), ".", start(BSobj)), chr = seqnames(BSobj), start = start(BSobj), end = end(BSobj))
rownames(p) <- posMeth$probeid

### make MatrixEQTL objects
meth = SlicedData$new(as.matrix(p))
meth$ResliceCombined(sliceSize = 5000)

rm(BSobj, p)
gc()

meMeth = Matrix_eQTL_main(snps=theSnps, gene = meth, 
                          cvrt = covs, output_file_name.cis =  "cis.txt" , output_file_name = "trans.txt",
                          pvOutputThreshold.cis = .01, pvOutputThreshold=0, #give this a value for trans
                          snpspos = snpspos, genepos = posMeth, 
                          useModel = modelLINEAR,	cisDist=20000)

# filter
meqtl = meMeth$cis$eqtl
meqtl$snps = as.character(meqtl$snps)
meqtl = meqtl[meqtl$FDR < 0.01,]
colnames(meqtl)[2] = "cpg"
meqtl$cpg = as.character(meqtl$cpg)

## annotate
m = match(meqtl$snps, snpMap$SNP)
meqtl$snpChr = snpMap$chr_hg38[m]
meqtl$snpPos = snpMap$pos_hg38[m]
meqtl$snpRsNum = snpMap$name[m]
meqtl$numImputed = snpMap$numImputed[m]
meqtl$snpCounted = snpMap$COUNTED[m]
meqtl$snpAlt = snpMap$ALT[m]


### see if SNP distrupts a CpG
library(BSgenome.Hsapiens.UCSC.hg38)
#meqtl$snpChr = gsub("chr23","chrX", meqtl$snpChr)
gr = GRanges(meqtl$snpChr, IRanges(meqtl$snpPos-1, meqtl$snpPos+1))
trio = getSeq(Hsapiens, gr)
meqtl$disruptCpG = vcountPattern("CG", trio)

meqtl$methChr = as.character(posMeth$chr[match(meqtl$cpg, posMeth$probeid)])
meqtl$methPos = posMeth$start[match(meqtl$cpg, posMeth$probeid)]
meqtl$distMethToSnp = meqtl$methPos - meqtl$snpPos

###PGC Analysis
#filter SNPs to only selected PGC SNPs
#meQTL analysis input
meMeth = Matrix_eQTL_main(snps=theSnps, gene = meth, 
                          cvrt = covs, output_file_name.cis =  "cis.txt" , output_file_name = "trans.txt",
                          pvOutputThreshold.cis = 1, pvOutputThreshold=1, #give this a value for trans
                          snpspos = snpspos, genepos = posMeth, 
                          useModel = modelLINEAR,	cisDist=250000)





