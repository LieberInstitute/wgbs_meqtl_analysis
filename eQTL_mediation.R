library('bsseq')
library('bumphunter')
library('devtools')
library("rtracklayer")
library('limma')

#load and arrange everything
#expression
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n453.rda")
#methylation
load("PC1filtsmoothDLPFC.rda")
#snps
load("meQTL/pgcsnps.rda")
#eqtls
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/DLPFC/PGCeqtls.rda")
#meqtls
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/DLPFC/meQTL/annotatedPGCmeqtl.rda")
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/Hippocampus/posMeth.rda")

#rename and rearrange
colnames(rse_gene) <- colData(rse_gene)$BrNum
rse_gene <- rse_gene[,colnames(rse_gene) %in% colnames(snp)]
mat <- match(colnames(rse_gene), colnames(snp))
snp <- snp[,mat]
mat <- match(colnames(rse_gene), colnames(BSobj))
BSobj <- BSobj[,mat]

explevel <- log2(recount::getRPKM(rse_gene, 'Length')+1)

###find meQTL-eQTL pairs to interrogate that have correlation between M and exp
fullhits <- data.frame(row = as.numeric(), col = as.numeric(), gene = as.character(), snp = as.character(), cpg = as.character())
for (i in 1:nrow(snp)){
  testsnp = rownames(snp)[i]
  if (testsnp %in% eqtl$snps & testsnp %in% cismeqtl$snps) {
    testeqtl <- eqtl[eqtl$snps == testsnp,]
    meqtl <- cismeqtl[cismeqtl$snps == testsnp,]
    gr <- makeGRangesFromDataFrame(meqtl, seqnames.field = 'methChr', start.field = 'methPos', end.field = 'methPos', keep.extra.columns=TRUE)
    oo <- findOverlaps(gr, BSobj)
    BS <- BSobj[subjectHits(oo),]
    meth <- assays(BS)$M
    sig <- testeqtl$FDR < 0.05
    exp <- explevel[rownames(rse_gene) %in% testeqtl$gene[sig],]
    ##if only one sig eqtl
    if (sum(rownames(rse_gene) %in% testeqtl$gene[sig]) == 1){
      exp <- as.data.frame(exp)
      colnames(exp) = rownames(rse_gene)[rownames(rse_gene) %in% testeqtl$gene[sig]]
      cor <- cor(t(meth), exp)
      rownames(cor) = meqtl$cpg
      if (sum(abs(cor) > 0.3) == 0){
        hits <- data.frame(row = as.numeric(), col = as.numeric(), gene = as.character(), snp = as.character(), cpg = as.character())
      } else {
        hits <- which(abs(cor) > 0.3, arr.ind = TRUE)
        hits <- as.data.frame(hits)
        hits$gene = colnames(cor)[hits$col]
        hits$snp = testsnp
        hits$cpg = rownames(cor)[hits$row]
      }
      ##if no sig eqtls
    } else if (sum(rownames(rse_gene) %in% testeqtl$gene[sig]) == 0) {
      hits <- data.frame(row = as.numeric(), col = as.numeric(), gene = as.character(), snp = as.character(), cpg = as.character())
    } else {
      cor <- cor(t(meth), t(exp))
      rownames(cor) = meqtl$cpg
      if (sum(abs(cor) > 0.3) == 0){
        hits <- data.frame(row = as.numeric(), col = as.numeric(), gene = as.character(), snp = as.character(), cpg = as.character())
      } else {
        hits <- which(abs(cor) > 0.3, arr.ind = TRUE)
        hits <- as.data.frame(hits)
        hits$gene = colnames(cor)[hits$col]
        hits$snp = testsnp
        hits$cpg = rownames(cor)[hits$row]
      }
    }
    fullhits <- rbind(fullhits,hits)
  }
}


###coefficient for expression~eQTL
fullco <- data.frame(Estimate = as.numeric(), 'Std. Error'= as.numeric(), `t value` = as.numeric(), `Pr(>|t|)` = as.numeric())
for (i in 1:nrow(fullhits)){
  testsnp <- fullhits$snp[i]
  genotype <- as.numeric(snp[rownames(snp) == testsnp,])
  mod = model.matrix(~genotype)
  exp = explevel[rownames(rse_gene) == fullhits$gene[i],]
  exp <- exp[as.numeric(rownames(mod))]
  genotype <- genotype[!is.na(genotype)]
  lm = lm(exp ~ genotype)
  coefs <- summary(lm)$coef[2,]
  fullco <- rbind(fullco, coefs)
}
snpmodel <- fullco

###add methylation variable
for (i in 1:nrow(fullhits)){
  testsnp <- fullhits$snp[i]
  genotype <- as.numeric(snp[rownames(snp) == testsnp,])
  mod = model.matrix(~genotype)
  exp = explevel[rownames(rse_gene) == fullhits$gene[i],]
  exp <- exp[as.numeric(rownames(mod))]
  genotype <- genotype[!is.na(genotype)]
  pos = posMeth[posMeth$probeid == fullhits$cpg[i],]
  gr = GRanges(pos$chr, IRanges(pos$start, pos$end))
  oo = findOverlaps(gr, BSobj)
  BS = BSobj[subjectHits(oo),]
  BS <- BS[,as.numeric(rownames(mod))]
  meth = assays(BS)$M[1,]
  newlm = lm(exp ~ genotype + meth)
  coefs <- summary(newlm)$coef[2,]
  fullco <- rbind(fullco, coefs)
}
duomodel <- fullco

###compare "Estimate" values
