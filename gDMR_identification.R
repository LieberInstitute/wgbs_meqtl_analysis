library(bumphunter)
library('bsseq')
library('bumphunter')
library('devtools')
library("rtracklayer")
library('limma')

#load all cis meQTL stats
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/Hippocampus/meQTL/FULLpgccis.rda")
cismeqtl = meqtl
#load posMeth
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/Hippocampus/posMeth.rda")
#load PGC snps
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/Hippocampus/meQTL/pgcsnps.rda")

#load trans meQTLs for each chromosome, order all meQTL statistics -- this is only so laborious due to a bug in regionFinder code
`%notin%` <- Negate(`%in%`)
full = c()
for (i in rownames(snpMap)){
  subcis = cismeqtl[cismeqtl$snps == i,]
  chr = snpMap[i,]$CHR
  load(paste0("/dcl01/lieber/ajaffe/Kira/SczWGBS/Hippocampus/meQTL/FULLpgc/FULLpgctrans_chr", chr, ".rda"))
  subtrans = transmeqtl[transmeqtl$snps == i,]
  colnames(subtrans)[2] = "cpg"
  subtrans$snpChr = NA
  subtrans$snpPos = NA
  subtrans$snpRsNum = NA
  subtrans$snpCounted = NA
  subtrans$snpAlt = NA
  subtrans$disruptCpG = NA
  subtrans$methChr = as.character(posMeth$chr[match(subtrans$cpg, posMeth$probeid)])
  subtrans$methPos = posMeth$start[match(subtrans$cpg, posMeth$probeid)]
  subtrans$distMethToSnp = NA
  subtrans = subtrans[subtrans$methPos %notin% subcis$methPos,]
  sub = rbind(subcis, subtrans)
  oo = order(sub$methPos)
  sub = sub[oo,]
  #perform region finding to identify DMRs
  regions = regionFinder(abs(sub$statistic), sub$methChr, sub$methPos, cutoff = 5)
  if (nrow(regions) > 0){
    regions$snp = i
    full = rbind(full, regions)
  }
  rm(transmeqtl)
  gc()
}

regions = full

#annotate all DMRs
regions$width = regions$end - regions$start
m = match(regions$snp, rownames(snpMap))
regions$snpPos = snpMap$pos_hg38[m]
regions$disruptCpG = snpMap$disruptCpG[m]
regions$distStart = regions$snpPos - regions$start
regions$distEnd = regions$snpPos - regions$end

#####gene annotation and ontology
gr = makeGRangesFromDataFrame(regions, keep.extra.columns=TRUE)

#gene annotation
gencode_v29_hg38 <- import('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz')
genes <- gencode_v29_hg38
genes <- genes[genes$type == "gene",]

annotations <- annotateNearest(gr, genes, annotate = TRUE)
annotations$gene_name <- genes$gene_name[annotations$subjectHits]
annotations$gene_id <- genes$ID[annotations$subjectHits]

annotations <- annotations[abs(annotations$distance) < 10000,]
uniquelist <- unique(annotations$gene_name)
uniqueids <- unique(annotations$gene_id)


#gene ontology
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
universe <- gencode_v29_hg38
universe <- universe[universe$type == "gene",]
universe <- universe$gene_id
universe <- jaffelab::ss(universe, "\\.")
gene <- uniqueids
gene <- jaffelab::ss(gene, "\\.")

ego <- enrichGO(gene          = gene,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)




