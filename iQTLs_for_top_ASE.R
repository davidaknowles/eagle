# genotype: a SNPs x individuals matrix of {0,1,2}
# ge: gene expression data, individuals x genes

max_genomic_distance=1e5 # look with 100kb of TSS

# quantile normalization 
geNorm=apply(ge,2,function(g) qqnorm(g,plot.it = F)$x)
dimnames(geNorm)=dimnames(ge)

# load TSSes
tss=read.table("tss.txt.gz", header=T, stringsAsFactors=F)

# aseGeneName: which gene to test (i.e. what gene is the locus in). Multiple genes should be separated by commas
# env: environment factor
mostCommonEnv=as.numeric(names(which.max(table(env))))
aseGeneName=geneNames[locus]

genesToTest=intersect( strsplit(aseGeneName,",")[[1]], colnames(ge))
genesToTest=intersect( genesToTest, as.character(tss$geneName))

pvals=list()
cisSnpList=list()

# cycle over genes (typically only one)
for (gni in seq_len(length(genesToTest))){
  gn=genesToTest[gni]
  pvals[[ gn ]] = c()
  ind=which(tss$geneName==gn)[1]
  chr=tss$chr[ ind ]
  gene_tss=tss$tss[ ind ]
  cisSnps=which(snpInfo$chrom==chr & (snpInfo$pos > 0) & (abs(snpInfo$pos-gene_tss) < max_genomic_distance)
  cisSnpList[[ gn ]]=cisSnps
  cat("num cis SNPs:",length(cisSnps),"\n")
  
  # cycle over cisSNPs (currently all on chromosome)
  for (cisSnpI in seq_len(length(cisSnps))){
    p=NA
    cisSnp=cisSnps[cisSnpI]
    if (cisSnpI %% 1000 == 0) cat("snp ",cisSnpI,"\n")
    snpData=as.numeric(t(genotype[cisSnp,]))
    expr=geNorm[,gn]
    
    # if environment factor is binary remove categories with <5 samples
    if (length(unique(env))<=2){
      toosmall=table(data.frame(snpData,env)) < 5
      for (wts in which(toosmall)){
        ai=arrayInd(wts,dim(toosmall))
        expr[ snpData==rownames(toosmall)[ai[1] ] & env==colnames(toosmall)[ai[2] ] ]=NA
      }
    }
    mostCommonGenotype=as.numeric(names(which.max(table(snpData))))
    temp=table(snpData==mostCommonGenotype,env==mostCommonEnv)
    if (all(temp>10)){ # check test will be well conditioned
      full=lm( expr ~ snpData*env) # fit null model
      null=lm( expr ~ snpData + env) # fit alternative model
      p=anova(full,null)[["Pr(>F)"]][2]
    }
    pvals[[ gn ]]=c(pvals[[gn]],p)
  }
}
