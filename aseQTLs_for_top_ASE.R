# example script to use EAGLE to look for candidate variants for specific environment - exonic locus pair. 

# env: environment
# chr: what chromosome the ASE locus is on
# gn: what gene the ASE locus is in. Seperated by "," if multiple. 
tss=read.table("tss.txt.gz", header=T, stringsAsFactors=F)
# genotype: # cisSnps x # individuals, 0,1,2 coding
# snpInfo: data.frame with chromosome and position of the SNPs in "genotype"
# hets, alt, ref: # individuals x # ASE loci
# locus: current locus to test

mostCommonEnv=as.numeric(names(which.max(table(env))))
tssid=which( as.character(tss$geneName) %in% strsplit(gn,",")[[1]] )
tssForThisLocus=if (length(tssid)>0) tss[ tssid[1] , "tss"] else locusinfo[locus,"pos"] # if can't find tss just use locus itself
cisSnps=which(snpInfo$chrom==chr & (snpInfo$pos > 0) & (abs(snpInfo$pos-tssForThisLocus) < 3e5))
pvals = numeric(length(cisSnps))*NA
cat("num cis SNPs:",length(cisSnps),"\n")
for (cisSnpI in seq_len(length(cisSnps))){
  cisSnp=cisSnps[cisSnpI]
  if (cisSnpI %% 100 == 0) cat("snp ",cisSnpI,"\n")
  hets=het[,locus]==1 & !is.na(data[,envName])
  heteq=as.numeric(t(genotype[cisSnp,hets]==1))
  if (sd(heteq)>0 & (min(table(heteq))>5)){
    x=env[hets]
    a=alt[hets,locus]
    r=ref[hets,locus]
    n=a+r
    a=pmin(a,r)
    num.samples=length(x)
    x.null=cbind(env=x,eqtl=heteq,intercept=numeric(num.samples)+1.0)
    x.full=cbind( x.null, interaction=heteq*x)
    if (length(unique(x))==2) if (any(table(data.frame(heteq,x))<5)) {
      cat("not enough samples\n")
      next
    }
    ei=eigen(t(x.full) %*% x.full) # used to check identifiable
    if (min(ei$values)/sum(ei$values) > 0.001){ 
      s=eagle.settings()
      s$normalised.depth=1
      s$max.iterations=5000
      
      s$convergence.tolerance=0.1
      s$traceEvery=100000
      s$learn.rev=F
      
      res=eagle.helper(list(a),list(n),list(x.full),list(x.null),s)
      pvals[cisSnpI]=res$p.values[1]
      cat(cisSnpI, "/", length(cisSnps), "\n")
    }
  }
}
