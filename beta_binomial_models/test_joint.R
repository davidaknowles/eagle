
require(MASS)
require(doMC)
rbb=function(n,p,conc) rbinom(length(p), n, rbeta(length(p), p*conc, (1-p)*conc))

N=100
K=3
mu=20.0 

betas=c(0.0, 2.0, 5.0, 10)

results=foreach (beta_h=betas) %do% {
  print(beta_h)
  res=as.data.frame( do.call( rbind, foreach(i=1:200) %dopar% { 

    library_size=runif(N)
    test_snp=matrix(runif(N*2),N)<.5
    exonic_snp=array(runif(N*K*2),c(N,K,2))<.5
   
    allelic_expression=mu+test_snp * beta_h
    
    total_exp=rowSums(allelic_expression)
    
    gene_counts=rnegbin(N, total_exp * library_size, theta=10)
    genotype=rowSums(test_snp)
    #boxplot(gene_counts/library_size ~ genotype)
    
    hets=apply(exonic_snp,c(1,2),sum)==1
    
    ns=hets * 20
    ys=do.call(cbind,foreach(i=1:K) %do% rbb( ns[,i], allelic_expression[,1] / total_exp , 10))
    
    xFull=list( cbind(1,test_snp[,1]), cbind(1,test_snp[,2]) )
    xNull=list( matrix(1,N,1), matrix(1,N,1) )
    
    ys[sample.int(N,N/10)]=0
    
    c( joint_linear_anova(ys, ns, gene_counts, library_size, xFull, xNull)$lrtp, joint_linear_anova_robust(ys, ns, gene_counts, library_size, xFull, xNull)$lrtp )
  } ) )
  colnames(res)=c("Default","Robust")
  res
}

# both look pretty calibrated

names(results)=betas

pdf("varying_effect_joint_robust_10pc_off.pdf", height=8, width=10)
do.call(grid.arrange,foreach(beta=betas) %do% { multiqq(results[[as.character(beta)]]) + ggtitle(paste("True effect:",beta)) + scale_color_discrete(guide=guide_legend(title="Link\nfunction"))})
dev.off()

