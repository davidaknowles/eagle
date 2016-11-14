require(rstan)

joint_model_linear_robust=stan_model("joint_model_linear_robust.stan")
gene_ase_linear_robust=stan_model("gene_ase_linear_robust.stan")
neg_bin_model=stan_model("neg_bin_glm.stan")


joint_linear_anova_robust=function(ys, ns, gene_counts, library_size, xFull, xNull, concShape=1.001, concRate=0.001, algorithm="Newton",  ...) {
  stopifnot(all(xNull[[1]][,1]==1) & all(xFull[[1]][,1]==1))
  df=ncol(xFull[[1]])-ncol(xNull[[1]])
  
  dat=list(N=nrow(ys), K=ncol(ys), P=ncol(xNull[[1]]), x_1=xNull[[1]], x_2=xNull[[2]], library_size=library_size, ys=ys, ns=ns, gene_counts=gene_counts, concShape=concShape, concRate=concRate, alpha=c(1.001,100,1.001), eps_a=1.001, eps_b=1000)
    
  model_choice= if (dat$K==0) neg_bin_model
      else if (all(library_size==0)) gene_ase_linear_robust
     else joint_model_linear_robust

  beta_init_null=as.array( c( if (all(library_size==0)) 1000 else mean(gene_counts/library_size)/2, numeric(dat$P-1) ) )
  init=list(beta=as.array(beta_init_null), nb_conc=100, conc=as.array(100+numeric(dat$K)), eps=0.001, pi=foreach(i=seq_len(dat$K)) %do% { dat$alpha/sum(dat$alpha)})
  fit_null=optimizing(model_choice, dat=dat, init=init, as_vector=F, algorithm="Newton", ...)
  
  dat$x_1=xFull[[1]]
  dat$x_2=xFull[[2]]
  dat$P=ncol(xFull[[1]])

  beta_init_full=NULL
  init$beta=c(beta_init_null,numeric(df))
  fit_full=optimizing(model_choice, dat=dat, init=init, as_vector=F, hessian=T, algorithm="Newton",  ...)
  
  loglr=fit_full$value - fit_null$value
  
  list( loglr=loglr, df=df, lrtp=pchisq(2*loglr, df, lower.tail = F), beta_init_null=beta_init_null, fit_null=fit_null, beta_init_full=beta_init_full, fit_full=fit_full )
}
