require(rstan)

joint_model_linear=stan_model("joint_model_linear.stan", auto_write = F, save_dso = F)
gene_ase_linear=stan_model("gene_ase_linear.stan", auto_write = F, save_dso = F)
#gene_ase_linear_no_overdispersion=stan_model("gene_ase_linear_no_overdispersion.stan")
neg_bin_model=stan_model("neg_bin_glm.stan", auto_write = F, save_dso = F)
#joint_model_linear_no_overdisp=stan_model("joint_model_linear_no_overdisp.stan")

robustSolve=function(si, eigenThreshold=0.01){
  svdSigma=eigen(si)
  svdSigma$values=pmax(svdSigma$values,eigenThreshold)
  svdSigma$vectors %*% diag(1/svdSigma$values) %*% t(svdSigma$vectors)
}

joint_linear_anova=function(ys, ns, gene_counts, library_size, xFull, xNull, concShape=1.001, concRate=0.001, algorithm="Newton",  ...) {
  stopifnot(all(xNull[[1]][,1]==1) & all(xFull[[1]][,1]==1))
  df=ncol(xFull[[1]])-ncol(xNull[[1]])
  
  dat=list(N=nrow(ys), K=ncol(ys), P=ncol(xFull[[1]]), x_1=xFull[[1]], x_2=xFull[[2]], library_size=library_size, ys=ys, ns=ns, gene_counts=gene_counts, concShape=concShape, concRate=concRate)
    
  model_choice= if (dat$K==0) neg_bin_model
      else if (all(library_size==0)) gene_ase_linear
     else joint_model_linear
  #no_overdispersion_model=if (all(library_size==0)) gene_ase_linear_no_overdispersion else no_overdispersion_model
  #beta_init_null=optimizing(no_overdispersion_model, dat, init=list(beta=as.array(beta_init)), as_vector=F, ... )$par$beta
  beta_init=as.array( c( if (all(library_size==0)) 1000 else mean(gene_counts/library_size)/2, numeric(dat$P-1) ) )
  init=list(beta=as.array(beta_init), nb_conc=100, conc=as.array(100+numeric(dat$K)))
  #init=optimizing(model_choice, dat=dat, init=init, as_vector=F, algorithm="LBFGS")$par
  #fit_null=optimizing(model_choice, dat=dat, init=init, as_vector=F, algorithm="Newton", ...)
  
  #dat$x_1=xFull[[1]]
  #dat$x_2=xFull[[2]]
  #dat$P=ncol(xFull[[1]])
  
  #beta_init_full=optimizing(no_overdispersion_model, dat, init=list(beta=as.array(c(beta_init_null,numeric(df)))), as_vector=F, ... )$par$beta
  #beta_init_full=NULL
  #init$beta=c(beta_init_null,numeric(df))
  #init=list(beta=beta_init_full, nb_conc=100, conc=100+numeric(dat$K))
  #init=optimizing(model_choice, dat=dat, init=init, as_vector=F,  algorithm="LBFGS")$par
  fit_full=optimizing(model_choice, dat=dat, init=init, as_vector=F, hessian=T, algorithm="Newton",  ...)
  
  dat$x_1=xNull[[1]]
  dat$x_2=xNull[[2]]
  dat$P=ncol(xNull[[1]])
  init_null=fit_full$par
  init_null$beta=as.array(init_null$beta[1:dat$P])
  fit_null=optimizing(model_choice, dat=dat, init=init_null, as_vector=F, algorithm="Newton",  ...)
  #if (refit_null$value > fit_null$value) fit_null=refit_null
  
  loglr=fit_full$value - fit_null$value
  refit_null_flag=F
  if (loglr>3) {
    init$beta=as.array(init$beta[1:dat$P])
    refit_null=optimizing(model_choice, dat=dat, init=init, as_vector=F, algorithm="Newton",  ...)
    if (refit_null$value > fit_null$value) {
      #cat("Using null model NOT initialized from full fit.\n")
      refit_null_flag=T
      fit_null=refit_null
      loglr=fit_full$value - fit_null$value
    }
  }
  
  list( loglr=loglr, df=df, lrtp=pchisq(2*loglr, df, lower.tail = F),  fit_null=fit_null, fit_full=fit_full,refit_null_flag=refit_null_flag )
}
