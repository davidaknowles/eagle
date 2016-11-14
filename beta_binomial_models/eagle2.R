require(rstan)

eagle2=stan_model("~/Dropbox/eagle/eagle/beta_binomial_models/eagle2.stan", auto_write = F, save_dso = F)

joint_linear_anova=function(ys, ns, gene_counts, library_size, xFull, xNull, concShape=1.001, concRate=0.001, algorithm="Newton",  ...) {
  stopifnot(all(xNull[[1]][,1]==1) & all(xFull[[1]][,1]==1))
  df=ncol(xFull[[1]])-ncol(xNull[[1]])
  
  N=dim(ys)[1]
  Ti=dim(ys)[2]
  K=dim(ys)[3]
  
  dat=list(N=N, K=K, P=ncol(xFull[[1]]), T=Ti, x_1=xFull[[1]], x_2=xFull[[2]], library_size=library_size, ys=ys, ns=ns, gene_counts=gene_counts, concShape=concShape, concRate=concRate)

  beta_init=as.array( 1000, numeric(dat$P-1) ) 
  init=list(beta=as.array(beta_init), conc=1)
  
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
