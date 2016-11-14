require(rstan)
# rstan_options(auto_write = TRUE)

#sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/betabinomial_glm_flips_rep.stan", save_dso=F, auto_write=F)
sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/bb_glm_flips_rep_prior_gene.stan", save_dso=F, auto_write=F)
#sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/betabinomial_glm_flips_rep_prior_0.stan", save_dso=F, auto_write=F)

# ys: numerator counts [n x T x K] where n are individuals, T are timepoints, K are SNPs
# ns: denominator counts [n x T x K]
# Prior on concentration parameter is Gamma(concShape,concRate)
betabinomial_glm_flips_rep_gene2=function(ys,ns,x,concShape=1.0001,concRate=1e-4,...) {
  
  N=dim(ys)[1]
  Ti=dim(ys)[2]
  K=dim(ys)[3]
  xNull=foreach(i=seq_len(Ti)) %do% matrix(1,N,1)
  dat=list(N=N,P=1,T=Ti,K=K,ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)
  
  # Fit null model
  fit_null <- optimizing(sm, data=dat,  algorithm="BFGS", hessian=T, as_vector=F)
  
  # Fit alternative model
  datFull=dat
  x=as.numeric(as.factor(x))
  stopifnot(length(x)==Ti)
  P=max(x)
  temp=matrix(0,N,P)
  temp[,1]=1
  
  # Initialize alternative model using null model
  initFull=fit_null$par
  initFull$beta=c(fit_null$par$beta,numeric(P-1))
  
  datFull$x=foreach(i=seq_len(Ti)) %do% { temp[,x[i]]=1; temp }
  datFull$P=P
  fit_full <- optimizing(sm, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F, verbose=F)
  
  loglr=fit_full$value - fit_null$value
  
  list( loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=P-1 ), fit_full=fit_full, fit_null=fit_null )
}
