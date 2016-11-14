require(rstan)
# rstan_options(auto_write = TRUE)

#sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/betabinomial_glm_flips_rep.stan", save_dso=F, auto_write=F)
sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/bb_glm_flips_rep_prior.stan", save_dso=F, auto_write=F)
#sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/bb_glm_flips_rep_prior_0.stan", save_dso=F, auto_write=F)

# ys: numerator counts [nx2]
# ns: denominator counts [nx2]
# Prior on concentration parameter is Gamma(concShape,concRate)
betabinomial_glm_flips_rep=function(ys,ns,concShape=1.0001,concRate=1e-4,...) {
  
  N=nrow(ys)
  xNull=matrix(1,N,1)
  xNull=list(xNull,xNull)
  dat=list(N=N,P=1,ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)
  
  # Fit null model
  fit_null <- optimizing(sm, data=dat,  algorithm="BFGS", hessian=T, as_vector=F)
  
  # Initialize alternative model using null model
  initFull=fit_null$par
  initFull$beta=c(fit_null$par$beta,0)
  
  # Fit alternative model
  datFull=dat
  datFull$x=list(cbind(xNull[[1]],0), cbind(xNull[[1]],1))
  datFull$P=2
  fit_full <- optimizing(sm, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F, verbose=F)
  
  loglr=fit_full$value - fit_null$value
  
  list( loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=1 ), fit_full=fit_full, fit_null=fit_null )
}

test_betabinomial_glm_flips_rep=function() {
  
  nsamp=50
  ns=matrix( rpois(nsamp*2,lambda=10) , nsamp, 2)
  temp=matrix(0,nsamp,2)
  temp[,2]=1
  phi=
  ys=matrix( rbinom(nsamp*2, as.numeric(ns),logistic( .3+as.numeric(temp)+.3*rnorm(nsamp))), nsamp, 2)
  m=melt(ys/ns)
  anova_result=betabinomial_glm_flips_rep(ys,ns) 
  anova_robust=betaBinomialGLM_robust(ys,ns,xFull,xNull) 
  c(anova_result$lrtp,anova_robust$lrtp)
  
}
