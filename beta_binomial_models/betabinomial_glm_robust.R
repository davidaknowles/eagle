require(rstan)
# rstan_options(auto_write = TRUE)

source("utils.R")

BETABINOMIAL_GLM_ROBUST=stan_model(file="betabinomial_glm_robust.stan", save_dso=F, auto_write=F)

# ys: numerator counts
# ns: denominator counts
# xFull: matrix of covariates for full model. First column must be ones. 
# xNull: matrix of covariates for null model. First column must be ones. 
# Prior on concentration parameter is Gamma(concShape,concRate)
betaBinomialGLM_robust=function(ys,ns,xFull,xNull,concShape=1.0001,concRate=1e-4,...) {
  stopifnot(all(xNull==xFull[,1:ncol(xNull)]))
  stopifnot(all(xNull[,1]==1))
  # try to get sensible initialization
  rat=ys/ns # ratios
  # moment estimator of concentration parameter
  conc=pmin( 1/var(rat, na.rm=T), 1000 ) 
  # thresholded moment estimator of the mean
  m=pmin( pmax( mean(rat, na.rm=T), 1/1000 ), 1-1/1000)
  
  betaInit=numeric(ncol(xNull))
  betaInit[1]=log(m/(1.0-m))
  init=list(conc=conc, beta=as.array(betaInit))
  
  dat=list(N=length(ys),P=ncol(xNull),ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate, alpha=c(1.001,100,1.001), eps_a=1.001, eps_b=1000)
  
  # Fit null model
  fit_null <- optimizing(BETABINOMIAL_GLM_ROBUST, data=dat, init=init, algorithm="BFGS", hessian=T, as_vector=F)
  
  # Initialize alternative model using null model
  betaInit=numeric(ncol(xFull))
  betaInit[1:ncol(xNull)]=fit_null$par$beta
  initFull=list(conc=fit_null$par$conc, beta=betaInit)
  
  # Fit alternative model
  datFull=dat
  datFull$x=xFull
  datFull$P=ncol(xFull)
  fit_full <- optimizing(BETABINOMIAL_GLM_ROBUST, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F, verbose=F)
  
  loglr=fit_full$value - fit_null$value
  
  list( loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=ncol(xFull)-ncol(xNull) ), fit_full=fit_full, fit_null=fit_null )
}

# Function to extract coefficients, standard errors, and Wald p-values
get_glm_coefs=function(res) {
  p=length(res$par$beta)
  variance=robustSolve(-res$hessian)
  dimnames(variance)=dimnames(res$hessian)
  betase=sqrt(diag(variance))[paste0("beta.",1:p)]
  beta=res$par[paste0("beta[",1:p,"]")]
  zscore=res$par$beta / betase
  data.frame(co=res$par$beta,se=betase,p=2.0*pnorm(-abs(zscore)))
}

testbbglm=function() {
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  ns[1:10]=0
  ys=rbinom(nsamp,ns,logistic(.3+.3*rnorm(nsamp)))
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,runif(nsamp))
  #anova_result=betaBinomialGLM(ys,ns,xFull,xNull)
  #anova_robust=betaBinomialGLM_robust(ys,ns,xFull,xNull)
  
  ys[ sample.int(100,10) ]=0
  anova_result=betaBinomialGLM(ys,ns,xFull,xNull) 
  anova_robust=betaBinomialGLM_robust(ys,ns,xFull,xNull) 
  c(anova_result$lrtp,anova_robust$lrtp)
}
