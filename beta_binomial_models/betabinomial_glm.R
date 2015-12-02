require(rstan)
# rstan_options(auto_write = TRUE)

source("utils.R")

BETABINOMIAL_GLM=stan_model(file="betabinomial_glm.stan", save_dso=F, auto_write=F)

# ys: numerator counts
# ns: denominator counts
# xFull: matrix of covariates for full model. First column must be ones. 
# xNull: matrix of covariates for null model. First column must be ones. 
# Prior on concentration parameter is Gamma(concShape,concRate)
betaBinomialGLM=function(ys,ns,xFull,xNull,concShape=1.0001,concRate=1e-4,...) {
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
  
  dat=list(N=length(ys),P=ncol(xNull),ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)
  
  # Fit null model
  fit_null <- optimizing(BETABINOMIAL_GLM, data=dat, init=init, algorithm="BFGS", hessian=T, as_vector=F)
  
  # Initialize alternative model using null model
  betaInit=numeric(ncol(xFull))
  betaInit[1:ncol(xNull)]=fit_null$par$beta
  initFull=list(conc=fit_null$par$conc, beta=betaInit)
  
  # Fit alternative model
  datFull=dat
  datFull$x=xFull
  datFull$P=ncol(xFull)
  fit_full <- optimizing(BETABINOMIAL_GLM, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F)
  
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
  ys=rbinom(nsamp,ns,logistic(.3+.3*rnorm(nsamp)))
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,runif(nsamp))
  anova_result=betaBinomialGLM(ys,ns,xFull,xNull)
  get_glm_coefs(anova_result$fit_null)
  get_glm_coefs(anova_result$fit_full)
}
