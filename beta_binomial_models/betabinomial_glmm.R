source("betabinomial_glm.R")
source("svem.R")

# compile stan models
BETABINOMIAL_GLMM=stan_model("betabinomial_glmm.stan", save_dso=F, auto_write=F)

betaBinomialGLMM=function(ys,ns,xFull,xNull,z,concShape=1.0001,concRate=1e-4,iterations=5000,elbo_samples=3000,...) {
  
  # check first ncol(xNull) columns of xFull and xNull match
  stopifnot(all(xNull==xFull[,1:ncol(xNull)]))
  
  # degress of freedom
  df=ncol(xFull)-ncol(xNull)
  
  # data for the null model
  dat=list(N=length(ys), P=ncol(xNull), K=max(z), ys=ys,ns=ns,x=xNull,z=z, concShape=concShape, concRate=concRate)
  
  # get stan_fit object so we can call sampler$grad_log_prob for the gradient
  sampler=get_sampler(BETABINOMIAL_GLMM, dat)
  
  # specific which parameters to optimize (rather than integrate over)
  sk=get_skeleton(sampler)
  sk$conc=T
  sk$logrev=T
  sk$beta[]=T
  sk$u[]=F
  to_optim=as.logical(unlist(sk))
  
  # we initialize using the model fit without the random effects
  o=optimizing(BETABINOMIAL_GLM, dat, as_vector=F)
  init=list(m=get_skeleton(sampler), s=get_skeleton(sampler))
  init$m$conc=log(o$par$conc)
  init$m$beta=o$par$beta
  for (n in names(init$s)) init$s[[n]][]=1 # set all standard deviations to 1 initially
  for (n in names(init)) init[[n]]=unlist(init[[n]])
  
  # set the seed so we can use the same seed for the alternative model
  set.seed(1)
  
  null_gradient_function=function(g) sampler$grad_log_prob(g,T)
  null_likelihood=function(g) sampler$log_prob(g,T,F)
  v_null=svem(null_gradient_function, to_optim, init, plot.elbo = F, iterations = iterations, samples_for_elbo = 0, log_prob = null_likelihood)
  
  # data for alternative (full) model
  dat_full=list(N=length(ys), P=ncol(xFull), K=max(z), ys=ys,ns=ns,x=xFull, z=z, concShape=concShape, concRate=concRate)
  
  # specify which parameters to optimize
  sampler2=get_sampler(BETABINOMIAL_GLMM, dat_full)
  sk_full=get_skeleton(sampler2)
  sk_full$conc=T
  sk_full$logrev=T
  sk_full$beta[]=T
  sk_full$u[]=F
  to_optim=as.logical(unlist(sk_full))
  
  # initialization choices: initializing from null fit seems best
  if (0) {
    o_init=o$par
    o_init$beta=c(o_init$beta, rep(0,df))
    o=optimizing(STANGLM_MV, dat_full, as_vector=F)
    init=list(m=get_skeleton(sampler2), s=get_skeleton(sampler2))
    init$m$conc=log(o$par$conc)
    init$m$beta=o$par$beta
    for (n in names(init$s)) init$s[[n]][]=1
  } else {
    init=list( m=rstan:::rstan_relist(v_null$m, sk), 
      s=rstan:::rstan_relist(v_null$s, sk) )
    init$m$beta=c(init$m$beta, rep(0,df))
    init$s$beta=c(init$s$beta, rep(1,df)) # optimized anyway so this value will be ignored
  }
  init=list(m=unlist(init$m),s=unlist(init$s))
  
  # fit alternative model
  set.seed(1)
  v_full=svem(function(g) sampler2$grad_log_prob(g,T), to_optim, init=init, plot.elbo = F, iterations = iterations, samples_for_elbo = 0, function(g) sampler2$log_prob(g,T,F))

  # using the same random draws to estimate the deviance is more efficient
  nintegrate=sum(!to_optim)
  loglr=mean(sapply(1:elbo_samples, function(g) {
    x=rnorm(nintegrate)
    v_full$elbo_func( x ) - v_null$elbo_func( x ) }))
  #loglr=mean(v_full$elbo_progress[(maxit/2):maxit] - v_null$elbo_progress[(maxit/2):maxit])
    
  list(loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ) )
}
