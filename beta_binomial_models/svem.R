get_sampler=function(object, data) {
  rstan:::prep_call_sampler(object)
  model_cppname <- object@model_cpp$model_cppname 
  mod <- get("module", envir = object@dso@.CXXDSOMISC, inherits = FALSE) 
  stan_fit_cpp_module <- eval(call("$", mod, paste('stan_fit4', model_cppname, sep = '')))
  stopifnot(is.list(data))
  parsed_data <- rstan:::parse_data(get_cppcode(object))
  for (i in seq_along(data)) parsed_data[[names(data)[i]]] <- data[[i]]
  parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
  new(stan_fit_cpp_module, parsed_data, object@dso@.CXXDSOMISC$cxxfun)
}

get_skeleton=function(sampler) {
  m_pars <- sampler$param_names()
  idx_wo_lp <- which(m_pars != "lp__")
  m_pars <- m_pars[idx_wo_lp]
  p_dims <- sampler$param_dims()[idx_wo_lp]
  rstan:::create_skeleton(m_pars, p_dims)
}

adagrad=function(grad, x, master_stepsize=0.01, eps=1e-6, iterations=300, verbose=F) {
  historical_grad=0
  progress=list()
  for (i in 1:iterations) {
    g=grad(x)
    historical_grad=historical_grad+g^2
    x=x-master_stepsize*g/(eps+sqrt(historical_grad))
    progress[[i]]=attr(g,"log_prob")
    if (verbose) cat(i, attr(g,"log_prob"), "\n")
  }
  list(x=x,log_prob=unlist(progress))
}

svem=function(grad_log_prob, to_optim, init=NULL, plot.elbo=F, samples_for_elbo=10000, log_prob=NULL, ...) {

  #to_optim=c(F,T,T)
  noptim=sum(to_optim)
  nintegrate=sum(!to_optim)
  P=length(to_optim)
  
  elbo_grad_mixed=function(temp0, eta=rnorm(nintegrate)) {
    x=numeric(P)
    x[to_optim]=temp0[1:noptim]
    temp=temp0[(noptim+1):length(temp0)]
    m=temp[1:nintegrate]
    logs=temp[(nintegrate+1):(2*nintegrate)]
    s=exp(logs)
    x[!to_optim]=m+s*eta
    g=grad_log_prob(x)
    stopifnot(!any(is.na(g)))
    grad_m=g[!to_optim]
    grad_logs=grad_m*eta*s+1 # plus 1 from entropy
    res=-c(g[to_optim],grad_m,grad_logs)
    attr(res,"log_prob")=attr(g,"log_prob") + sum(logs)
    res
  }
  
  if (is.null(init))
    init=numeric(noptim+2*nintegrate) else
      if (is.list(init)) 
        init=c( init$m[to_optim], init$m[!to_optim], log(init$s[!to_optim]))
  
  stopifnot(length(init)==noptim+2*nintegrate)
  
  adagrad_fit=adagrad(elbo_grad_mixed, init, master_stepsize = 1, ...) 
  
  elbo_func=function(eta=rnorm(nintegrate)) attr( elbo_grad_mixed(adagrad_fit$x, eta), "log_prob" )
  
  elbo_estimate=if (samples_for_elbo>=1) mean(unlist(foreach(i=1:samples_for_elbo) %do% { elbo_func() }), na.rm=T) else NA
  
  if (plot.elbo) {
    require(stats)
    ma <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)}
    elbo_progression=adagrad_fit$log_prob
    elbo_ma=ma(elbo_progression,100)
    plot(elbo_progression,pch=16,col=rgb(.5,.5,.5,.5), ylim=c(min(elbo_progression[100:length(elbo_progression)]),max(elbo_progression)))
    lines(elbo_ma, col="red", lwd=2)
  }
  
  temp0=adagrad_fit$x
  
  if (!is.null(log_prob)) {
    x=numeric(P)
    x[to_optim]=temp0[1:noptim]
    temp=temp0[(noptim+1):length(temp0)]
    mh=temp[1:nintegrate]
    logs=temp[(nintegrate+1):(2*nintegrate)]
    sh=exp(logs)
    sumlogs=sum(logs)
    elbo_func=function(eta=rnorm(nintegrate)) {
      x[!to_optim]=mh+sh*eta
      log_prob(x) + sumlogs
    }
  }
  
  m=numeric(P)
  s=numeric(P)
  m[to_optim]=temp0[1:noptim]
  s[to_optim]=0
  temp=temp0[(noptim+1):length(temp0)]
  m[!to_optim]=temp[1:nintegrate]
  s[!to_optim]=exp(temp[(nintegrate+1):(2*nintegrate)])
  
  list(m=m, s=s, elbo=elbo_estimate, elbo_func= elbo_func, elbo_progress=adagrad_fit$log_prob)
}

log1pe=function(i) {
  lp=numeric(length(i))
  isLinear=i>10
  lp[isLinear]=i[isLinear]
  lp[!isLinear]=log1p(exp(i[!isLinear]))
  lp
}

log1pe_inv=function(x) log(exp(x)-1.0)

get_gamma_params=function(a,b,eta,eps=1e-10) {
  P=length(a)
  minz=pgamma(eps*(numeric(P)+1.0),shape=a,rate=b)
  z=minz+eta*(1.0-minz)
  
  x=numeric(P)
  dfda=numeric(P)
  
  useApprox=(a<1) & ((24-22.6*z)*log(a) < -10)
  zua=z[useApprox]
  aua=a[useApprox]
  bua=b[useApprox]
  
  logx=(log(zua)+log(aua)+lgamma(aua))/aua-log(bua)
  x[useApprox]=exp(logx)
  dlogxda=-(log(zua)+log(aua)+lgamma(aua))/aua^2+(1/aua+digamma(aua))/aua
  dfda[useApprox]=dlogxda * x[useApprox]
  
  mgi=mygaminv(z[!useApprox],a[!useApprox],b[!useApprox]) 
  x[!useApprox]=mgi$x
  dfda[!useApprox]=mgi$g
  
  dfdb=-x/b
  
  list(x=x, dfda=dfda, dfdb=dfdb)
}

gamma_svi=function(grad_log_prob, P, init=NULL, plot.elbo=F, iterations=1000, samples_for_elbo=10000) {

  elboHistory=list()
  
  # l=log(x)
  # glp=df/dl, dl/dx=1/x
  # f_x(x) = f_l(log(x)) dl/dx = f_l(log(x)) / x
  # log f_x(x) = log f_l(log(x)) - log(x)
  # d_x log f_x(x) = glp(log(x))/x - 1/x = (glp(l) - 1)/x
  
  elbo_grad=function(temp0, eta=runif(P)) {
    P=length(temp0)/2
    a=log1pe(temp0[1:P])
    b=log1pe(temp0[(P+1):(2*P)])
    
    gamma_params=get_gamma_params(a,b,eta)
    
    glp=(grad_log_prob(log(gamma_params$x))-1)/gamma_params$x
    logl=attr(glp,"log_prob")
    elbo=logl +sum( a - log(b) + lgamma(a) + (1-a)*digamma(a) )
    elboHistory[[ length(elboHistory)+1 ]]=elbo
    if (length(elboHistory)>10)
      elboHistory[[1]]=NULL 
    g_gamma=(g[modes==GAMMA_MODE] - 1.0) / gamma_params$x
    g_a = gamma_params$dfda * g_gamma + ( 1.0 + (1-a) * trigamma(a) ) 
    g_b = gamma_params$dfdb * g_gamma - 1.0/b
    g_a = g_a / (1+exp(-temp_a))
    g_b = g_b / (1+exp(-temp_b))  # account for transfomed variables
    ghat=c(g_a,g_b)
    attr(ghat,"log_prob")=elbo
    -ghat # adagrad does minimization usually, so flip the sign
  }
  
  if (is.null(init))
    init=numeric(P*2) else
      if (is.list(init)) 
        init=c( log1pe_inv(init$a), log1pe_inv(init$b) )
  
  stopifnot(length(init)==2*P)
  
  adagrad_fit=adagrad(elbo_grad, init, master_stepsize = 1, iterations = iterations) 
  
  elbo_func=function(eta=runif(P)) attr( elbo_grad(adagrad_fit$x, eta), "log_prob" )
  
  elbo_estimate=if (samples_for_elbo>=1) mean(unlist(foreach(i=1:samples_for_elbo) %do% { elbo_func() }), na.rm=T) else NA
  
  temp0=adagrad_fit$x
  
  a=log1pe(temp0[1:P])
  b=log1pe(temp0[(P+1):(2*P)])
  
  list(a=a, b=b, elbo=elbo_estimate, elbo_func=elbo_func, elbo_progress=adagrad_fit$log_prob)
}

# modes: 0=optimize, 1=gaussian, 2=gamma
OPTIMIZE_MODE=0
GAUSSIAN_MODE=1
GAMMA_MODE=2

svem_all=function(grad_log_prob, modes, init=NULL, plot.elbo=F, iterations=1000, samples_for_elbo=10000, log_prob=NULL) {
  
  noptim=sum(modes==OPTIMIZE_MODE)
  ngaussian=sum(modes==GAUSSIAN_MODE)
  ngamma=sum(modes==GAMMA_MODE)
  P=length(modes)
  
  elbo_grad_mixed=function(temp0, eta_normal=rnorm(ngaussian), eta_gamma=runif(ngamma)) {
    x=numeric(P)
    x[modes==OPTIMIZE_MODE]=temp0[1:noptim]
    
    m=temp0[noptim+seq_len(ngaussian)]
    logs=temp0[noptim+ngaussian+seq_len(ngaussian)]
    s=exp(logs)
    x[modes==GAUSSIAN_MODE]=m+s*eta_normal
    
    temp_a=temp0[noptim+ngaussian*2+seq_len(ngamma)]
    temp_b=temp0[noptim+ngaussian*2+ngamma+seq_len(ngamma)]
    a=log1pe(temp_a)
    b=log1pe(temp_b)
    gamma_params=get_gamma_params(a,b,eta_gamma)
    
    x[modes==GAMMA_MODE]=log(gamma_params$x)
    
    g=grad_log_prob(x)
    grad_m=g[modes==GAUSSIAN_MODE]
    grad_logs=grad_m*eta_normal*s+1.0 # plus from entropy
    
    g_gamma=(g[modes==GAMMA_MODE] - 1.0) / gamma_params$x
    g_a = gamma_params$dfda * g_gamma + ( 1.0 + (1-a) * trigamma(a) ) 
    g_b = gamma_params$dfdb * g_gamma - 1.0/b
    g_a = g_a / (1+exp(-temp_a))
    g_b = g_b / (1+exp(-temp_b))
    
    res=-c(g[modes==OPTIMIZE_MODE],grad_m,grad_logs,g_a,g_b)
    attr(res,"log_prob")=attr(g,"log_prob") + sum(logs) + # normal entropy
      sum( a - log(b) + lgamma(a) + (1-a)*digamma(a) ) + # gamma entropy
      - sum(x[modes==GAMMA_MODE]) # cancelling jacobian
    res
  }
  
  # m=a/b, s=a^.5/b, v=a/b^2. m/v=b, a=m*b
  
  if (is.null(init))
    init=numeric(noptim+2*ngaussian+2*ngamma) else
      if (is.list(init)) {
        b=init$m[modes==GAMMA_MODE]/(init$s[modes==GAMMA_MODE]^2)
        a=b*init$m[modes==GAMMA_MODE]
        init=c( init$m[modes==OPTIMIZE_MODE], init$m[modes==GAUSSIAN_MODE], log(init$s[modes==GAUSSIAN_MODE]), log1pe_inv(a), log1pe_inv(b) )
      }
  
  stopifnot(length(init)==noptim+2*ngaussian+2*ngamma)
  
  adagrad_fit=adagrad(elbo_grad_mixed, init, master_stepsize = 1, iterations = iterations) 
  
  elbo_func=function(eta_normal=rnorm(ngaussian),eta_gamma=runif(ngamma)) attr( elbo_grad_mixed(adagrad_fit$x, eta_normal,eta_gamma), "log_prob" )
  
  elbo_estimate=if (samples_for_elbo>=1) mean(unlist(foreach(i=1:samples_for_elbo) %do% { elbo_func() }), na.rm=T) else NA
  
  if (plot.elbo) {
    require(stats)
    ma <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)}
    elbo_progression=adagrad_fit$log_prob
    elbo_ma=ma(elbo_progression,100)
    plot(elbo_progression,pch=16,col=rgb(.5,.5,.5,.5), ylim=c(min(elbo_progression[100:length(elbo_progression)]),max(elbo_progression)))
    lines(elbo_ma, col="red", lwd=2)
  }
  
  temp0=adagrad_fit$x
  
  m=numeric(P)
  s=numeric(P)
  m[modes==OPTIMIZE_MODE]=temp0[seq_len(noptim)]
  s[modes==OPTIMIZE_MODE]=0
  
  temp=temp0[noptim+seq_len(ngaussian*2)]
  m[modes==GAUSSIAN_MODE]=temp[seq_len(ngaussian)]
  s[modes==GAUSSIAN_MODE]=exp(temp[ngaussian+seq_len(ngaussian)])
  
  a=log1pe(temp0[noptim+ngaussian*2+seq_len(ngamma)])
  b=log1pe(temp0[noptim+ngaussian*2+ngamma+seq_len(ngamma)])
  
  m[modes==GAMMA_MODE]=a/b
  s[modes==GAMMA_MODE]=sqrt(a)/b 
  
  list(m=m, s=s, elbo=elbo_estimate, elbo_func= elbo_func, elbo_progress=adagrad_fit$log_prob)
}
