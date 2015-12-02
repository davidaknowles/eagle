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

adagrad=function(grad, x, master_stepsize=0.01, eps=1e-6, iterations=300) {
  historical_grad=0
  progress=list()
  for (i in 1:iterations) {
    g=grad(x)
    historical_grad=historical_grad+g^2
    x=x-master_stepsize*g/(eps+sqrt(historical_grad))
    progress[[i]]=attr(g,"log_prob")
    #print(x)
  }
  list(x=x,log_prob=unlist(progress))
}

svem=function(grad_log_prob, to_optim, init=NULL, plot.elbo=F, iterations=1000, samples_for_elbo=10000, log_prob=NULL) {

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
    #print(x)
    g=grad_log_prob(x)
    grad_m=g[!to_optim]
    grad_logs=grad_m*eta*s+1 # plus from entropy
    res=-c(g[to_optim],grad_m,grad_logs)
    attr(res,"log_prob")=attr(g,"log_prob") + sum(logs)
    res
  }
  
  if (is.null(init))
    init=numeric(noptim+2*nintegrate) else
      if (is.list(init)) 
        init=c( init$m[to_optim], init$m[!to_optim], log(init$s[!to_optim]))
  
  stopifnot(length(init)==noptim+2*nintegrate)
  
  adagrad_fit=adagrad(elbo_grad_mixed, init, master_stepsize = 1, iterations = iterations) 
  
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


