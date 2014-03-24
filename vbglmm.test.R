# it: 99 re_var: 0.00119508 flips_prior_logodds 3.21112 lb: -1637.28
require(vbglmm)
set.seed(1)
n.loci=1
alt=list()
n=list()
xFull=list()
xNull=list()
fullToFlip=list()
nullToFlip=list()
for (i in 1:n.loci){
    fullToFlip[[i]]=c(1.0,1.0)
    nullToFlip[[i]]=c(1.0)
    n.samples=rpois(1,10)+2
    xFull[[i]]=cbind(rnorm(n.samples),numeric(n.samples)+1.0)
    xNull[[i]]=cbind(numeric(n.samples)+1.0)
    n[[i]]=rpois(n.samples,100)
    alt[[i]]=rbinom(n.samples,n[[i]],numeric(n.samples)+0.4)
}
s=default.settings()
flips="structured"
rev.model="global" 
  if (rev.model=="global"){
      s$rev.model=as.integer(0)
  } else if (rev.model=="regression" || rev.model=="local.regression") {
      s$rev.model=as.integer(1)
      s$rep.slope=.6
      s$rep.intercept=2.7
  } else if (rev.model=="local") {
      s$rev.model=as.integer(2)
  } else {
      error("Invalid random effect variance model: options are global, local, regression, local.regression")
  }
  if (rev.model=="local.regression") s$rev.model=as.integer(3)
  if (flips == "none"){
      s$flips.setting=0
  } else if (flips == "hard"){
      s$flips.setting=1
  } else if (flips == "soft"){
      s$flips.setting=2
  } else if (flips == "structured") {
      s$flips.setting=3
  } else {
      error("Invalid setting of flips: options are none, hard, soft")
  }
  fullToFlip=lapply(xFull,function(xHere) numeric(ncol(xHere)) + if (flips == "none") 0 else 1) 
  s$coeff.regulariser=0.0
  s$normalised.depth=scale(log10(unlist(lapply(n,sum))))
  s$max.iterations=100
  s$allow.flips=flips
  s$convergence.tolerance=0.0
  s$debug=T
  s$trace=T
  s$learn.flips.prior=T
  s$learn.rev=T
  s$learnBetas=T
  s$random.effect.variance=.1
  s$toFlip=fullToFlip
  res.first = run.vb(alt,n,xFull,s)


