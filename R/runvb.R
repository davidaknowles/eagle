
run.vb = function(alt,n,x,settings)
{
  .Call("runvb", alt, n, x, settings, PACKAGE = "vbglmm" )
}

fit.null = function(alt,n,max.its=1000,tol=10.0,random.effect.var=1.0,debug=F)
{
  .Call("fitnull", alt, n, max.its, tol, random.effect.var, debug, PACKAGE = "vbglmm" )
}

default.settings = function(){
  list(debug=F,
       max.iterations=1000,
       convergence.tolerance=10.0,
       random.effect.variance=1.0,
       learn.rev=T,
       rev.model="global",
       rep.global.rate=1.0,
       rep.global.shape=1.0,
       rep.slope=0.0,
       rep.intercept=0.0,
       rep.rep=1.0,
       flips.setting=0,
       flips.logodds.prior=-2.0,
       learn.flips.prior=T, 
       learn.coeffs=T,
       trace=T)
}

run.all = function(alt,n,x,max.its=1000,tol=10.0,debug=F,flips="none",learn.rev=T,rev=1.0,trace=T,rev.model="global",null.first=F)
{
  s=default.settings()
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
  } else {
      error("Invalid setting of flips: options are none, hard, soft")
  }
      
  s$normalised.depth=scale(log10(unlist(lapply(n,sum))))
  s$max.iterations=max.its
  s$allow.flips=flips
  s$null.first=null.first
  s$convergence.tolerance=tol
  s$debug=debug
  s$trace=trace
  s$learn.coeffs=!null.first
  s$learn.flips.prior=s$learn.coeffs
  s$learn.rev=learn.rev
  s$random.effect.variance=rev
  res.first = run.vb(alt,n,x,s)
  s$learn.coeffs=null.first
  s$learn.flips.prior=null.first
  s$flips.logodds.prior=res.first$flips.logodds.prior
  s$learn.rev=F
  s$random.effect.variance=res.first$random.effect.var
  s$rep.slope=res.first$rep.slope
  s$rep.rep=res.first$rep.rep
  s$rep.intercept=res.first$rep.intercept
  s$rep.global.rate=res.first$rep.global.rate
  s$rep.global.shape=res.first$rep.global.shape
  res.second = run.vb(alt,n,x,s)
  if (null.first){
      res.full=res.second
      res.null=res.first
  } else {
      res.full=res.first
      res.null=res.second
  }
  log.like.ratios=2.0*(res.full$log.likelihoods-res.null$log.likelihoods)
  p=1.0-pchisq(log.like.ratios,df=1)
  q=p.adjust(p,method="fdr")
  list(p.values=p,q.values=q,res.full=res.full,res.null=res.null)
}


run.perms = function(alt,n,x,max.its=1000,tol=10.0,debug=F,flips=F,n.perms=10,learn.rev=T,rev=1.0)
{
  res=list()
  for (perm in 1:n.perms){
      cat("Running permutation number ",perm," of ",n.perms," ---------------\n")
      for (i in 1:length(x)) x[[i]]=sample(x[[i]])
      res[[perm]]=run.all(alt,n,x,max.its=max.its,tol=tol,debug=debug,flips=flips,learn.rev=learn.rev,rev=rev)
  }
  res
}

## run.perms = function(alt,n,x,max.its=1000,tol=10.0,debug=F,nperm=1000,fdr=0.05)
## {
##   res = run.vb(alt,n,x,max.its=max.its,tol=tol,debug=debug)
##   res.null = fit.null(alt,n,max.its=max.its,tol=tol,random.effect.var=res$random.effect.var,debug=debug)
##   log.like.ratios=2.0*(res$log.likelihoods-res.null$log.likelihoods)
##   p=1.0-pchisq(log.like.ratios,df=1)
##   q=p.adjust(p,method="fdr")
##   hits=q<0.05
##   alt.hits=alt[hits]
##   n.hits=n[hits]
##   x.hits=x[hits]
##   res.perm=list()
##   for (i in 1:nperm){
##       cat("perm ",i,"\n")
##       for (j in 1:length(x.hits))
##           x.hits[[j]]=sample(x.hits[[j]])
##       res.perm[[i]] = run.vb(alt.hits,n.hits,x.hits,rev=res$random.effect.var,max.its=max.its,tol=tol,debug=debug)
##   }
##   list(p.values=p,q.values=q,res.full=res,res.null=res.null,hits=hits,res.perm=res.perm)
## }
