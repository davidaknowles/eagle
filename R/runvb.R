
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
       dependent.rev=F,
       allow.flips=F,
       learn.coeffs=T,
       trace=T)
}

run.all = function(alt,n,x,max.its=1000,tol=10.0,debug=F,flips=F,learn.rev=T,rev=1.0,trace=T,dependent.rev=F)
{
  s=default.settings()
  if (dependent.rev){
      s$normalised.depth=scale(unlist(lapply(n,sum)))
      s$rep.slope=.6
      s$rep.intercept=2.7
  }
  s$dependent.rev=dependent.rev
  s$max.iterations=max.its
  s$allow.flips=flips
  s$convergence.tolerance=tol
  s$debug=debug
  s$trace=trace
  s$learn.coeffs=F
  s$learn.rev=learn.rev
  s$random.effect.variance=rev
  res.null = run.vb(alt,n,x,s)
  s$learn.coeffs=T
  s$learn.rev=F
  s$random.effect.variance=res.null$random.effect.var
  s$rep.slope=res.null$rep.slope
  s$rep.intercept=res.null$rep.intercept
  res.full = run.vb(alt,n,x,s)
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
