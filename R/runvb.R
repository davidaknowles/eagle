
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
       learnBetas=T,
       rev.model="global",
       rep.global.rate=1.0,
       rep.global.shape=1.0,
       rep.slope=0.0,
       rep.intercept=0.0,
       rep.rep=1.0,
       flips.setting=0,
       flips.logodds.prior=-2.0,
       learn.flips.prior=T, 
	coeff.regulariser=0.0,
	init.flips=F,
       trace=T)
}

run.all = function(alt,n,xFull,xNull,max.its=1000,tol=10.0,debug=F,flips="none",learn.rev=T,rev=1.0,trace=T,rev.model="global",null.first=T,coeff.reg=0.0,fullToFlip=NA,nullToFlip=NA)
{
  s=default.settings()
  if (rev.model=="global"){
      s$rev.model=as.integer(0)
  } else if (rev.model=="regression" || rev.model=="local.regression") {
      s$rev.model=as.integer(1)
      s$rep.slope=.6
      s$rep.intercept=4.5
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
  if (!is.list(nullToFlip))
      nullToFlip=lapply(xNull,function(xHere) numeric(ncol(xHere)) + if (flips == "none") 0 else 1) 
  if (!is.list(fullToFlip))
      fullToFlip=lapply(xFull,function(xHere) numeric(ncol(xHere)) + if (flips == "none") 0 else 1) 
  
  s$coeff.regulariser=coeff.reg
  s$normalised.depth=scale(log10(unlist(lapply(n,mean))))
  s$max.iterations=max.its
  s$allow.flips=flips
  s$null.first=null.first
  s$convergence.tolerance=tol
  s$debug=debug
  s$trace=trace
  s$learn.flips.prior=T
  s$learn.rev=learn.rev
  s$random.effect.variance=rev
  s$toFlip=if (null.first) nullToFlip else fullToFlip
  # run first model -------------
  res.first = run.vb(alt,n,if (null.first) xNull else xFull,s)
  
  s$learn.flips.prior=F
  s$flips.logodds.prior=res.first$flips.logodds.prior
  s$learn.rev=F
  s$init.flips=T
  s$flips=res.first$flips
  s$random.effect.variance=res.first$random.effect.var
  s$rep.slope=res.first$rep.slope
  s$rep.rep=res.first$rep.rep
  s$rep.intercept=res.first$rep.intercept
  s$rep.global.rate=res.first$rep.global.rate
  s$rep.global.shape=res.first$rep.global.shape
  s$toFlip=if (null.first) fullToFlip else nullToFlip
  # run second model --------------
  res.second = run.vb(alt,n,if (null.first) xFull else xNull,s)
  if (null.first){
      res.full=res.second
      res.null=res.first
  } else {
      res.full=res.first
      res.null=res.second
  }
  log.like.ratios=2.0*(res.full$log.likelihoods-res.null$log.likelihoods)
  df=mapply(FUN=function(a,b) ncol(a)-ncol(b),xFull,xNull)
  p=1.0-pchisq(log.like.ratios,df=df)
  q=p.adjust(p,method="fdr")
  list(p.values=p,q.values=q,res.full=res.full,res.null=res.null)
}

run.perms = function(alt,n,xFull,xNull,max.its=1000,tol=10.0,debug=F,flips="none",learn.rev=T,rev=1.0,trace=T,rev.model="global",null.first=F,n.perms=10,coeff.reg=0.0,fullToFlip=NA,nullToFlip=NA)
{
  res=list()
  for (perm in 1:n.perms){
      cat("Running permutation number ",perm," of ",n.perms," ---------------\n")
      for (i in 1:length(xFull)) {
          to.permute=if (ncol(xFull[[i]])>=4) c(1,4) else 1
          xFull[[i]][,to.permute]=xFull[[i]][sample.int(nrow(xFull[[i]])),to.permute]
      }
      res[[perm]]=run.all(alt,n,xFull,xNull,max.its=max.its,tol=tol,debug=debug,flips=flips,learn.rev=learn.rev,rev=rev,trace=trace,rev.model=rev.model,null.first=null.first,coeff.reg=coeff.reg,nullToFlip=nullToFlip,fullToFlip=fullToFlip)
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
