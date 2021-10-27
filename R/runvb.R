packageName="eagle"


# alt: list (over exonic SNPs) of alternative read counts
# n: list (over exonics SNPs) of total read counts
# x: list of design matrices for each exonic SNP
# settings: list of settings. Create using default.settings() and then customize. 
eagle.vem <- function(alt,n,x,settings)
{
    .Call("runvb", alt, n, x, settings, PACKAGE = packageName )
}

eagle.settings <- function(){
  list(debug=F, # output debugging information? Also performs additional checks. 
       max.iterations=1000, # maximum iterations of EM to run
       convergence.tolerance=0.1, # consider EM to have converged if the change in the lower bound is less than this
       random.effect.variance=1.0, # initial random effect variance
       learn.rev=F, # whether to learn the overdispersion hyperparameters a,b
       rev.model="local", # rev.model: one of c("global","regression","local","local.regression"). global: single fixed random effect variance across all exonic SNPs (not recommended). regression: random effect variance is log-linear in the log average read depth (not recommended). local: random effect variance at exonic SNP s is v_s ~ InverseGamma(a,b) [recommended]. local.regression: log(rev)  ~ N(mt+c,v) [gives similar results to "local" but not as well calibrated]
       rep.global.shape=1.1, # the shape parameter (a in the manuscript) in the default "local" overdispersion prior model
       rep.global.rate=0.0033, # the rate parameter (b in the manuscript) in the default "local" overdispersion prior model
       rep.slope=0.0, # m in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)
       rep.intercept=0.0, # c in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)
       rep.rep=1.0, # 1/v in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)
       coeff.regulariser=0.0,#  whether to regularize the regression coefficents with a term -coef.reg |beta|^2
       return.aux.variables=F, # whether to return the auxiliary variables g (only used for debugging)
       storeAllCoeffs=F, # whether to store coefficients throughout the EM algorithm (for debugging)
       null.first=T, # whether to run the null or alternative model first. Only matters if learning the overdispersion hyperparamters, in which case these are learnt on whichever is run first, and fixed for the second. 
       rerunFirst=F, # whether to rerun the first model after learning hyperparameters (not used)
       traceEvery=1, # traceEvery: how often to output convergence info to stdout
       learnRepRep=T, # whether to learn v in the local.regression model
       learnBetas=T # whether to learn the regression coefficients)
  )
}

# alt: list (over exonic SNPs) of alternative read counts
# n: list (over exonics SNPs) of total read counts
# xFull: list of design matrices for the alternative hypothesis (e.g. including environment)
# xNull: list of design matrices for the null hypothesis
# max.its: maximum iterations of EM to run
# tol: consider EM to have converged if the change in the lower bound is less than this
# debug: output debugging information? Also performs additional checks. 
# learn.rev: whether to learn the overdispersion hyperparameters (e.g. a,b) or hold them fixed
# rev: random effect variance
# traceEvery: how often to output convergence info to stdout
# rev.model: one of c("global","regression","local","local.regression"). global: single fixed random effect variance across all exonic SNPs (not recommended). regression: random effect variance is log-linear in the log average read depth (not recommended). local: random effect variance at exonic SNP s is v_s ~ InverseGamma(a,b) [recommended]. local.regression: log(rev)  ~ N(mt+c,v) [gives similar results to "local" but not as well calibrated]
# null.first: whether to run the null or alternative model first. Only matters if learning the overdispersion hyperparamters, in which case these are learnt on whichever is run first, and fixed for the second. 
# coeff.reg: whether to regularize the regression coefficents with a term -coef.reg |beta|^2
# return.aux.variables: whether to return the auxiliary variables g (only used for debugging)
# storeAllCoeffs: whether to store coefficients throughout the EM algorithm (for debugging)
# repRep: for local.regression the inital value of 1/v
eagle <- function(alt,n,xFull,xNull,max.its=1000,tol=10.0,debug=F,learn.rev=T,rev=1.0,traceEvery=1,rev.model="global",null.first=T,coeff.reg=0.0,return.aux.variables=F,storeAllCoeffs=F,repRep=1){
  s=eagle.settings()
  s$return.aux.variables=return.aux.variables
  s$storeAllCoeffs=storeAllCoeffs
  if (rev.model=="global"){
      s$rev.model=as.integer(0)
  } else if (rev.model=="regression" || rev.model=="local.regression") {
      s$rev.model=as.integer(1)
      s$rep.slope=.6
      s$rep.intercept=4.5
      s$rep.rep=repRep # note: not used for "regression", equivalent to repRep=Inf
  } else if (rev.model=="local") {
      s$rev.model=as.integer(2)
  } else {
      error("Invalid random effect variance model: options are global, local, regression, local.regression")
  }
  if (rev.model=="local.regression") s$rev.model=as.integer(3)
  s$coeff.regulariser=coeff.reg
  s$normalised.depth=scale(log10(unlist(lapply(n,mean))))
  s$max.iterations=max.its
  s$null.first=null.first
  s$convergence.tolerance=tol
  s$debug=debug
  s$traceEvery=traceEvery
  s$learn.rev=learn.rev
  s$random.effect.variance=rev

  eagle.helper(alt,n,xFull,xNull,s)
}

# alt: list (over exonic SNPs) of alternative read counts
# n: list (over exonics SNPs) of total read counts
# xFull: list of design matrices for the alternative hypothesis (e.g. including environment)
# xNull: list of design matrices for the null hypothesis
# s: list of settings. Create using default.settings() and then customize. 
eagle.helper = function(alt,n,xFull,xNull,s){

    if (is.null(s$normalised.depth))
        s$normalised.depth=scale(log10(unlist(lapply(n,mean))))

  stopifnot( !is.nan(s$normalised.depth))

  stopifnot( is.list(alt), is.list(n), is.list(xFull), is.list(xNull) )
    
  if (any(sapply(alt,function(g) any(is.na(g))))) stop("NAs not allowed in alt, please remove these elements.")
  if (any(sapply(n,function(g) any(is.na(g))))) stop("NAs not allowed in n, please remove these elements.")
  if (any(sapply(xFull,function(g) any(is.na(g))))) stop("NAs not allowed in xFull, please remove these elements.")
  if (any(sapply(xNull,function(g) any(is.na(g))))) stop("NAs not allowed in xNull, please remove these elements.")

  stopifnot( length(alt)==length(n), length(n)==length(xFull), length(xFull)==length(xNull))
  stopifnot( sapply(alt,length)==sapply(n,length), sapply(n,length)==sapply(xFull,nrow), sapply(xFull,nrow)==sapply(xNull,nrow) )

  stopifnot( all(sapply(alt,class)  %in% c("integer","numeric") ), all(sapply(n,class) %in% c("integer","numeric") ) )
#  stopifnot( all(sapply(xFull,class)=="matrix" ), all(sapply(xNull,class)=="matrix") )

  if (any( sapply(xFull,ncol) <= sapply(xNull,ncol) )) warning("Some xFull matrices have less than or equal columns compared to xNull") 
  
  stopifnot( all(sapply(alt,mode)=="numeric") )
  stopifnot( all(sapply(n,mode)=="numeric") )
  stopifnot( all(sapply(xFull,mode)=="numeric") )
  stopifnot( all(sapply(xNull,mode)=="numeric") ) 

  problem_loci=which( sapply(xNull, function(g) det( t(g) %*% g )==0.0) )
  if (length(problem_loci)>0 & s$coeff.regulariser==0.0) stop( paste( "Det(xNull)=0 at", paste(problem_loci, collapse=" ") ) )
  if (length(problem_loci)>0 & s$coeff.regulariser>0.0) warning( paste( "Det(xNull)=0 at", paste(problem_loci, collapse=" ")  ) )    

  problem_loci=which( sapply(xFull, function(g) det( t(g) %*% g )==0.0) )
  if (length(problem_loci)>0 & s$coeff.regulariser==0.0) stop( paste( "Det(xFull)=0 at", paste(problem_loci, collapse=" ") ) )
  if (length(problem_loci)>0 & s$coeff.regulariser>0.0) warning( paste( "Det(xFull)=0 at", paste(problem_loci, collapse=" ")  ) )    

  xToListList = function(x) lapply( x, function(g) lapply( as.list(1:nrow(g)), function(h) g[h,] ))
    
  xFullList=xToListList(xFull)
  xNullList=xToListList(xNull)
    
  # run first model -------------
  timeFirst = system.time( res.first <- eagle.vem(alt,n,if (s$null.first) xNullList else xFullList,s) )[1]

  # hold global dispersion parameters fixed
  s$learn.rev=F
  s$random.effect.variance=res.first$random.effect.var
  s$rep.slope=res.first$rep.slope
  s$rep.rep=res.first$rep.rep
  s$rep.intercept=res.first$rep.intercept
  s$rep.global.rate=res.first$rep.global.rate
  s$rep.global.shape=res.first$rep.global.shape
  
  # rerun first model with fixed rev model? (doesn't make much difference)
  res.firstOld=NULL
  if (s$rerunFirst){
      res.firstOld=res.first
      res.first = eagle.vem(alt,n,if (s$null.first) xNullList else xFullList,s)
  }
  
  # run second model --------------
  timeSecond = system.time( res.second <- eagle.vem(alt,n,if (s$null.first) xFullList else xNullList,s) )[1]
  if (s$null.first){
      res.full=res.second
      res.null=res.first
  } else {
      res.full=res.first
      res.null=res.second
  }
  
  like.ratios.statistics=2.0*(res.full$log.likelihoods-res.null$log.likelihoods)
  # degress of freedom calculation
  df=mapply(FUN=function(a,b) ncol(a)-ncol(b),xFull,xNull)
  # calculate p-values
  p=pchisq(like.ratios.statistics,df=df,lower.tail=F)
  # calculate q-values
  q=p.adjust(p,method="fdr")
  list(p.values=p,q.values=q,res.full=res.full,res.null=res.null,settings=s,timeFirst=timeFirst,timeSecond=timeSecond,res.firstOld=res.firstOld)
}


