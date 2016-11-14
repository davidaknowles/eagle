# This is a test script for the EAGLE, the Environment-Ase Generalized Linear modEl
# Data sampled from the generative model and then we run inference to attempt to
# recover the true "hits"

require(eagle) # load the package
set.seed(1) # for reproducibility

#-------------- generate some synthetic data -------------------------
n.loci=100 # number of loci
n.samples=200 # number of individuals
alt=list() # list of alternative counts
n=list() # list of total counts
xFull=list() # list of design matrices for the alternative models
xNull=list() # list of design matrices for the null models
environmentVar=rnorm(n.samples) # values of the environment variable, e.g. age
# make "true" regression coefficients, most of which are 0
trueBeta=ifelse(runif(n.loci) < 0.05, rgamma(n.loci,shape=2,rate=1), 0.0) 
logistic=function(x) 1/(1+exp(-x))
for (i in 1:n.loci){
    maf=runif(1)*.4+.1 # MAF between .1 and .5 for this locus
    # sample which individuals are heterozygous at this locus assuming HWE
    hap1=runif(n.samples)<maf 
    hap2=runif(n.samples)<maf
    hets=xor(hap1,hap2)
    numHets=sum(hets)
    ones=numeric(numHets)+1.0 
    x=environmentVar[hets]
    xFull[[i]]=cbind(x,ones) # design matrix at this locus for alternative ("full") model
    xNull[[i]]=cbind(ones) # # design matrix at this locus for null model
    n[[i]]=rpois(numHets,100*rgamma(numHets,shape=2,rate=2)) # sample read depth from overdispersed Poisson (NB?)
    n[[i]][ n[[i]]<20 ]=20
    p=logistic(trueBeta[i]*x+.3*rnorm(numHets)) # sample underlying probabilities, with overdispersion
    alt[[i]]=rbinom(numHets,n[[i]],p) # sample alternative counts
    # "min" model: doesn't work well with the way I sample betas (and having no intercept)
    # alt[[i]]=pmin(alt[[i]],n[[i]]-alt[[i]])
}

#---------------- run the model --------------------------
s=eagle.settings()
s$debug=F
#s$rev.model=as.integer(3) # local regression
s$rev.model=2
s$normalised.depth=scale(log10(unlist(lapply(n,mean)))) # only required for rev.model=3
s$max.iterations=10000
s$convergence.tolerance=.001
s$coeff.regulariser=0.0
s$learn.rev=T

# initial hyperparamters
s$rep.global.shape=1.0
s$rep.global.rate=0.0033

s$traceEvery=1

system.time( res <- eagle.helper(alt,n,xFull,xNull,s) ) # 4s

# save(file="testNoRerun.RData",res,trueBeta)
cat("p-values for true hits:",res$p.values[trueBeta!=0],"\n")
cat("true hit betas: ",trueBeta[trueBeta!=0],"\n")

#hist( res$p.values[ trueBeta==0.0 ], main="p values for null sites")

s$learn.rev=F

totest=which(trueBeta==0)[1]

scales=c(1,2,5,10,20,50,100,200,500,1000,2000,5000, 1e4)

ps=unlist(foreach(scal=scales) %do% {
  xTemp=xFull[[totest]]
  xTemp[,1]=xTemp[,1]*scal
  eagle.helper(list( alt[[totest]] ), list( n[[totest]] ), list( xTemp ), list( xNull[[totest]]),s)$p.values[1]
})
qplot(scales,ps) + scale_x_log10()+  scale_y_log10( limits=c(1e-14,1e-13))

intercepts=scales
ps=unlist(foreach(inter=intercepts) %do% {
  xTemp=xFull[[totest]]
  xTemp[,1]=xTemp[,1]+ inter
  eagle.helper(list( alt[[totest]] ), list( n[[totest]] ), list( xTemp ), list( xNull[[totest]]),s)$p.values[1]
})
qplot(intercepts,ps) + scale_x_log10()+ scale_y_log10( limits=c(1e-14,1e-13))


xTemp=xFull[[totest]]
xTemp[,1]=xTemp[,1]*0.00001 + 10
eagle.helper(list( alt[[totest]] ), list( n[[totest]] ), list( xTemp ), list( xNull[[totest]]),s)$p.values[1]
qplot(xFull[[totest]][,1], alt[[totest]]/n[[totest]] ) + ylim(0,.5)


unlist(foreach(i=1:100) %dopar% {
  maf=runif(1)*.4+.1 # MAF between .1 and .5 for this locus
  # sample which individuals are heterozygous at this locus assuming HWE
  hap1=runif(n.samples)<maf 
  hap2=runif(n.samples)<maf
  hets=xor(hap1,hap2)
  numHets=sum(hets)
  ones=numeric(numHets)+1.0 
  x=environmentVar[hets]
  xFull=cbind(x,ones) # design matrix at this locus for alternative ("full") model
  xNull=cbind(ones) # # design matrix at this locus for null model
  n=rpois(numHets,100*rgamma(numHets,shape=2,rate=2)) # sample read depth from overdispersed Poisson (NB?)
  n[ n<20 ]=20
  p=logistic(5+.3*rnorm(numHets)) # sample underlying probabilities, with overdispersion
  alt=rbinom(numHets,n,p) 
  alt=pmin(alt,n-alt)
  qplot(xFull[,1], alt/n ) + ylim(0,.5)
  
  xFull[,1]=xFull[,1]*0.00001 + 10
  eagle.helper(list( alt ), list( n ), list( xFull ), list( xNull),s)$p.values[1]
}) -> p