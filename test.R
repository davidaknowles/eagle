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
.tmp = for (i in 1:n.loci){
    maf=runif(1)*.4+.1 # MAF between .1 and .5 for this locus
    # sample which individuals are heterozygous at this locus assuming HWE
    hap1=runif(n.samples)<maf 
    hap2=runif(n.samples)<maf
    hets=xor(hap1,hap2)
    numHets=sum(hets)
    ones=numeric(numHets)+1.0 
    x=environmentVar[hets]
    xFull[[i]]=cbind(x,ones) # design matrix at this locus for alternative ("full") model
    xNull[[i]]= cbind(ones) # matrix(0,ncol=0,nrow=numHets) # # design matrix at this locus for null model
    n[[i]]=rpois(numHets,100*rgamma(numHets,shape=2,rate=2)) # sample read depth from overdispersed Poisson (NB?)
    n[[i]][ n[[i]]<20 ]=20
    p=logistic(trueBeta[i]*x+.3*rnorm(numHets)) # sample underlying probabilities, with overdispersion
    alt[[i]]=rbinom(numHets,n[[i]],p) # sample alternative counts
    #for (j in 1:numHets) # "min" model: doesn't work well with the way I sample betas here
    #    alt[[i]][j]=min(alt[[i]][j],n[[i]][j]-alt[[i]][j])
}

#---------------- run the model --------------------------
s=default.settings()
s$debug=F
#s$rev.model=as.integer(3) # local regression
s$rev.model=as.integer(2)
s$normalised.depth=scale(log10(unlist(lapply(n,mean)))) # only required for rev.model=3
s$max.iterations=10000
s$convergence.tolerance=.001
s$coeff.regulariser=0.0
s$learn.rev=T

# initial hyperparamters
s$rep.global.shape=1.0
s$rep.global.rate=0.0033

s$traceEvery=1

system.time( res <- run.helper(alt,n,xFull,xNull,s) ) # 4s

# save(file="testNoRerun.RData",res,trueBeta)
cat("p-values for true hits:",res$p.values[trueBeta!=0],"\n")
cat("true hit betas: ",trueBeta[trueBeta!=0],"\n")

#hist( res$p.values[ trueBeta==0.0 ], main="p values for null sites")
